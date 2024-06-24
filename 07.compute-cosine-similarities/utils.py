import pandas as pd
import numpy as np
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_array, as_float_array
from sklearn.base import TransformerMixin, BaseEstimator
import scipy
from sklearn.metrics import average_precision_score
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import random
import itertools
import numpy as np
from typing import List, Tuple
import re
from copairs import compute
from copairs.matching import Matcher, UnpairedException
import logging
from scipy.sparse import coo_matrix

logger = logging.getLogger("copairs")


def get_metacols(df):
    """return a list of metadata columns"""
    return [c for c in df.columns if c.startswith("Metadata_")]


def get_featurecols(df):
    """returna  list of featuredata columns"""
    return [c for c in df.columns if not c.startswith("Metadata")]


def get_metadata(df):
    """return dataframe of just metadata columns"""
    return df[get_metacols(df)]


def get_featuredata(df):
    """return dataframe of just featuredata columns"""
    return df[get_featurecols(df)]


def remove_negcon_empty_wells(df):
    """return dataframe of non-negative control wells"""
    df = (
        df.query('Metadata_pert_type!="negcon"')
        .dropna(subset=["Metadata_Symbol"])
        .reset_index(drop=True)
    )
    return df


def remove_low_infection_efficiency_wells(df):
    """return dataframe of non-low infection efficiency wells"""

    reagents_to_remove_df = pd.read_csv(
        "../00.download-and-process-annotations/output/orf-reagents-low-infection-efficiency-and-outliers.csv.gz",
        usecols=["Metadata_plate_map_name", "Metadata_Well", "Metadata_broad_sample"],
    )

    df = (
        df.merge(
            reagents_to_remove_df,
            on=["Metadata_plate_map_name", "Metadata_Well", "Metadata_broad_sample"],
            how="left",
            indicator=True,
        )
        .query('_merge=="left_only"')
        .drop(columns=["_merge"])
    )

    return df


def remove_empty_wells(df):
    """return dataframe of non-empty wells"""
    df = df.dropna(subset=["Metadata_Symbol"]).reset_index(drop=True)
    return df


def concat_profiles(df1, df2):
    """Concatenate dataframes"""
    if df1.shape[0] == 0:
        df1 = df2.copy()
    else:
        frames = [df1, df2]
        df1 = pd.concat(frames, ignore_index=True, join="inner")

    return df1


def consensus(profiles_df, group_by_feature):
    """
    Computes the median consensus profiles.
    Parameters:
    -----------
    profiles_df: pandas.DataFrame
        dataframe of profiles
    group_by_feature: str
        Name of the column
    Returns:
    -------
    pandas.DataFrame of the same shape as `plate`
    """

    metadata_df = get_metadata(profiles_df).drop_duplicates(subset=[group_by_feature])

    feature_cols = [group_by_feature] + get_featurecols(profiles_df)
    profiles_df = (
        profiles_df[feature_cols].groupby([group_by_feature]).median().reset_index()
    )

    profiles_df = metadata_df.merge(profiles_df, on=group_by_feature)

    return profiles_df


def remove_nan_features(df):
    """
    Remove NaN features.
    Parameters
    ----------
    df: pandas.DataFrame
        DataFrame of profiles.
    Returns
    -------
    df: pandas.DataFrame
        DataFrame of profiles after dropping NaN features.
    """

    r, c = np.where(df.isna())
    features_to_remove = [
        _ for _ in list(df.columns[list(set(c))]) if not _.startswith("Metadata_")
    ]
    print(f"Removed nan features: {features_to_remove}")
    df.drop(features_to_remove, axis=1, inplace=True)
    return df


def flatten_str_list(*args):
    """create a single list with all the params given"""
    columns = set()
    for col in args:
        if isinstance(col, str):
            columns.add(col)
        elif isinstance(col, dict):
            columns.update(itertools.chain.from_iterable(col.values()))
        else:
            columns.update(col)
    columns = list(columns)
    return columns


def evaluate_and_filter(df, columns) -> Tuple[pd.DataFrame, List[str]]:
    """Evaluate the query and filter the dataframe"""
    parsed_cols = []
    for col in columns:
        if col in df.columns:
            parsed_cols.append(col)
            continue

        column_names = re.findall(r"(\w+)\s*[=<>!]+", col)
        valid_column_names = [col for col in column_names if col in df.columns]
        if not valid_column_names:
            raise ValueError(f"Invalid query or column name: {col}")

        try:
            df = df.query(col)
            parsed_cols.extend(valid_column_names)
        except:
            raise ValueError(f"Invalid query expression: {col}")

    return df, parsed_cols


def cosine_similarity(id, feats, batch_size=20000):
    n_id = len(id)
    indices = [_ for _ in range(n_id)]

    pos_pairs = [list(itertools.combinations(indices, 2))]
    pos_total = sum(len(p) for p in pos_pairs)

    if pos_total == 0:
        raise UnpairedException("Unable to find positive pairs.")
    pos_pairs = np.fromiter(
        itertools.chain.from_iterable(pos_pairs),
        dtype=np.dtype((np.int32, 2)),
        count=pos_total,
    )

    logger.info("Computing positive similarities...")
    pos_sims = compute.pairwise_cosine(feats, pos_pairs, batch_size)

    rows, cols = zip(*pos_pairs)
    cos_sim = coo_matrix((pos_sims, (rows, cols)), shape=(n_id, n_id)).toarray()

    cos_sim += cos_sim.T
    np.fill_diagonal(cos_sim, 1)

    cos_sim_df = pd.DataFrame(
        cos_sim,
        index=id,
        columns=id,
    )

    return cos_sim_df


# Copied from https://github.com/cytomining/pycytominer/blob/90c9a899cf634d96da3a325fb07c4cfdb9640c58/pycytominer/cyto_utils/parse_cp_features.py
def parse_cp_features(
    feature: str, channels: list = ["DNA", "RNA", "AGP", "Mito", "ER", "mito_tubeness"]
):
    """Parses a CellProfiler feature string into its semantic components.

    This function will take a feature string and return a dictionary containing its semantic components,
    specifically: the compartment, feature group, feature type, and channel.
    If the feature string is not in a recognized format, the function will assign 'Unknown' to the non-comprehensible components.
    Channel information will be returned as 'None' where it's not applicable.

    Parameters
    ----------
    feature : str
        The CellProfiler feature string to parse.

    channels : list, optional
        A list of channel names to use when parsing the feature string. The default is ['DNA', 'RNA', 'AGP', 'Mito', 'ER', "mito_tubeness"].

    Returns
    -------
    dict
        A dictionary with the following keys: 'feature', 'compartment', 'feature_group', 'feature_type', 'channel'.
        Each key maps to the respective component of the feature string.

    Raises
    ------
    ValueError
        Raised if the input is not a string.
    """

    if not isinstance(feature, str):
        raise ValueError(f"Expected a string, got {type(feature).__name__}")

    if not isinstance(channels, list):
        raise ValueError(f"Expected a list, got {type(channels).__name__}")

    def channel_standardizer(channel):
        channel = channel.replace("Orig", "")
        return channel

    unique_token = "XUNIQUEX"
    tokenized_feature = feature
    for channel in channels:
        tokenized_channel = channel.replace("_", unique_token)
        tokenized_feature = tokenized_feature.replace(channel, tokenized_channel)

    parts = tokenized_feature.split("_")

    feature_group = parts[1]
    if parts[0] not in ["Cells", "Cytoplasm", "Nuclei", "Image"]:
        compartment = "XUNKNOWN"
        feature_group = "XUNKNOWN"
        feature_type = "XUNKNOWN"
        channel = "XUNKNOWN"
    else:
        compartment = parts[0]
        feature_group = parts[1]
        feature_type = "XNONE"  # default value
        channel = "XNONE"  # default value

        if feature_group in [
            "AreaShape",
            "Neighbors",
            "Children",
            "Parent",
            "Number",
            "Threshold",
            "ObjectSkeleton",
        ]:
            # Examples:
            # Cells,AreaShape,Zernike_2_0
            # Cells,AreaShape,BoundingBoxArea
            # Cells,Neighbors,AngleBetweenNeighbors_Adjacent
            # Nuclei,Children,Cytoplasm_Count
            # Nuclei,Parent,NucleiIncludingEdges
            # Nuclei,Number,ObjectNumber
            # Image,Threshold,SumOfEntropies_NucleiIncludingEdges
            # Nuclei,ObjectSkeleton,NumberTrunks_mito_skel

            feature_type = parts[2]

        elif feature_group == "Location":
            # Examples:
            # Cells,Location_CenterMassIntensity_X_DNA
            # Cells,Location_Center_X

            feature_type = parts[2]
            if feature_type != "Center":
                channel = parts[4]

        elif feature_group == "Count":
            # Examples:
            # Cells,Count,Cells
            pass

        elif feature_group == "Granularity":
            # Examples:
            # Cells,Granularity,15_ER
            channel = parts[3]

        elif feature_group in ["Intensity", "ImageQuality"]:
            # Examples:
            # Cells,Intensity,MeanIntensity_DNA
            # Image,ImageQuality,MaxIntensity_OrigAGP
            feature_type = parts[2]
            channel = parts[3]

        elif feature_group == "Correlation":
            # Examples:
            # Cells,Correlation,Correlation_DNA_ER
            feature_type = parts[2]
            channel = [parts[3], parts[4]]
            channel.sort()
            channel = "_".join(channel)

        elif feature_group in ["Texture", "RadialDistribution"]:
            # Examples:
            # Cells,Texture,SumEntropy_ER_3_01_256
            # Cells,RadialDistribution,FracAtD_mito_tubeness_2of16
            feature_type = parts[2]
            channel = parts[3]

        else:
            feature_group = "XUNKNOWN"
            feature_type = "XUNKNOWN"
            channel = "XUNKNOWN"

    channel = "_".join(list(map(channel_standardizer, channel.split("_"))))

    channel = channel.replace(unique_token, "_")

    return {
        "feature": feature,
        "compartment": compartment,
        "feature_group": feature_group,
        "feature_type": feature_type,
        "channel": channel,
    }
