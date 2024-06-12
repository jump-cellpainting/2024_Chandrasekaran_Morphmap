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
from typing import List, Tuple
from copairs import compute
from copairs.matching import MatcherMultilabel, UnpairedException


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
        .dropna(subset=["Metadata_broad_sample"])
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
    df = df.query("Metadata_pert_type!='empty'").reset_index(drop=True)
    df = df.dropna(subset=["Metadata_broad_sample"]).reset_index(drop=True)
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
    pandas.DataFrame
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


def compute_pairwise_cosine_similarities(
    meta,
    feats,
    pos_sameby,
    pos_diffby,
    neg_sameby,
    neg_diffby,
    multilabel_col,
    batch_size=20000,
):
    similarity_df = pd.DataFrame()

    columns = flatten_str_list(pos_sameby, pos_diffby, neg_sameby, neg_diffby)

    meta = meta.reset_index(drop=True).copy()

    matcher = MatcherMultilabel(
        *evaluate_and_filter(meta, columns), multilabel_col=multilabel_col, seed=0
    )

    pos_pairs = matcher.get_all_pairs(sameby=pos_sameby, diffby=pos_diffby)
    pos_keys = pos_pairs.keys()
    pos_counts = np.fromiter(map(len, pos_pairs.values()), dtype=np.int32)
    pos_total = sum(pos_counts)
    if pos_total == 0:
        raise UnpairedException("Unable to find positive pairs.")
    pos_pairs = np.fromiter(
        itertools.chain.from_iterable(pos_pairs.values()),
        dtype=np.dtype((np.int32, 2)),
        count=pos_total,
    )

    neg_pairs = matcher.get_all_pairs(sameby=neg_sameby, diffby=neg_diffby)
    neg_total = sum(len(p) for p in neg_pairs.values())
    if neg_total == 0:
        raise UnpairedException("Unable to find any negative pairs.")
    neg_pairs = np.fromiter(
        itertools.chain.from_iterable(neg_pairs.values()),
        dtype=np.dtype((np.int32, 2)),
        count=neg_total,
    )

    neg_pairs = np.unique(neg_pairs, axis=0)

    pos_sims = compute.pairwise_cosine(feats, pos_pairs, batch_size)

    neg_sims = compute.pairwise_cosine(feats, neg_pairs, batch_size)

    similarity_df = pd.concat(
        [
            pd.DataFrame(
                {
                    "Cosine_similarity": pos_sims,
                    "Pairs": ["Matches"] * len(pos_sims),
                }
            ),
            pd.DataFrame(
                {
                    "Cosine_similarity": neg_sims,
                    "Pairs": ["Non-matches"] * len(neg_sims),
                }
            ),
        ]
    )

    return similarity_df
