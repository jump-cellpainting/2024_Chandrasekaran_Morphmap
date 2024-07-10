import pandas as pd
import numpy as np
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_array, as_float_array
from sklearn.base import TransformerMixin, BaseEstimator
import scipy
from sklearn.metrics import average_precision_score
from sklearn.metrics.pairwise import cosine_similarity
from tqdm import tqdm
import skimage.util
import numpy as np
import matplotlib.pyplot as plt
import os
import boto3
from botocore import UNSIGNED
from botocore.client import Config


row_name_map = {
    "row": {
        "A": "01",
        "B": "02",
        "C": "03",
        "D": "04",
        "E": "05",
        "F": "06",
        "G": "07",
        "H": "08",
        "I": "09",
        "J": "10",
        "K": "11",
        "L": "12",
        "M": "13",
        "N": "14",
        "O": "15",
        "P": "16",
    }
}

orf_channel_map = {
    "DNA": "5",
    "AGP": "2",
    "Mito": "1",
    "ER": "4",
    "RNA": "3",
}

crispr_channel_map = {
    "DNA": "1",
    "AGP": "5",
    "Mito": "2",
    "ER": "4",
    "RNA": "3",
}


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
    print(f'Removed nan features: {features_to_remove}')
    df.drop(features_to_remove, axis=1, inplace=True)
    return df

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

def find_s3_subfolders(bucket_name, prefix, subfolder_name):
    # Create an S3 client without authentication for public buckets
    s3 = boto3.client("s3", config=Config(signature_version=UNSIGNED))

    # List objects in the bucket with the specified prefix
    paginator = s3.get_paginator("list_objects_v2")
    pages = paginator.paginate(Bucket=bucket_name, Prefix=prefix, Delimiter="/")

    subfolders = pd.Series()

    for page in pages:
        if "CommonPrefixes" in page:
            for common_prefix in page["CommonPrefixes"]:
                # Extract the subfolder name
                subfolder = common_prefix["Prefix"][len(prefix) :].rstrip("/")
                if subfolder:  # Ignore empty strings
                    subfolders = pd.concat([subfolders, pd.Series([subfolder])])

    subfolder_correct_name = subfolders[subfolders.str.contains(subfolder_name)][0]
    return subfolder_correct_name


def download_images(metadata):
    # Example
    # aws s3 ls --no-sign-request s3://cellpainting-gallery/cpg0016-jump/source_13/images/20221109_Run5/images/CP-CC9-R5-22/CP-CC9-R5-22_C08_T0001F008L01A02Z01C03.tif
    # aws s3 ls --no-sign-request s3://cellpainting-gallery/cpg0016-jump/source_4/images/2021_05_31_Batch2/images/BR00123530__2021-05-11T17_11_46-Measurement1/Images/r01c14f08p01-ch6sk1fk1fl1.tiff

    aws_s3_prefix = f"aws s3 sync"
    sign_request = "--no-sign-request"
    bucket = "cellpainting-gallery"
    accession_number = "cpg0016-jump"

    metadata = (
        metadata.assign(row=lambda x: x.apply(lambda y: y.Metadata_Well[0:1], axis=1))
        .assign(col=lambda x: x.apply(lambda y: y.Metadata_Well[1:], axis=1))
        .replace(row_name_map)
    )

    plane = f"p01"

    for row in tqdm(metadata.itertuples(), total=metadata.shape[0]):
        batch = row.Metadata_Batch
        perturbation = row.Metadata_JCP2022
        modality = row.modality
        source = row.Metadata_Source

        if modality == "ORF":
            row_name = f"r{row.row}"
            col_name = f"c{row.col}"
            field = f"f05"
            images_folder = f"{accession_number}/{source}/images/{batch}/images/"
            plate = find_s3_subfolders(bucket, images_folder, row.Metadata_Plate)
            channel = f"ch{orf_channel_map[row.channel]}"
            download_cmd = f'{aws_s3_prefix} {sign_request} s3://{bucket}/{images_folder}{plate}/Images/ images/{perturbation}/ --exclude "*" --include "{row_name}{col_name}{field}{plane}-{channel}*.tiff"'
            os.system(download_cmd)
            mv_cmd = f"mv images/{perturbation}/{row_name}{col_name}{field}{plane}-{channel}*.tiff images/{perturbation}/{perturbation}_{row.channel}.tiff"
            os.system(mv_cmd)
        else:
            well = f"{row.Metadata_Well}"
            field = f"F005"
            channel = f"C0{crispr_channel_map[row.channel]}"
            plate = row.Metadata_Plate
            images_folder = f"{accession_number}/{source}/images/{batch}/images/"
            download_cmd = f'{aws_s3_prefix} {sign_request} s3://{bucket}/{images_folder}{plate}/ images/{perturbation}/ --exclude "*" --include "{plate}_{well}_T0001{field}*{channel}.tif"'
            os.system(download_cmd)
            mv_cmd = f"mv images/{perturbation}/{plate}_{well}_T0001{field}*{channel}.tif images/{perturbation}/{perturbation}_{row.channel}.tiff"
            os.system(mv_cmd)

def brighten_contrast_stretch(image, low_percentile=2, high_percentile=98):
    p_low, p_high = np.percentile(image, (low_percentile, high_percentile))
    return skimage.exposure.rescale_intensity(image, in_range=(p_low, p_high))


def standardize_image(img, target_size=(1080, 1080)):
    # Resize image to fit within target size while maintaining aspect ratio
    h, w = img.shape[:2]
    scale = min(target_size[0] / h, target_size[1] / w)
    new_h, new_w = int(h * scale), int(w * scale)
    resized = skimage.transform.resize(img, (new_h, new_w), anti_aliasing=True)

    # Pad to reach target size
    pad_h = (target_size[0] - new_h) // 2
    pad_w = (target_size[1] - new_w) // 2
    pad_h_top, pad_h_bottom = pad_h, target_size[0] - new_h - pad_h
    pad_w_left, pad_w_right = pad_w, target_size[1] - new_w - pad_w

    # Use numpy.pad instead of skimage.util.pad
    padding = ((pad_h_top, pad_h_bottom), (pad_w_left, pad_w_right))
    if resized.ndim == 3:
        padding += ((0, 0),)  # Add padding for color channel if present

    padded = np.pad(resized, padding, mode="constant")

    return padded


def create_facet_grid_montage(
    images, row_labels, col_labels, grid_shape, image_labels=None, label_font_size=14
):
    # Create the montage
    m = skimage.util.montage(images, grid_shape=grid_shape)
    m_rescaled = skimage.exposure.rescale_intensity(m, out_range=(0, 1))

    # Calculate sizes
    n_rows, n_cols = grid_shape
    height, width = m_rescaled.shape[:2]
    img_height, img_width = height // n_rows, width // n_cols

    # Create a new figure
    fig, ax = plt.subplots(figsize=(12, 10))  # Adjust size as needed

    # Display the montage
    ax.imshow(m_rescaled, cmap="gray")

    # Add row labels
    for i, label in enumerate(row_labels):
        ax.text(
            -0.1,
            (i + 0.5) * img_height,
            label,
            ha="right",
            va="center",
            transform=ax.transData,
            fontsize=label_font_size,
        )

    # Add column labels
    for i, label in enumerate(col_labels):
        ax.text(
            (i + 0.5) * img_width,
            -0.1 * img_height,
            label,
            ha="center",
            va="bottom",
            transform=ax.transData,
            fontsize=label_font_size,
        )

    # Add image labels on top
    if image_labels is not None:
        for i in range(n_rows):
            for j in range(n_cols):
                ax.text(
                    (j + 0.5) * img_width,
                    i * img_height + 0.15 * img_height,
                    image_labels[i * n_cols + j],
                    ha="center",
                    va="bottom",
                    transform=ax.transData,
                    fontsize=10,
                    fontdict=dict(color="red"),
                    fontweight="bold",
                    bbox=dict(facecolor="black", edgecolor="none", pad=2),
                )

    # Add grid lines
    for i in range(1, n_rows):
        ax.axhline(y=i * img_height, color="white", linewidth=2)
    for i in range(1, n_cols):
        ax.axvline(x=i * img_width, color="white", linewidth=2)

    # Remove axes ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # Add a border around the entire plot
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("white")
        spine.set_linewidth(2)

    # Adjust layout
    plt.tight_layout()

    return fig