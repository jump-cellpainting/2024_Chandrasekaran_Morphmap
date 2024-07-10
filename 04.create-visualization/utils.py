import pandas as pd
import numpy as np

def infer_platemap(df):
    plates = df.Plate.unique().tolist()
    platemap_df = pd.DataFrame({"Plate": plates, "Platemap": np.nan})
    id = 1
    for i in range(len(plates) - 1):
        plate_1_name = plates[i]
        if (
            platemap_df.query("Plate == @plate_1_name")["Platemap"]
            .isnull()
            .values.any()
        ):
            platemap_name = f"platemap_{id}"
            id += 1
            platemap_df.loc[platemap_df.Plate == plate_1_name, "Platemap"] = (
                platemap_name
            )
            plate_1 = df.query("Plate == @plate_1_name")[["Metadata_Well", "Symbol"]]
            for j in range(i + 1, len(plates)):
                plate_2_name = plates[j]
                if (
                    platemap_df.query("Plate == @plate_2_name")["Platemap"]
                    .isnull()
                    .values.any()
                ):
                    plate_2 = df.query("Plate == @plate_2_name")[
                        ["Metadata_Well", "Symbol"]
                    ]
                    merged_plate = pd.merge(
                        plate_1, plate_2, on=["Metadata_Well", "Symbol"], how="inner"
                    )
                    if len(merged_plate) > 0.95 * len(plate_2):
                        platemap_df.loc[
                            platemap_df.Plate == plate_2_name, "Platemap"
                        ] = platemap_name

    return platemap_df

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