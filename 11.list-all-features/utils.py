import pandas as pd

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
