#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import utils
import numpy as np
import warnings
from tqdm import tqdm
from copairs.map import mean_average_precision, multilabel
import logging


# In[ ]:


all_profiles = {
    "ORF": "wellpos_cc_var_mad_outlier_featselect_sphering_harmony",
    "CRISPR": "wellpos_cc_var_mad_outlier_featselect_sphering_harmony_PCA_corrected",
}

retrieval_label = "wikipathway"
multi_label_col = "Metadata_wikipathway_list"
annotation_col = "Metadata_wikipathway"


batch_size = 20000
null_size = 20000
fdr = 0.05

pos_sameby = [f"{multi_label_col}"]
pos_diffby = []
neg_sameby = []
neg_diffby = [f"{multi_label_col}"]


# Read protein class annotations

# In[ ]:


orf_annotation_df = (
    pd.read_csv(
        "../00.download-and-process-annotations/output/orf_metadata.tsv.gz", sep="\t"
    )[
        [
            "Metadata_JCP2022",
            f"{annotation_col}",
            "Metadata_pert_type",
            "Metadata_NCBI_Gene_ID",
        ]
    ]
    .dropna()
    .assign(col=lambda x: list(x[f"{annotation_col}"].str.split("|")))
    .rename(columns={"col": f"{multi_label_col}"})
)


# In[ ]:


crispr_controls_df = pd.DataFrame(
    {
        "Metadata_JCP2022": ["JCP2022_805264", "JCP2022_800001", "JCP2022_800002"],
        "Metadata_pert_type": ["poscon", "negcon", "negcon"],
    },
    index=[0, 1, 2],
)

crispr_annotation_df = (
    pd.read_csv(
        "../00.download-and-process-annotations/output/crispr_metadata.tsv.gz", sep="\t"
    )[["Metadata_JCP2022", f"{annotation_col}", "Metadata_NCBI_Gene_ID"]]
    .merge(crispr_controls_df, on="Metadata_JCP2022", how="left")
    .fillna(value={"Metadata_pert_type": "trt"})
    .dropna()
    .assign(col=lambda x: list(x[f"{annotation_col}"].str.split("|")))
    .rename(columns={"col": f"{multi_label_col}"})
)


# In[ ]:


map_df = pd.DataFrame()

for modality in all_profiles:
    print(f"Modality: {modality}")
    parquet_file_name = f"profiles_{all_profiles[modality]}.parquet"
    phenotypic_activity_file_name = (
        f"phenotypic-activity-{all_profiles[modality]}.csv.gz"
    )

    df = pd.read_parquet(f"../profiles/{parquet_file_name}")
    phenotypic_activity_df = pd.read_csv(
        f"output/{phenotypic_activity_file_name}"
    ).query("below_corrected_p==True")
    df = df.merge(phenotypic_activity_df, on="Metadata_JCP2022", how="inner")

    if modality == "ORF":
        df = df.merge(
            orf_annotation_df,
            on="Metadata_JCP2022",
            how="inner",
        )
    elif modality == "CRISPR":
        df = df.merge(
            crispr_annotation_df,
            on="Metadata_JCP2022",
            how="inner",
        )

    consensus_df = utils.consensus(df, "Metadata_NCBI_Gene_ID")

    metadata_df = utils.get_metadata(consensus_df)
    feature_df = utils.get_featuredata(consensus_df)
    feature_values = feature_df.values

    result = multilabel.average_precision(
        metadata_df,
        feature_values,
        pos_sameby,
        pos_diffby,
        neg_sameby,
        neg_diffby,
        batch_size=batch_size,
        multilabel_col=multi_label_col,
    )

    agg_result = mean_average_precision(
        result, pos_sameby, null_size, threshold=0.05, seed=12527
    )

    agg_result["Modality"] = modality

    map_df = pd.concat([map_df, agg_result], ignore_index=True)

    map_df.to_parquet(f"output/{retrieval_label}_retrieval.parquet")

