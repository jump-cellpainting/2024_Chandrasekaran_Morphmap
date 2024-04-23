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
    "ORF": {
        "ORF": "wellpos_cc_var_mad_outlier_featselect_sphering_harmony",
        "ORF-CRISPR-pipeline": "wellpos_var_mad_int_featselect_harmony_PCA",
    },
    "CRISPR": {
        "CRISPR": "wellpos_var_mad_int_featselect_harmony_PCA_corrected",
        "CRISPR-ORF-pipeline": "wellpos_cc_var_mad_outlier_featselect_sphering_harmony_PCA_corrected",
    },
}

retrieval_label = "disease_association"
multi_label_col = "Metadata_disease_list"
annotation_col = "Metadata_disease_involvement"


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
    )[["Metadata_NCBI_Gene_ID", f"{annotation_col}"]]
    .dropna()
    .assign(
        col=lambda x: list(x[f"{annotation_col}"].str.split("|"))
    ).rename(columns={"col": f"{multi_label_col}"})
)

orf_annotation_df.Metadata_NCBI_Gene_ID = (
    orf_annotation_df.Metadata_NCBI_Gene_ID.astype(int).astype(str)
)


# In[ ]:


crispr_annotation_df = (
    pd.read_csv(
        "../00.download-and-process-annotations/output/crispr_metadata.tsv.gz", sep="\t"
    )[["Metadata_NCBI_Gene_ID", f"{annotation_col}"]]
    .dropna()
    .assign(col=lambda x: list(x[f"{annotation_col}"].str.split("|")))
    .rename(columns={"col": f"{multi_label_col}"})
)

crispr_annotation_df.Metadata_NCBI_Gene_ID = (
    crispr_annotation_df.Metadata_NCBI_Gene_ID.astype(int).astype(str)
)


# In[ ]:


map_df = pd.DataFrame()

for modality in all_profiles:
    for name, pipeline in all_profiles[modality].items():
        print(f"Profile type: {name}")
        parquet_file_name = f"profiles_{pipeline}.parquet"
        phenotypic_activity_file_name = f"phenotypic-activity-{pipeline}.csv.gz"

        df = pd.read_parquet(f"../profiles/{parquet_file_name}")
        phenotypic_activity_df = pd.read_csv(
            f"output/{phenotypic_activity_file_name}"
        ).query("below_corrected_p==True")
        df = df.merge(phenotypic_activity_df, on="Metadata_JCP2022", how="inner").query(
            "Metadata_pert_type!='control'"
        )
        consensus_df = utils.consensus(df, "Metadata_JCP2022")

        if modality == "ORF":
            consensus_df["Metadata_NCBI_Gene_ID"] = (
                consensus_df["Metadata_NCBI_Gene_ID"].astype(str)
            )
            consensus_df = consensus_df.merge(
                orf_annotation_df,
                on="Metadata_NCBI_Gene_ID",
                how="inner",
            )
        elif modality == "CRISPR":
            consensus_df["Metadata_NCBI_Gene_ID"] = (
                consensus_df["Metadata_NCBI_Gene_ID"].astype(int).astype(str)
            )
            consensus_df = consensus_df.merge(
                crispr_annotation_df,
                on="Metadata_NCBI_Gene_ID",
                how="inner",
            )

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

        agg_result["Profile_type"] = name

        map_df = pd.concat([map_df, agg_result], ignore_index=True)

        map_df.to_parquet(f"output/{retrieval_label}_retrieval.parquet")

