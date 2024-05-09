#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import utils
from copairs.map import mean_average_precision, multilabel
import numpy as np


# In this notebook, I want to compare six different scenarios for either computing the consensus ORF profile or how to handle duplicate profiles for downstream analysis in copairs and identify the best scenario or the one that makes the most sense. For doing this, I will be using the `ORF` profiles and two sets of annotations, CORUM complex ID and HGNC gene group ID. Here are the six scenarios:
# - s1: group profiles by JCPID
# - s2: group profiles by GeneID
# - s3: group profiles by JCPID but use gene ID as pos_diffby in copairs
# - s4: Keep the most replicable reagent when there are many targeting the same gene
# - s5: Keep a random ORF reagent
# - s6: Keep the reagent with the longest insert and highest number of proteins matching and then randomly choose one reagent if there are multiple ones with the same number of protein matching and the same insert length

# In[ ]:


gene_ID_col = "Metadata_NCBI_Gene_ID"
reagent_ID_col = "Metadata_JCP2022"

scenarios = {
    "s1": f"{reagent_ID_col}",
    "s2": f"{gene_ID_col}",
    "s3": f"{reagent_ID_col}",
    "s4": f"{reagent_ID_col}",
    "s5": f"{reagent_ID_col}",
    "s6": f"{reagent_ID_col}",
}

annotations = {
    "corum": {
        "annotation_col": "Metadata_complexname",
        "multi_label_col": "Metadata_corum_complex_list",
    },
    "gene-group": {
        "annotation_col": "Metadata_gene_group_id",
        "multi_label_col": "Metadata_gene_group_list",
    },
}

profiles = {
    "ORF": "wellpos_cc_var_mad_outlier_featselect_sphering_harmony",
}

batch_size = 20000
null_size = 20000
fdr = 0.05


# In[ ]:


orf_metadata_df = pd.read_csv(
    "../00.download-and-process-annotations/output/orf_metadata.tsv.gz", sep="\t"
)[[f"{reagent_ID_col}", f"{gene_ID_col}", "Metadata_Insert_Length", "Metadata_Prot_Match"]]


# In[ ]:


ap_df = pd.DataFrame()
map_df = pd.DataFrame()

for scenario in scenarios:
    grouping_col = scenarios[scenario]
    for profile in profiles:
        profile_filename = f"profiles_{profiles[profile]}.parquet"
        phenotypic_activity_filename = f"phenotypic-activity-{profiles[profile]}.csv.gz"
        df = pd.read_parquet(f"../profiles/{profile_filename}")
        phenotypic_activity_df = (
            pd.read_csv(
                f"output/{phenotypic_activity_filename}",
                usecols=[
                    f"{reagent_ID_col}",
                    "below_corrected_p",
                    "mean_average_precision",
                ],
            )
            .rename(
                columns={"mean_average_precision": "Metadata_mean_average_precision"}
            )
            .query("below_corrected_p==True")
            .reset_index(drop=True)
            .drop(columns="below_corrected_p")
        )
        df = df.merge(phenotypic_activity_df, on=f"{reagent_ID_col}", how="inner")
        df = df.merge(orf_metadata_df, on=f"{reagent_ID_col}", how="inner")

        for annotation in annotations:
            annotation_col = annotations[annotation]["annotation_col"]
            multi_label_col = annotations[annotation]["multi_label_col"]
            annotation_df = (
                pd.read_csv(
                    f"../00.download-and-process-annotations/output/orf_metadata.tsv.gz",
                    sep="\t",
                    usecols=[
                        f"{reagent_ID_col}",
                        f"{annotation_col}",
                    ],
                )
                .dropna()
                .assign(col=lambda x: x[f"{annotation_col}"].str.split("|"))
                .rename(columns={"col": f"{multi_label_col}"})
                .drop(columns=f"{annotation_col}")
            )

            annotated_df = df.merge(annotation_df, on=f"{reagent_ID_col}", how="inner")

            pos_sameby = [f"{multi_label_col}"]
            neg_sameby = []
            neg_diffby = [f"{multi_label_col}"]

            if not scenario == "s3":
                pos_diffby = []
            else:
                pos_diffby = [f"{gene_ID_col}"]

            consensus_df = utils.consensus(annotated_df, grouping_col)

            if scenario == "s4":
                consensus_df = consensus_df.sort_values(
                    "Metadata_mean_average_precision", ascending=False
                ).drop_duplicates(subset=[gene_ID_col])
            elif scenario == "s5":
                consensus_df = consensus_df.drop_duplicates(subset=[gene_ID_col])
            elif scenario == "s6":
                consensus_df = consensus_df.sort_values(
                    by=["Metadata_Insert_Length", "Metadata_Prot_Match"],
                    ascending=False,
                ).drop_duplicates(subset=[gene_ID_col])

            print(
                f"Scenario: {scenario} | Profile: {profile} | Annotation: {annotation} | Reagents: {consensus_df[reagent_ID_col].nunique()} | Genes: {consensus_df[gene_ID_col].nunique()}"
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
                result, pos_sameby, null_size, threshold=fdr, seed=12527
            )

            result["scenario"] = scenario
            result["annotation"] = annotation
            result["profile"] = profile

            agg_result["scenario"] = scenario
            agg_result["annotation"] = annotation
            agg_result["profile"] = profile

            ap_df = pd.concat([ap_df, result], ignore_index=True)
            map_df = pd.concat([map_df, agg_result], ignore_index=True)

ap_df.to_parquet(f"output/find-consensus-profiles-ap.parquet")
map_df.to_parquet(f"output/find-consensus-profiles-map.parquet")

