We use the same pipeline to process both CRISPR and ORF profiles. The following are the steps performed

- Mean aggregate single cell profiles to well-level profiles
- wellpos: Well position correction
- cc: Cell count regression
- var: variance thresholding
- mad: Robust MAD normalization
- outlier: Outlier removal
- featselect: Feature selection
- sphering: Sphering
- harmony: Harmony correction

In addition to these steps, the following steps are performed on the CRISPR profiles
- PCA: Principal component analysis
- corrected: Chromosome arm correction

TODO describe each step in more detail.

The processed profiles were downloaded using the following commands

```bash
aws s3 cp s3://cellpainting-gallery/cpg0016-jump-assembled/source_all/workspace/profiles/jump-profiling-recipe_2024_a917fa7/ORF/profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony/profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony.parquet ../profiles/

aws s3 cp cpg0016-jump-assembled/source_all/workspace/profiles/jump-profiling-recipe_2024_a917fa7/CRISPR/profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony_PCA_corrected/profiles_wellpos_cc_var_mad_outlier_featselect_sphering_harmony_PCA_corrected.parquet ../profiles/
```