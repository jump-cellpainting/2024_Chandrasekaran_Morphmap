### Human Protein Atlas data 
A number of gene annotations were downloaded from the Human Protein Atlas webpage (https://www.proteinatlas.org/download/proteinatlas.tsv.zip) (v23.0). The output was saved in the folder `input/proteinatlas.tsv.zip`. The file was processed in the notebook `00.0.process-human-protein-atlas-annotations.ipynb`. Four output files were then written to the `output` folder.

- `human_protein_atlas_class_annotations.tsv.gz`: Protein class annotations.
- `human_protein_atlas_location_annotations.tsv.gz`: Subcellular location annotations.
- `human_protein_atlas_disease_annotations.tsv.gz`: Disease annotations.
- `human_protein_atlas_go_annotations.tsv.gz`: `Biological Process` and `Molecular Function` GO terms.

### CORUM protein complexes
CORUM protein complex annotations were downloaded from the [PERISCOPE project repo](https://github.com/broadinstitute/2022_PERISCOPE/tree/main/common_files/CORUM_humanComplexes.txt) and stored as `input/CORUM_humanComplexes.txt`. The annotations were processed in the notebook `01.process-corum-data.ipynb` and the processed file was written to `output/CORUM_proteincomplex_annotations.tsv`. The file contains the following columns

- `gene`: Gene symbol.
- `ComplexName`: Name of the protein complex.

### Infection efficiency metadata
Infection effiency data was obtained from GPP (`input/JUMP-ORF-Infection-Efficiencies.xlsx`) and processed in the notebook `02.0.process-viral-infection-efficiency-data.ipynb`. The processed file was written to the file `output/gpp-infection-efficiency.csv.gz`. The file contains the following columns

- `Batch`: Batch number.
- `Batch_name`: Batch name.
- `Plate_name`: Plate name.
- `Well`: Well name.
- `Row`: Row name.
- `Column`: Column number.
- `broad_sample`: Broad sample name.
- `Plus_Blast`: Plus blasticidin count.
- `Minus_Blast`: Minus blasticidin count.
- `I.E_percentage`: Infection efficiency percentage (`Plus_Blast`/`Minus_Blast`).
- `Minus_Avg_Background`: `I.E_percentage` - Background (value changes depending on the batch).

### HGNC Approved symbol
HGNC approved symbols were downloaded from https://genenames.org/download/custom/. The file was stored in `input/hgnc_approved_symbol.tsv`. The file was then processed in the notebook `03.process-hgnc-approved-symbol.ipynb` and then stored in `output/hgnc_approved_symbol_processed.tsv`. The file contains the following columns

-`Approved_symbol`: HGNC approved symbol
-`Previous_symbols`: `|` separatd previous symbols
-`Gene_group_ID`: Gene group ID
-`NCBI_Gene_ID`: NCBI gene ID

### NCBI dataset
NCBI dataset was downloaded from https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/9606/?gene_type=protein-coding,small-rnas,pseudogenes,non-coding. The files were stored in `input/ncbi_dataset.tsv`. The file was then processed in the notebook `04.process-ncbi-dataset.ipynb` and then stored in `output/ncbi_dataset_processed.tsv`. The file contains the following columns

- `NCBI_Gene_ID`: NCBI gene ID
- `Approved_Symbol`: HGNC approved symbol

### ORF Experiment metadata
This file was created in a private repository (https://github.com/jump-cellpainting/morphmap/blob/main/00.0.explore-data/0.create-experiment-metadata.ipynb). The file contains the following columns

- `Batch`: batch name
- `Batch_name`: Full batch name
- `Plate_Map_Name`: Platemap name
- `Assay_Plate_Barcode`: Plate name
- `Type`: Type of plate - ORF, Untreated, Target2, Target1
- `Anomaly`: Any anamalous behavior observed during data generation or analysis
- `Density`: Cell seeding density
- `Polybrene`: Indicates whether Polybrene is present or absent

### Gene ontology annotation
The annotation files were downloaded as described in `06.process-gene-ontology-annotations.ipynb` and stored in the `input/` folder. Then they were processed in the notebook `06.process-gene-ontology-annotations.ipynb` and the processed file was written to `output/gene-ontology-annotations.tsv.gz`. The file contains the following columns.

- `gene`: Gene symbol.
- `go_bp_list`: A list of `|` separated Biological Process GO terms.
- `go_cc_list`: A list of `|` separated Cellular Component GO terms.
- `go_mf_list`: A list of `|` separated Molecular Function GO terms.

### Homologue data
Homologue data was downloaded from https://ftp.ncbi.nih.gov/pub/HomoloGene/build68/homologene.data and stored in the `input/` folder. The was processed in the notebook `07.process-homologue-data.ipynb` and the processed file was written to `output/homologue_data.tsv.gz`. The file contains the following columns (the columns names are not in the file).

- HID (HomoloGene group id)
- Taxonomy ID
- Gene ID
- Gene Symbol
- Protein gi
- Protein accession

### Wikipathway annotations
Wikipathway annotations were downloaded as described in `08.process-wikipathway-annotations.ipynb` and stored in the `input/` folder. Then they were processed in the notebook `08.process-wikipathway-annotations.ipynb` and the processed file was written to `output/wikipathway-annotations.tsv.gz`. The file contains the following columns.

- `gene`: HGNC gene symbol
- `wikipathway`: Wikipathway ID

### U2OS gene essentiality data form DepMap
Essentiality data is downloaded and stored outside the repository, as described in `09.download-and-process-depmap-data.ipynb`. Then the data was processed using the same notebook to filter out all cell lines that are not U2OS. The processed file was written to `output/CRISPRGeneDependency_U2OS.csv.gz`. The file contains the following columns.

- `geneID`: NCBI gene ID.
- `Symbol`: Gene symbol.
- `probability`: Probability to determine essentiality from DepMap.

### U2OS expression data from DepMap
Expression data is downloaded and stored outside the repository, as described in `09.download-and-process-depmap-data.ipynb`. The processed file was written to `output/OmicsExpressionProteinCodingGenesTPMLogp1_U2OS.csv.gz`. The file contains the following columns.

- `geneID`: NCBI gene ID.
- `Symbol`: Gene symbol.
- `TPM`: Transcripts per million (TPM) normalized expression value.

### U2OS gene effect data  from DepMap
Gene effect data is downloaded and stored outside the repository, as described in `09.download-and-process-depmap-data.ipynb`. The processed file was written to `output/CRISPRGeneEffect_U2OS.csv.gz`. The file contains the following columns.

- `geneID`: NCBI gene ID.
- `Symbol`: Gene symbol.
- `effect`: Gene effect
