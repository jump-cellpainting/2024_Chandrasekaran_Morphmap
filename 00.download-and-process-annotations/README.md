### Human Protein Atlas data 
A number of gene annotations were downloaded from the Human Protein Atlas webpage (https://www.proteinatlas.org/download/proteinatlas.tsv.zip) (v23.0). The output was saved in the folder `input/proteinatlas.tsv.zip`. The file was processed in the notebook `00.process-human-protein-atlas-annotations.ipynb`. Four output files were then written to the `output` folder.

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