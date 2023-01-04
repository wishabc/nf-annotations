# nf-annotations
Repository with nextflow scripts for motif enrichment calculation and phenotypes annotation<br><br>
- Global cluster and executor parameters are specified in ```nextflow.config``` file.
- All external files and parameters are specified in ```params.config```
## Args in params.config
### Common params
- `outdir` - directory to save results into
- `conda` - path to existing conda environment
- `pval_file_dir` - directory with BED pvalue files files (output of [nf-babachi pipeline](https://github.com/wishabc/nf-babachi)) 
- `genome_fasta_file` - reference genome fasta
- `moods_scans_dir` - directory to save moods motif hits to. If exists, the pipeline will look for ```<motif_id>.moods.log.bed.gz``` files in this directory
### Moods params
 - `alt_fasta_file` - Fasta file with all SNPs encoded as IUPAC symbols (see [MOODS README](https://github.com/jhkorhonen/MOODS/wiki/Getting-Started))
- `bg_file` - file with background frequencies of nucleotides (pA pC pG pT), new line separated. Ignored if file does not exist
- `motifs_list` -  tsv file listing all the motifs, expected to have `motif` and `motif_file` fields, other fields are ignored.
- `motif_pval_tr` - threshold to scan motifs at

### Index motif enrichment params
- `index_file` - path to index file of all the DHSs
- `sample_names` - sample names, new line separated
- `step` - # of samples in a chunk, should be more than 1
- `binary_matrix` - binary matrix of size DHS_count x #_of_samples

### LDSC params
- `base_ann_path`, `frqfiles`, `gtfiles`, `weights` - input files for LDSC regresion, see [here](https://github.com/bulik/ldsc)
- `ldsc_scripts_path` - path to cloned LDSC repo (won't be required once dockerized)
- `annotations_dir` - directory with BED annotations files. Expected to have `#chr`, `start`, `end` fields. Should contain only SNPs to include in LDSC.
- `phenotypes_meta` - tsv file listing all phenotypes to test. Should have `phen_id` and `sumstats_file` (path to sumstats file) fields, other fields are ignored.
- `ldsc_conda` - conda with installed LDSC (won't be required once dockerized)
- `tested_snps` - Intersection of SNPs from HM3 and data on which sumstats are calculated

## Motif enrichment
### Moods scans
Moods scans are automatically created if you provide non-existent directory in ```moods_scans_dir``` parameter. The resulting moods files are saved into the ```moods_scans_dir```.<br> If directory exists, the pipeline will look for ```<motif_id>.moods.log.bed.gz``` files in this directory.<br><br> Before running Moods, please fill in parameters in ```params.config``` file, section ```Moods params```.

### CAVs motif enrichment
1) Fill in parameters in ```params.config``` file, section ```Motif enrichment for index```
2) Run motif enrichment on CAV p-value files (see [here](https://github.com/wishabc/nf-babachi)) with 
```nextflow run <src-dir>/motif_enrichment.nf -profile Altius -resume```

### Index motif enrichment
1) Fill in parameters in ```params.config``` file, section ```Motif enrichment for index```
2) Run motif enrichment on index file with
```nextflow run <src-dir>/motif_enrichment.nf -profile Altius -entry indexEnrichment -resume```

## Phenotypes associations analysis
### UKBB associations
1) Fill in parameters in ```params.config``` file, section ```LDSC data```
2) Run LDSC on annotation files with
```nextflow run <src-dir>/phenotype_annotation.nf -profile Altius -resume```


### ClinVAR, Finemapping, GRASP, GTEx, EBI-GWAS, pheWAS associations annotation
❗Input files are set up only at Altius cluster.<br><br>
Run phenotype annotation on the BED files in `pval_file_dir` with
```nextflow run <src-dir>/phenotype_annotation.nf -profile Altius -entry annotateWithPheno -resume```
<br>The files for this step should have `#chr`, `start`, `end`, `ID` columns
