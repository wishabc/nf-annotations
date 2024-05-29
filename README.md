# nf-annotations
Repository with nextflow scripts for motif and GWAS analysis<br><br>
- Global cluster and executor parameters are specified in ```nextflow.config``` file.
- All external files and parameters are specified in ```params.config```
## Args in params.config
### Common params
- `outdir` - directory to save results into
- `conda` - path to existing conda environment

# Quick start
- preprocess_ebi_sumstats.nf - download, reformat and run munge_sumstats on GWAS data downloaded from ebi portal. **You'd need to annotate metadata with `n_samples` field after downloading**
- preprocess_ukbb.nf - download, reformat and run munge_sumstats on GWAS data downloaded from [https://www.nealelab.is/uk-biobank](https://www.nealelab.is/uk-biobank). **You'd need to annotate metadata with `n_samples` field after downloading**
- Create a metadata file for all phenotypes of interest, add paths to files reformatted with munge_sumstats
- Run `ldsc.nf -entry calcBaseline` to get h^2 estimates for phenotypes. The pipeline will create `output/ldsc_enrichments_results.tsv` file with h^2 estimates.
- (optional) Filter your metadata by h^2 estimates
- Fill in `baseline_ld` option in params.config with path to `output/ldsc/l2/baseline.`
- Run LDSC (`ldsc.nf` script) using custom annotations or matrix (use respective entry).

### LDSC params
- `base_ann_path`, `frqfiles`, `gtfiles`, `weights` - input files for LDSC regresion, see [here](https://github.com/bulik/ldsc)
- `ldsc_scripts_path` - path to cloned LDSC repo (won't be required once dockerized)
- `phenotypes_meta` - tsv file listing all phenotypes to test. Should have `phen_id` and `sumstats_file` (path to sumstats file) fields, other fields are ignored.
- `ldsc_conda` - conda with installed LDSC (won't be required once dockerized)
- `tested_snps` - Intersection of SNPs from HM3 and data on which sumstats are calculated

- `annotations_dir` - directory with BED annotations files. Expected to have `#chr`, `start`, `end` fields. Should contain only SNPs to include in LDSC.

### Motif enrichment
1) Fill in parameters in ```params.config``` file, section ```Motif enrichment for index```
2) Run motif enrichment on CAV p-value files (see [here](https://github.com/wishabc/nf-babachi)) with 
```nextflow run <src-dir>/motif_enrichment.nf -profile Altius -resume -entry scanWithMoods``` (just once for motif collection)
```nextflow run <src-dir>/motif_enrichment.nf -profile Altius -resume``` (just once for set of unique variants)
```nextflow run <src-dir>/motif_enrichment.nf -profile Altius -resume -entry cavsEnrichment```

# Workflows:
## Phenotypes associations
### LDSC
1) Fill in parameters in ```params.config``` file, section ```LDSC data```
2) (optional) If baseline LD scores are not calculated, run ```nextflow run <src-dir>/phenotype_annotation.nf -profile Altius -resume -entry calcBaseline```. They will be saved to the output directory.
2) Run LDSC on annotation files with
```nextflow run <src-dir>/phenotype_annotation.nf -profile Altius -resume```
