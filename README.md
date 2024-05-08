# nf-annotations
Repository with nextflow scripts for motif enrichment calculation and phenotypes annotation<br><br>
- Global cluster and executor parameters are specified in ```nextflow.config``` file.
- All external files and parameters are specified in ```params.config```
## Args in params.config
### Common params
- `outdir` - directory to save results into
- `conda` - path to existing conda environment


### LDSC params
- `base_ann_path`, `frqfiles`, `gtfiles`, `weights` - input files for LDSC regresion, see [here](https://github.com/bulik/ldsc)
- `ldsc_scripts_path` - path to cloned LDSC repo (won't be required once dockerized)
- `annotations_dir` - directory with BED annotations files. Expected to have `#chr`, `start`, `end` fields. Should contain only SNPs to include in LDSC.
- `phenotypes_meta` - tsv file listing all phenotypes to test. Should have `phen_id` and `sumstats_file` (path to sumstats file) fields, other fields are ignored.
- `ldsc_conda` - conda with installed LDSC (won't be required once dockerized)
- `tested_snps` - Intersection of SNPs from HM3 and data on which sumstats are calculated

### CAVs motif enrichment
1) Fill in parameters in ```params.config``` file, section ```Motif enrichment for index```
2) Run motif enrichment on CAV p-value files (see [here](https://github.com/wishabc/nf-babachi)) with 
```nextflow run <src-dir>/motif_enrichment.nf -profile Altius -resume -entry scanWithMoods``` (just once for motif collection)
```nextflow run <src-dir>/motif_enrichment.nf -profile Altius -resume``` (just once for set of unique variants)
```nextflow run <src-dir>/motif_enrichment.nf -profile Altius -resume -entry cavsEnrichment```

## Phenotypes associations analysis
### LDSC
1) Fill in parameters in ```params.config``` file, section ```LDSC data```
2) (optional) If baseline LD scores are not calculated, run ```nextflow run <src-dir>/phenotype_annotation.nf -profile Altius -resume -entry calcBaseline```. They will be saved to the output directory.
2) Run LDSC on annotation files with
```nextflow run <src-dir>/phenotype_annotation.nf -profile Altius -resume```
