# nf-annotations
Repository with nextflow scripts for motif enrichment calculation and phenotypes annotation<br><br>
- Global cluster and executor parameters are specified in ```nextflow.config``` file.
- All external files and parameters are specified in ```params.config```
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
2) Run LDSC on index file with
```nextflow run <src-dir>/phenotype_annotation.nf -profile Altius -resume```
