#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"


// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"

process extract_context {
    conda params.conda
    scratch true

    input:
        path(variants)

    script:
    name = "variants_context.bed"
    """
    cat ${variants} | awk -v OFS='\t' '{ print \$1,\$2,\$3 }' > variants.bed 
    bedtools getfasta -fi ${params.genome_fasta_file} -bed variants.bed -tab > ${name}
    """
}

workflow {
    out =  Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
        | collect(sort: true)
        | filterUniqVariants
        | extract_context
}