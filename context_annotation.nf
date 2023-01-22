#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"


// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"
params.window = 20


process extract_context {
    conda params.conda
    scratch true
    publishDir "${params.outdir}"

    input:
        path variants
    output:
        path name

    script:
    name = "variants_context.bed"
    """
    cat ${variants} | awk -v OFS='\t' '{ print \$1,\$2-${params.window},\$3+${params.window} }' > variants.bed 
    bedtools getfasta -fi ${params.genome_fasta_file} -bed variants.bed -bedOut \
        | awk -v OFS='\t' '{ print \$1,\$2+${params.window},\$3-${params.window},\$4 }' > ${name}
    """
}

workflow {
    out =  Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
        | collect(sort: true)
        | filterUniqVariants
        | extract_context
}