#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"


process process_mutation_rates {
    tag "${vcf.simpleName}"
    scratch true
    conda params.conda
    label "highmem"

    input:
        path vcf
        path variants_file

    output:
        path name

    script:
    name = "${vcf.simpleName}.bed"
    """
    echo -e "#chr\tstart\tend\tID\tref\talt\tchr\tstart_mr\tend_mr\tref_mr\talt_mr\tmut_rates_roulette\tmut_rates_gnomad" > tmp.bed
    
    bcftools query -f"%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%INFO/MR\t%INFO/MG\n" \
        ${vcf} | awk '{print "chr"\$0}' | bedtools intersect \
        -a ${variants_file} -b stdin -sorted -wa -wb >> tmp.bed
    
    python3 $moduleDir/bin/filter_variants.py tmp.bed ${name}
    """
}

process merge_and_sort {
    publishDir "${params.outdir}/mutation_rates"
    conda params.conda
    scratch true

    input:
        path bed_files
    
    output:
        path name

    script:
    name = "mut_rates.annotation.bed"
    """
    for file in ${bed_files}; do
        awk 'NR>1' \$file >> tmp.bed
    done
    head -1 ${bed_files[0]} > ${name}
    sort-bed tmp.bed >> ${name}
    """
}

workflow annotateMutationRates {
    take:
        vcfs
        pval_file
    main:
        out = process_mutation_rates(vcfs, pval_file) 
            | collect(sort: true)
            | merge_and_sort
    emit:
        out
}

workflow {
    pval_file = Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
        | collect(sort: true)
        | filterUniqVariants
    mut_rate_vcfs = Channel.fromPath("${params.vcfs_dir}/*.vcf.gz")
    annotateMutationRates(mut_rate_vcfs, pval_file)
}