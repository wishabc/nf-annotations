#!/usr/bin/env nextflow

process process_mutation_rates {
    tag "${vcf.simpleName}"
    scratch true
    conda params.conda

    input:
        path vcf

    output:
        path name

    script:
    name = "${vcf.simpleName}.bed"
    """
    echo -e "`head -1 ${params.variants}`\tchr\tstart_mr\tend_mr\tref_mr\talt_mr\tmut_rates_roulette\tmut_rates_gnomad" > tmp.bed
    
    bcftools query -f"%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%INFO/MR\t%INFO/MG\n" \
        ${vcf} | awk '{print "chr"\$0}' | bedtools intersect \
        -a ${params.variants} -b stdin -wa -wb >> tmp.bed
    
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
        awk 'NR>1' \$file > tmp.bed
    done
    head -1 ${bed_files[0]} > ${name}
    sort-bed tmp.bed >> ${name}
    """
}

workflow annotateMutationRates {
    take:
        vcfs
    main:
        out = process_mutation_rates(vcfs) | collect(sort: true) | merge_and_sort
    emit:
        out
}

workflow {
    Channel.fromPath("${params.vcfs_dir}*.vcf.gz")
        | annotateMutationRates
}