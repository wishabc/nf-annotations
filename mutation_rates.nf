
process process_mutation_rates {

    tag "${vcf.simpleName}"
    scratch true

    input:
        path vcf

    output:
        path name

    script:
    name = "${vcf.simpleName}.bed"
    """
    echo -e "#chr\tstart\tend\tref\talt\tmut_rates_roulette\tmut_rates_gnomad\t`head -1 ${params.variants}`" > tmp.bed
    
    bcftools query -f"%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%INFO/MR\t%INFO/MG\n" \
        ${vcf} | awk '{print "chr"\$0}' | bedtools intersect \
        -a ${params.variants} -b stdin -wa -wb >> tmp.bed
    
    python3 $moduleDir/bin/filter_variants.py tmp.bed ${name}
    """
}