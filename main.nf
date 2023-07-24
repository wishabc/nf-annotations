#!/usr/bin/env nextflow
//include { annotateLD } from "./ld_scores"
// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"

params.phenotypes_data = "/home/sabramov/phenotypes_data"


process process_mutation_rates {
    tag "${vcf.simpleName}"
    scratch true
    conda params.conda
    label "highmem"

    input:
        tuple  path(vcf), path(variants_file)

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
    publishDir "${params.outdir}"
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
    cat ${variants} \
        | awk -v OFS='\t' '{ print \$1,\$2-${params.window},\$3+${params.window} }' \
        | uniq > variants.bed 
    bedtools getfasta -fi ${params.genome_fasta_file} -bed variants.bed -bedOut \
        | awk -v OFS='\t' '{ print \$1,\$2+${params.window},\$3-${params.window},\$4 }' > ${name}
    """
}

// Annotates with pheWAS, clinvar, finemapping, grasp, ebi-gwas phenotypes
process annotate_with_phenotypes {
    conda params.conda
    publishDir "${params.outdir}"

    input:
        path pval_file

    output:
        path name

    script:
    name = "phenotypes_ann.bed"
    """
    python3 $moduleDir/bin/annotate_with_phenotypes.py ${params.phenotypes_data} ${pval_file} ${name}
    """
}


process filter_tested_variants {
    conda params.conda
    scratch true

    input:
        path pval_files

    output:
        path name

    script:
    // Expected all files to be in the same format
    command = pval_files[0].extension == 'gz' ? 'zcat' : 'cat'
    name = pval_files.size() > 1 ? "unique_variants.bed" : "${pval_files[0].simpleName}.bed"
    """
    ${command} ${pval_files} \
        | awk -v OFS='\t' -v col='is_tested' \
            'NR==1 {
                for(i=1;i<=NF;i++){
                    if (\$i==col){
                        c=i;
                        break
                    }
                }
            }
            ((NR>1) && (\$c == "True")) {
                print \$1,\$2,\$3,\$4,\$5,\$6
            }' \
        | sort-bed - \
        | uniq > ${name}
    """
}

workflow filterTestedVariants {
    take:
        pval_file
    main:
        out = filter_tested_variants(pval_file)
    emit:
        out
}


workflow mutationRates {
    take:
        data
    main:
        out = Channel.fromPath("${params.vcfs_dir}/*.vcf.gz")
            | combine(data)
            | process_mutation_rates
            | collect(sort: true)
            | merge_and_sort
    emit:
        out
}

workflow {
    sample_wise_pvals = Channel.fromPath("${params.by_sample_pval_files}/*.bed")
    
    sample_wise_pvals \
        | annotateLD

    sample_wise_pvals
        | filterTestedVariants
        | (extract_context & mutationRates & annotate_with_phenotypes)


}
