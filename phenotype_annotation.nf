#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"


// Annotates with pheWAS, clinvar, finemapping, grasp, ebi-gwas phenotypes
process annotate_with_phenotypes {
    conda params.conda
    tag "${sample_id}"
    publishDir "${params.outdir}/phenotypes"

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


// TODO wrap in apptainer
process find_ld {
    publishDir "${params.outdir}/l2_logs", pattern: "${name}.log"
    publishDir "${params.outdir}/l2", pattern: "${name}.l2.ldscore.gz"
    publishDir "${params.outdir}/l2", pattern: "${name}.l2.M*"
    tag "chr${chrom}"
    scratch true
    conda params.ldsc_conda

    input:
        val chrom
    
    output:
        tuple path("${name}*"), path(ann_path)
    
    script:
    prefix = file(params.ann_path).name
    name = "result/${prefix}${chrom}"
    ann_path = "${params.ann_path}${chrom}.annot.gz"
    """
    mkdir result
    ${ldsc_scripts_path}/ldsc.py \
        --print-snps ${params.ukbb_snps} \
        --ld-wind-cm 1.0 \
        --out ${name} \
        --bfile ${params.gtfiles}${chrom} \
        --annot ${ann_path} \
        --l2
    """
}
// TODO wrap in apptainer
process run_ldsc {
    conda params.ldsc_conda
    publishDir "${params.outdir}/ldsc", pattern: "${name}.results"
    publishDir "${params.outdir}/ldsc_logs", pattern: "${name}.logs"
    publishDir "${params.outdir}/ldsc_logs", pattern: "${name}.part_delete"
    tag "${phen_name}"
    scratch true

    input:
        tuple val(phen_id), val(phen_name), path(sumstats_file)
        path "ld_files/*"
    
    output:
        tuple val(phen_id), val(phen_name), path("${name}*")

    script:
    name = "${phen_id}"
    ld_prefix = file(params.ann_path).name
    """
    ${ldsc_scripts_path}/ldsc.py \
        --h2 ${sumstats_file} \
        --ref-ld-chr ${ld_prefix} \
        --frqfile-chr ${params.frqfiles} \
        --w-ld-chr ${params.weights} \
        --overlap-annot \
        --print-coefficients \
        --print-delete-vals \
        --out ${name}
    """
}

workflow annotateWithPheno {
    
    pvals = Channel.fromPath("${params.pval_file_dir}/*.bed")
        .map(it -> file(it))
    annotate_with_phenotypes(pvals)
}


workflow LDSC {
    take:
        ld_data
    main:
        phens = Channel.fromPath(params.phenotypes_meta)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.phen_id, row.phen_name, file(row.sumstats_file)))
        
        run_ldsc(phens, ld_data)
    emit:
        run_ldsc.out
}

workflow regressionOnly {
    // Remove second part from concat once M logs being processed
    ld_data = Channel.fromPath("${params.ann_path}*")
        .concat(
            Channel.fromPath("/net/seq/data2/projects/sabramov/LDSC/test_ldsc/output/l2/result/baselineLD.*"),
            Channel.fromPath("/net/seq/data2/projects/sabramov/LDSC/test_ldsc/output/l2_logs/result/baselineLD.*.M*")
        ).collect()
    LDSC(ld_data)
}

workflow {
    chroms = Channel.of(1..22)
    LDSC(find_ld(chroms).collect())
}
