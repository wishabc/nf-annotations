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

process make_ldsc_annotation {
    conda params.conda
    tag "chr${chrom}:${annotation.simpleName}"
    scratch true

    input:
        tuple val(chrom), path(annotation)

    output:
        tuple val(chrom), path(name)
    
    script:
    suffix = "${chrom}.annot.gz"
    baseannotation = "${params.base_ann_path}${suffix}"
    name = "${annotation.simpleName}.${suffix}"
    """
    python3 $moduleDir/bin/make_annotation.py ${baseannotation} ${annotation} ${name}
    """
}

// TODO wrap in apptainer
process calc_ld {
    publishDir "${params.outdir}/l2_logs", pattern: "${name}.log"
    publishDir "${params.outdir}/${annotation_file.simpleName}/l2", pattern: "${name}.l2.ldscore.gz"
    publishDir "${params.outdir}/${annotation_file.simpleName}/l2", pattern: "${name}.l2.M*"
    publishDir "${params.outdir}/${annotation_file.simpleName}/l2/result", pattern: "${annotation_file}"
    tag "chr${chrom}:${annotation.simpleName}"
    scratch true
    conda params.ldsc_conda

    input:
        tuple val(chrom), path(annotation_file)
    
    output:
        tuple val(annotation_file.simpleName), path("${name}*"), path(annotation_file)
    
    script:
    name = "result/${annotation_file.simpleName}.${chrom}"
    """
    mkdir result
    # Check if --print-snps parameter is needed
    ${params.ldsc_scripts_path}/ldsc.py \
        --print-snps ${params.tested_snps} \
        --ld-wind-cm 1.0 \
        --out ${name} \
        --bfile ${params.gtfiles}${chrom} \
        --annot ${annotation_file} \
        --l2
    """
}
// TODO wrap in apptainer
process run_ldsc {
    conda params.ldsc_conda
    publishDir "${params.outdir}/${prefix}/ldsc", pattern: "${name}.results"
    publishDir "${params.outdir}/ldsc_logs", pattern: "${name}.logs"
    publishDir "${params.outdir}/ldsc_logs", pattern: "${name}.part_delete"
    tag "${prefix}:${phen_name}"
    scratch true
    cache false

    input:
        tuple val(phen_id), val(phen_name), path(sumstats_file), val(prefix), path(ld_files)
    
    output:
        tuple val(phen_id), val(phen_name), path("${name}*")

    script:
    name = "${prefix}_${phen_id}"
    pref = "${prefix}."
    """
    ${params.ldsc_scripts_path}/ldsc.py \
        --h2 ${sumstats_file} \
        --ref-ld-chr ${params.base_ann_path},${pref} \
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
        d = phens.combine(ld_data)
        run_ldsc(d)
    emit:
        run_ldsc.out
}


workflow calcBaseline {
    data = Channel.of(1..22).map(
        it -> tuple(it, file("${params.base_ann_path}/${it}.annot.gz"))
    )
    calc_ld(data)
}
workflow {
    custom_annotations = Channel.fromPath(
        "${params.annotations_dir}/*"
    )
    data = Channel.of(1..22).combine(custom_annotations)
    lds = make_ldsc_annotation(data) | calc_ld
    ldsc_data = lds.map(it -> tuple(it[0], [it[1], it[2]].flatten()))
        .groupTuple(size: 22)
        .map(
            it -> tuple(it[0], it[1].flatten())
        )
    LDSC(ldsc_data)
}
