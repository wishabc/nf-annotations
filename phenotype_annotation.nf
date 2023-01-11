#!/usr/bin/env nextflow
include { filterUniqPvals } from "./motif_enrichment"

params.conda = "$moduleDir/environment.yml"

params.phenotypes_data = "/home/sabramov/phenotypes_data"
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
    echo "CHR\tBP\tSNP\tCM\t${annotation.simpleName}" > ${name}
    cat ${baseannotation} \
        | awk -v OFS='\t' -F'\t' '{ print \$1,\$2-1,\$2,\$3,\$4 }'\
        | bedtools intersect -wa -c -a - -b ${annotation} \
        | awk -v OFS='\t' '{ print \$1,\$3,\$4,\$5}' >> ${name}
    """
}

// TODO wrap in apptainer
process calc_ld {
    publishDir "${params.outdir}/l2_logs", pattern: "${name}.log"
    publishDir "${params.outdir}/${annotation_file.simpleName}/", pattern: "${name}.l2.ldscore.gz"
    publishDir "${params.outdir}/${annotation_file.simpleName}/", pattern: "${name}.l2.M*"
    publishDir "${params.outdir}/${annotation_file.simpleName}/l2", pattern: "${annotation_file}"
    tag "chr${chrom}:${annotation_file.simpleName}"
    scratch true
    conda params.ldsc_conda

    input:
        tuple val(chrom), path(annotation_file)
    
    output:
        tuple val(annotation_file.simpleName), path("${name}*"), path(annotation_file)
    
    script:
    name = "l2/${annotation_file.simpleName}.${chrom}"
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
    tag "${prefix}:${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file), val(prefix), path(ld_files)
    
    output:
        tuple val(phen_id), path("${name}*")

    script:
    name = "${prefix}.${phen_id}"
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


workflow LDSC {
    take:
        ld_data
    main:
        phens = Channel.fromPath(params.phenotypes_meta)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.phen_id, file(row.sumstats_file)))
        d = phens.combine(ld_data)
        run_ldsc(d)
    emit:
        run_ldsc.out
}

workflow calcBaseline {
    data = Channel.of(1..22).map(
        it -> tuple(it, file("${params.base_ann_path}${it}.annot.gz", checkIfExists: true))
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

workflow annotateWithPheno {
    out = filterUniqPvals() | annotate_with_phenotypes
}
