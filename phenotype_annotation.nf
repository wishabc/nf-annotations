#!/usr/bin/env nextflow
include { filterUniqPvals } from "./motif_enrichment"

params.conda = "$moduleDir/environment.yml"

params.phenotypes_data = "/home/sabramov/phenotypes_data"

is_baseline = false
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
    echo "CHR\tBP\tSNP\tCM\t${annotation.simpleName}" | gzip > ${name}
    zcat ${baseannotation} \
        | awk -v OFS='\t' -F'\t' '(NR > 1) { print \$1,\$2-1,\$2,\$3,\$4 }'\
        | bedtools intersect -wa -c -a stdin -b ${annotation} \
        | awk -v OFS='\t' '{ print \$1,\$3,\$4,\$5}' | gzip >> ${name}
    """
}

// TODO wrap in apptainer
process calc_ld {
    publishDir "${outdir}/l2_logs", pattern: "${name}.log", enabled: !is_baseline
    publishDir "${outdir}", pattern: "${name}.l2.ldscore.gz"
    publishDir "${outdir}", pattern: "${name}.l2.M*"
    publishDir "${outdir}/l2", pattern: "${annotation_file}", enabled: !is_baseline
    tag "chr${chrom}:${annotation_file.simpleName}"
    scratch true
    conda params.ldsc_conda

    input:
        tuple val(chrom), path(annotation_file)
    
    output:
        tuple val(annotation_file.simpleName), path("${name}*"), path(annotation_file)
    
    script:
    if (is_baseline) {
        outdir = file(params.base_ann_path).parent
        name = "${annotation_file.simpleName}.${chrom}"
    } else {
        outdir = "${params.outdir}/${annotation_file.simpleName}/"
        name = "l2/${annotation_file.simpleName}.${chrom}"
    }
    """
    mkdir l2
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
    is_baseline = true
    calc_ld(data)
}

workflow {
    custom_annotations = Channel.fromPath(
        "${params.annotations_dir}/*"
    )
    data = Channel.of(1..22).combine(custom_annotations)
    anns = make_ldsc_annotation(data) 
    lds = calc_ld(anns)
    ldsc_data = lds.map(it -> tuple(it[0], [it[1], it[2]].flatten()))
        .groupTuple(size: 22)
        .map(
            it -> tuple(it[0], it[1].flatten())
        )
    LDSC(ldsc_data)
}

workflow annotateWithPheno {
    out =  Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
        | filterUniqPvals
        | annotate_with_phenotypes
}
