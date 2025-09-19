#!/usr/bin/env nextflow
include { matricesListFromMeta; splitMatrices } from './helpers'
params.conda = "$moduleDir/environment.yml"


process overlap_and_sample {
    conda params.conda
    tag "${annotation_name}"
    publishDir "${params.outdir}/motif_enrichment/per_category_samples/"
    label "med_mem"

    input:
        tuple val(annotation_name), path(annotation), path(annotation_coordinates), path(sampled_regions_pool), path(masterlist)

    output:
        tuple val(annotation_name), path(name)

    script:
    name = "${annotation_name}.sampled_regions.bed"
    """
    python3 $moduleDir/bin/motif_enrichment/sample_regions.py \
        ${masterlist} \
        ${sampled_regions_pool} \
        ${annotation} \
        ${annotation_coordinates} \
        ${name}
    """
}

process motif_hits_intersect {
    tag "${motif_id}"
    conda params.conda
    publishDir { prefix == "keep" ? "${params.outdir}/motif_hits" : null }

    input:
        tuple val(prefix), path(bed_file), val(motif_id), path(moods_file)

    output:
        tuple val(prefix), val(motif_id), path(bed_file), path(name)

    script:
    name = "${motif_id}.motif_hits_mask.txt"
    """
    zcat ${moods_file} \
        | bedmap --indicator --sweep-all \
        --fraction-map 1 <(grep -v '#' ${bed_file}) - > ${name}
    
    """
}

process count_entries {
    conda params.conda
    tag "${prefix}"
    publishDir "${params.outdir}/motif_enrichment/"
    label "high_mem"

    input:
        tuple val(prefix), val(motif_id), path(bed_file), path(motif_hits_mask)

    output:
        tuple val(prefix), val(motif_id), path(name)

    script:
    name = "${prefix}.${motif_id}.stats.tsv"
    """
    python3 $moduleDir/bin/motif_enrichment/count_entries.py \
        ${bed_file} \
        ${motif_hits_mask} \
        ${motif_id} \
        ${prefix} \
        ${name}
    """
}


process calc_pvals {
    conda params.conda
    tag "${prefix}"
    publishDir "${params.outdir}/motif_enrichment/"
    label "high_mem"

    input:
        tuple val(prefix), path(sampling_results)

    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.pvals.tsv"
    """
    python3 $moduleDir/bin/motif_enrichment/calc_pvals.py \
        ${sampling_results} \
        ${name}
    """
}


workflow annotationEnrichment {
    take:
        annotations // matrix_name, prefix, annotation_mask, dhs_coordinates
    main:
        motifs_meta = Channel.fromPath(params.motifs_metadata)
            | splitCsv(header: true, sep: "\t")
            | map(row -> tuple(
                row.motif_id, 
                row.moods_scans_ref // result of nf-genotyping scan_motifs pipeline
                )
            )
            //| filter { it[0] in ["M02739_2.00"] }

        sampled = Channel.fromPath("${params.template_run}/motif_enrichment/sampled_regions_pool.parquet")
         
        annotated_masterlist = Channel.fromPath("${params.template_run}/motif_enrichment/index.annotated.bed")

        result = annotations
            | map(it -> tuple(it[1], it[2], it[3]))
            | combine(sampled)
            | combine(annotated_masterlist)
            | overlap_and_sample
            | combine(motifs_meta)
            | motif_hits_intersect
            | count_entries // prefix, stats
            | combine(
                annotations.map(it -> tuple(it[1], it[0])),
                by: 0
            ) // prefix, stats, matrix_name
            | collectFile(
                skip: 1,
                keepHeader: true
            ) {
                [
                    "${it[3]}.sampled.bed",
                    it[2].text
                ]
            }
            | map(it -> tuple(it.name.replaceAll('.sampled.bed', ''), it))
            | calc_pvals
    emit:
        result
}


workflow fromMatricesList {
    matricesListFromMeta()
        | splitMatrices
        | annotationEnrichment
}



workflow motifIndexOverlap {
    motifs_meta = Channel.fromPath(params.motifs_metadata)
        | splitCsv(header: true, sep: "\t")
        | map(row -> tuple(
            row.motif_id, 
            row.moods_scans_ref // result of nf-genotyping scan_motifs pipeline
            )
        )
        | map(
            it -> tuple("keep", file(params.index_bed), it[0], it[1])
        )
        | motif_hits_intersect
}
