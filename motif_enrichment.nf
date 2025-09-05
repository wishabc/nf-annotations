#!/usr/bin/env nextflow
include { matricesListFromMeta; splitMatrices } from './helpers'
params.conda = "$moduleDir/environment.yml"


process split_masterlist_in_chunks {
    conda params.conda

    input:
        tuple val(prefix), path(masterlist_file)

    output:
        path "${prefix}*.bed"

    script:
    """
    grep -v '#' ${masterlist_file} \
        | cut -f 1-3 \
        | split \
            -l 100000 \
            --suffix-length=4 \
            -d \
            --additional-suffix=.bed \
            - \
            ${prefix}
    """
}
process sample_matching_bg {

    conda '/home/afathul/miniconda3/envs/r-kernel'
    tag "${prefix}:${iter}"

    input:
        tuple val(iter), val(prefix), path(bed_file)
    
    output:
        tuple val(iter), path(name)

    script:
    name = "${prefix}.${iter}.bed"
    """
    Rscript $moduleDir/bin/motif_enrichment/delta_svm_match_bg.R \
        ${bed_file} \
        tmp.bed
    
    grep -w -F -f ${params.nuclear_chroms} tmp.bed > ${name}
    """
}

process annotate_regions {

    conda '/home/sabramov/miniconda3/envs/super-index'
    publishDir "${params.outdir}/motif_enrichment/", pattern: "index.annotated.bed"
    tag "${prefix}"
    label "high_mem"
    scratch true

    input:
        tuple val(prefix), path(bed_file)

    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.annotated.bed"
    """
    cat ${bed_file} \
        | cut -f 1-3 \
        | sort-bed - > tmp.bed

    faidx -i nucleotide \
        -b tmp.bed \
        ${params.genome_fasta} \
        | awk -v OFS="\t" \
            'NR>1 {  \
                total=\$4+\$5+\$6+\$7+\$8; \
                cg=\$6+\$7; \
                print \$1, \$2-1, \$3, \$3-\$2, cg, cg/total; }' > annotated.bed
    
    python3 $moduleDir/bin/motif_enrichment/assign_bins.py annotated.bed ${params.length_bins_bounds} ${name}
    """
}


process to_parquet {
    conda params.conda
    tag "${prefix}"
    publishDir "${params.outdir}/motif_enrichment/"
    label "high_mem"

    input:
        tuple val(prefix), path(bed_file)

    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.parquet"
    """
    python3 $moduleDir/bin/motif_enrichment/convert_to_parquet.py ${bed_file} ${name}
    """
}

workflow getRegionsPool {
    masterlist = Channel.fromPath(params.masterlist_file)
        | map(it -> tuple("index", it))

    chunks = masterlist
        | split_masterlist_in_chunks
        | flatten()
        | map(it -> tuple(it.baseName, it)) // prefix, bed_file
    
    sampled_bg = Channel.of(1..100)
        | combine(chunks)
        | sample_matching_bg
        | mix(chunks)
        | mix(masterlist)
        | annotate_regions
        | filter { it[0] != 'index' }
        | map(it -> it[1])
        | collectFile(
            sort: true,
            name: 'sampled_regions_pool.bed',
        )
        | map(it -> tuple('sampled_regions_pool', it))
        | to_parquet
}



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
    echo 1
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

workflow fromMatricesList {
    matricesListFromMeta()
        | splitMatrices
        | annotationEnrichment
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
            | count_entries
            | combine(
                annotations.map(it -> tuple(it[1], it[0])),
                by: 0
            )
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
