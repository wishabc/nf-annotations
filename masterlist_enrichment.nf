#!/usr/bin/env nextflow
include { extract_from_anndata } from './helpers'
params.conda = "$moduleDir/environment.yml"


process calc_prop_accessibility {
    conda params.conda
    publishDir "${params.outdir}"
    label "high_mem"
    
    input:
        tuple val(id), path(binary_matrix), path(sample_names), path(masterlist_file)

    output:
        path name
    
    script:
    name = "proportion_accessibility.tsv"
    """
    python $moduleDir/bin/calc_prop_accessibility.py \
        ${binary_matrix} \
        ${masterlist_file} \
        ${name} \
        --samples_weights ${params.sample_weights}
    """
}

process motif_hits_intersect {
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/motif_hits", pattern: "${motif_id}.hits.bed"

    input:
        tuple val(motif_id), path(moods_file), path(masterlist_file)

    output:
        tuple val(motif_id), path(indicator_file)

    script:
    indicator_file = "${motif_id}.hits.bed"
    """
    zcat ${moods_file} \
        | bedmap --indicator --sweep-all \
        --fraction-map 1 <(grep -v '#' ${masterlist_file}) - > ${indicator_file}
    """
}

process motif_enrichment_z_score {
    conda params.conda
    tag "${motif_id}:${matrix_type}"
    publishDir "${params.outdir}/motif_enrichment/${matrix_type}"
    label "med_mem"

    input:
        tuple val(motif_id), path(indicator_file), val(matrix_type), path(binary_matrix), path(sample_names), path(accessibility_proportion)
    
    output:
        tuple val(matrix_type), val(motif_id), path(name)
    
    script:
    name = "${motif_id}.${matrix_type}.z_score.tsv"
    """
    python $moduleDir/bin/subsample_proportion.py \
        ${motif_id} \
        ${indicator_file} \
        ${name} \
        ${binary_matrix} \
        ${accessibility_proportion} \
        --sample_names ${sample_names} \
        --n_bins ${params.matching_bins}
    """
}

workflow motifEnrichment {
    take:
        data
    main: 
        out = data
            | motif_enrichment_z_score
            | collectFile(
                storeDir: params.outdir,
                skip: 1,
                sort: true,
                keepHeader: true
            ) { it -> [ "${it[0]}.z_score_stats.tsv", it[2].text ] }
    emit:
        out
}

workflow fromMatrix {
    take:
        matrices // prefix, matrix, names
    main:
        accessibility = Channel.fromPath("${params.template_run}/proportion_accessibility.tsv", checkIfExists: true)
        Channel.fromPath("${params.template_run}/motif_hits/*.hits.bed")
            | map(it -> tuple(it.name.replaceAll('.hits.bed', ''), it)) // motif_id, motif_hits
            | combine(matrices) // motif_id, motif_hits, prefix, matrix, names
            | combine(accessibility) // motif_id, motif_hits, prefix, matrix, names, accessibility
            | motifEnrichment
}

workflow fromMatricesList {
    Channel.fromPath(params.matrices_list)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.matrix_name, file(row.matrix), file(row.sample_names)))
        | fromMatrix
}


workflow fromBinaryMatrix {
    matrices = Channel.of(tuple("DHS_binary", file(params.index_anndata), file(params.peaks_mask))) // kinda hotfixes
        | extract_from_anndata // prefix, matrix, sample_names, dhs_names
    
    prop_accessibility = matrices
        | calc_prop_accessibility
    
    Channel.fromPath("${params.moods_scans_dir}/*") // result of nf-genotyping scan_motifs pipeline
        | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))
        | combine(matrices.map(it -> it[3]))
        | motif_hits_intersect // motif_id, indicator
        | combine(matrices)
        | combine(prop_accessibility)
        | motifEnrichment
}


process sample_matching_bg {

    conda '/home/afathul/miniconda3/envs/r-kernel'
    tag "${iter}"

    input:
        tuple val(iter), path(masterlist_file)
    
    output:
        tuple val(iter), path(name)

    script:
    name = "masterlist.sampled.${iter}.bed"
    """
    grep -v '#' ${masterlist_file} \
        | cut -f 1-3 > tmp.bed

    Rscript $moduleDir/bin/motif_enrichment/delta_svm_match_bg.R \
        tmp.bed \
        ${name}
    """
}

workflow getRegionsSamplingPool {
    masterlist = Channel.fromPath(params.masterlist_file)
    
    sampled_bg = Channel.of(1..100)
        | combine(masterlist)
        | sample_matching_bg
        | map(it -> it[1])
        | collectFile(
            sort: true,
            name: 'sampled_regions_pool.bed',
        )
        | combine(masterlist)
        //| merge_annotations
}



// DEFUNC rn
// process choose_bins {
//     conda params.conda

//     input:
//         path sampled_regions_pool

//     output:
//         path name

//     script:
//     name = "sampled_regions.with_bins.bed"
//     """
//     Rscript $moduleDir/bin/motif_enrichment/genome_background.R \
//         tmp.bed \
//         ${name}
//     """
// }

// process overlap_and_sample {
//         conda params.conda

//     input:
//         tuple val(motif_id), path(motif_indicator), path(sampled_regions_pool), path(masterlist)

//     output:
//         path name

//     script:
//     name = "${motif_id}.sampled_regions.bed"
//     """
//     python3 $moduleDir/bin/sample_regions.py \
//         ${sampled_regions_pool} \
//         ${masterlist} \
//         ${motif_indicator} \
//         ${name}
//     """
// }

// workflow matchBackground {

//     Channel.fromPath("${params.template_run}/motif_hits/*.hits.bed")
//         | map(it -> tuple(it.name.replaceAll('.hits.bed', ''), it))
//         | filter { it[0] == "M02739_2.00" }
//         | combine(
//             masterlist
//         ) // motif_id, motif_hits, masterlist
//         | motif_hits_intersect
//         | combine(sampled_bg)
//         | 
        
//         | generate_bed

// }