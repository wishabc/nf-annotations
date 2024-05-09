
process calc_prop_accessibility {
    conda params.conda
    publishDir "${params.outdir}"
    
    output:
        path name
    
    script:
    name = "proportion_accessibility.tsv"
    """
    python $moduleDir/bin/calc_prop_accessibility.py \
        ${params.binary_matrix} \
        ${params.masterlist_file} \
        ${params.dhs_annotations} \
        ${name} \
        --samples_weights ${params.sample_weights} \
        --sampling_n 500 \
        --n 10
    """
}

process motif_hits_intersect {
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/motif_hits", pattern: "${motif_id}.hits.bed"

    input:
        tuple val(motif_id), path(moods_file)

    output:
        tuple val(motif_id), path(indicator_file)

    script:
    indicator_file = "${motif_id}.hits.bed"
    """
    zcat ${moods_file} \
        | bedmap --indicator --sweep-all \
        --fraction-map 1 ${params.masterlist_file} - > ${indicator_file}
    """
}

process motif_enrichment_z_score {
    conda params.conda
    tag "${motif_id}:${matrix_type}"
    publishDir "${params.outdir}/${matrix_type}"

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
        ${params.samples_file} \
        ${accessibility_proportion} \
        --sample_names ${sample_names}
    """
}

workflow motifEnrichment {
    take:
        data
    main: 
        out = Channel.fromPath("${params.moods_scans_dir}/*") // result of nf-genotyping scan_motifs pipeline
            | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))
            | motif_hits_intersect // motif_id, indicator
            | combine(data)
            | motif_enrichment_z_score
            | groupTuple()
            | collectFile(
                {
                    it -> [ "${it[0]}.z_score_stats.tsv", it[2].text ]
                },
                storeDir: "${params.outdir}",
                skip: 1,
                sort: true,
                keepHeader: true
            )
    emit:
        out
}


workflow categoryEnrichment {
    params.template_run = "${params.outdir}"
    accessibility = Channel.fromPath("${params.template_run}/proportion_accessibility.tsv")
    matrices = Channel.fromPath(params.matrices_list)
       | splitCsv(header:true, sep:'\t')
       | map(row -> tuple(row.matrix_name, file(row.matrix), file(row.sample_names)))
       | combine(accessibility)
       | motifEnrichment
}


workflow {
    matrices = Channel.fromPath(params.binary_matrix)
        | map(it -> tuple("DHS_Binary", it, file(params.sample_names)))
        | combine(calc_prop_accessibility())
        | motifEnrichment
}



// DEFUNC
// process tf_by_components {
//     conda params.conda
//     publishDir "${params.outdir}/plot"

//     input:
//         path(all_coefs)
    
//     output:
//         path("${prefix}*.pdf")
    
//     script:
//     prefix = "components"
//     """
//     python $moduleDir/bin/plot_tf_by_components.py \
//         ${all_coefs} \
//         ${params.metadf} \
//         ${prefix}
//     """
// }