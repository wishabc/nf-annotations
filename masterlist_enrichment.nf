
params.conda = "/home/afathul/miniconda3/envs/motif_enrichment"

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
        tuple val(motif_id), path(indicator_file), val(matrix_type), path(binary_matrix), path(accessibility_proportion)
    
    output:
        tuple val(motif_id), val(matrix_type), path(name)
    
    script:
    name = "${motif_id}.${matrix_type}.z_score.tsv"
    """
    python $moduleDir/bin/subsample_proportion.py \
        ${motif_id} \
        ${indicator_file} \
        ${name} \
        ${binary_matrix} \
        ${params.samples_file} \
        ${accessibility_proportion}
    """
}

workflow motifEnrichment {
    take:
        matrices
        prop_accessibility
    main: 
        out = Channel.fromPath("${params.moods_scans_dir}/*") // result of nf-genotyping scan_motifs pipeline
            | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))
            | motif_hits_intersect // motif_id, indicator
            | combine(matrices)
            | combine(prop_accessibility)
            | motif_enrichment_z_score
            | map(it -> it[2])
            | collectFile(name: 'all.z_score.stats.tsv',
                storeDir: "${params.outdir}",
                skip: 1,
                sort: true,
                keepHeader: true
            )
    emit:
        out
}

workflow {
    // matrices = Channel.fromPath(params.matrices_list)
    //     | splitCsv(header:true, sep:'\t')
    //     | map(row -> tuple(row.matrix_id, file(row.matrix)))
    matrices = Channel.of(
        //tuple("NMF", file(params.binary_nmf)),
        tuple("DHS_Binary", file(params.binary_matrix))
    )
    acc = calc_prop_accessibility()

    motifEnrichment(matrices, acc)

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