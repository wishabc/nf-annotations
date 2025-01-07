
process extract_from_anndata {
    conda params.conda
    label "high_mem"

    input:
        tuple val(prefix), path(index_anndata), path(peaks_mask)
    
    output:
        tuple val(prefix), path(name), path(sample_names), path(masterlist_file)
    
    script:
    name = "binary_matrix.npy"
    sample_names = "sample_names.txt"
    masterlist_file = "masterlist.bed"
    """
    python $moduleDir/bin/extract_from_anndata.py \
        ${index_anndata} \
        ${peaks_mask} \
        ${name} \
        ${sample_names} \
        ${masterlist_file} \
    """
}


process calc_prop_accessibility {
    conda params.conda
    publishDir "${params.outdir}"
    label "high_mem"
    
    input:
        tuple val(id), path(binary_matrix), path(sample_names), path(masterlist_file)

    output:
        tuple val(id), path(binary_matrix), path(sample_names), path(name)
    
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
        --fraction-map 1 ${masterlist_file} - > ${indicator_file}
    """
}

process motif_enrichment_z_score {
    conda params.conda
    tag "${motif_id}:${matrix_type}"
    publishDir "${params.outdir}/motif_enrichment/${matrix_type}"
    label "high_mem"

    input:
        tuple val(motif_id), path(indicator_file), val(matrix_type), path(binary_matrix), path(sample_names), path(dhs_masterlist), path(accessibility_proportion)
    
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
        matrices // prefix, matrix, names, dhs_names
    main:
        accessibility = Channel.fromPath("${params.template_run}/proportion_accessibility.tsv", checkIfExists: true)
        Channel.fromPath("${params.template_run}/motif_hits/*.hits.bed")
            | map(it -> tuple(it.name.replaceAll('.hits.bed', ''), it)) // motif_id, motif_hits
            | combine(matrices) // motif_id, motif_hits, prefix, matrix, names,         matrices // prefix, matrix, names, dhs_names
            | combine(accessibility) // motif_id, motif_hits, prefix, matrix, names, dhs_names, accessibility
            | motifEnrichment
}

workflow categoryEnrichment {
    Channel.fromPath(params.matrices_list)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.matrix_name, file(row.matrix), file(row.sample_names)))
        | fromMatrix
}


workflow fromBinaryMatrix {
    matrices = Channel.of(tuple("DHS_binary", file(params.index_anndata), file(params.peaks_mask))) // kinda hotfixes
        | extract_from_anndata
    
    prop_accessibility = matrices
        | calc_prop_accessibility
    
    Channel.fromPath("${params.moods_scans_dir}/*") // result of nf-genotyping scan_motifs pipeline
        | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))
        | combine(matrices.map(it -> it[3]))
        | motif_hits_intersect // motif_id, indicator
        | combine(prop_accessibility)
        | motifEnrichment
}



// TESTING //
process generate_bed {
    publishDir "${params.outdir}/matched_bg"

    tag "${motif_id}:${iter}"
    // scratch true
    conda params.conda

    input:
        tuple val(motif_id), path(indicator_file), val(iter)

    output:
        tuple val(motif_id), val(iter), path(name)

    script:
    name = "${motif_id}.${iter}.bed"
    """
    paste ${params.masterlist_file} ${indicator_file} \
        | awk \
            -v OFS='\t' \
            -F'\t' \
            '\$NF == 1' \
        | cut -f1-3 > tmp.bed
    
    Rscript $moduleDir/bin/genome_background.R \
        tmp.bed \
        ${name}
    """
}


workflow matchBackground {
    Channel.fromPath("${params.template_run}/motif_hits/*.hits.bed")
        | map(it -> tuple(it.name.replaceAll('.hits.bed', ''), it))
        | filter { it[0] == "M02739_2.00" }
        | combine(
            Channel.of(1..100)
        ) // motif_id, indicator, iter
        
        | generate_bed

}