include { readMotifsList } from "./motif_enrichment"


process cut_matrix {
    conda params.conda
    tag "samples ${interval}"
    scratch true

    input:
        val sample_id

    output:
        tuple val(interval), path(name)

    script:
    end = sample_id + params.step - 1
    interval = "${sample_id}-${end}"
    name = "${interval}.cut_matrix.npy"
    """
    zcat ${params.binary_matrix} | cut -f${interval} > tmp.txt
    python3 $moduleDir/bin/convert_to_numpy.py tmp.txt ${name}
    """
}

process motif_hits_intersect {
    publishDir "${params.outdir}/counts", pattern: "${counts_file}"
    tag "${motif_id}"
    conda params.conda

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

process calc_index_motif_enrichment {
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/motif_stats_chunks"
    scratch true

    input:
        tuple val(motif_id), path(indicator_file), val(chunk_n), path(binary_matrix_chunk)
    
    output:
        path name

    script:
    name = "${motif_id}.${chunk_n}.enrichment.tsv"
    """
    python3 $moduleDir/bin/index_motif_enrichment.py  \
        ${binary_matrix_chunk} ${indicator_file} ${motif_id} \
        ${params.sample_names} ${chunk_n} > ${name}
    """
}


workflow indexEnrichment {
    chunks_count = file(params.sample_names).countLines().intdiv(params.step)

    c_mat = Channel.of(0..chunks_count) // 0, 1, 2, 3 ,4 , 5.. 9
        | map(it -> it * params.step + 1) // 1 201 401 601 1801
        | cut_matrix // [1, 1-200.cut_matrix.npy], [201, 201-400.cut_matrix.npy]

    moods_scans = Channel.fromPath("${moods_scans_dir}/*") // motif_id, moods_path
        | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))
        | motif_hits_intersect // motif_id, indicator_file
        | combine(c_mat) // motif_id, indicator_file, chunk_n, binary_matrix_chunk
        | calc_index_motif_enrichment // enrichment_file
        | collectFile(storeDir: "$launchDir/${params.outdir}", name: "enrichments.tsv")
}

workflow {
    indexEnrichment()
}