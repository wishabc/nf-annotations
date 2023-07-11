include { readMotifsList } from "./motif_enrichment"


process cut_matrix {
    conda params.conda
    tag "samples ${interval}"
    scratch true

    input:
        tuple val(sample_id), path(binary_matrix)

    output:
        tuple val(interval), path(name)

    script:
    end = sample_id + params.step - 1
    interval = "${sample_id}-${end}"
    name = "${interval}.cut_matrix.npy"
    """
    zcat ${binary_matrix} | cut -f${interval} > tmp.txt
    python3 $moduleDir/bin/convert_to_numpy.py tmp.txt ${name}
    """
}

process motif_hits_intersect {
    publishDir "${params.outdir}/counts", pattern: "${counts_file}"
    tag "${motif_id}"
    conda params.conda

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
        | combine(file(params.binary_matrix)) // [chunk_n, binary_matrix]
        | cut_matrix // [1, 1-200.cut_matrix.npy], [201, 201-400.cut_matrix.npy]

    moods_scans = readMotifsList() // motif_id, motif_path
        | combine(file(params.masterlist_file)) // motif_id, motif_path, masterlist
        | motif_hits_intersect // motif_id, indicator_file
        | combine(c_mat) // motif_id, indicator_file, chunk_n, binary_matrix_chunk
        | calc_index_motif_enrichment // enrichment_file
        | collectFile(storeDir: "$launchDir/${params.outdir}", name: "enrichments.tsv")
}

workflow {
    indexEnrichment()
}


// Workflow 
process logistic_regression {
    conda params.r_conda
    publishDir "${params.outdir}/logreg_results"

    input:
        tuple val(motif_id), path(indicator_file)
        path matrix
    
    output:
        tuple val(motif_id), path("${prefix}.metrics.tsv"), path("${prefix}.coeff.tsv")
    
    script:
    prefix = "${motif_id}"
    """
    Rscript $moduleDir/bin/motif_enruchment.R \
        ${matrix} \
        ${indicator_file} \
        ${prefix}
    """

}

workflow {
    params.r_conda = "/home/afathul/miniconda3/envs/r-kernel"
    params.samples_file = ""
    params.matrix = ""
    motifs = Channel.fromPath(params.samples_file)
		| splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.motif_id, file(row.indicator_file)))
    
    logistic_regression(motifs, params.matrix)
}