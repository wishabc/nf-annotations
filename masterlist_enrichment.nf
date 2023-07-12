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
    tag "${motif_id}"
    publishDir "${params.outdir}/metrics", pattern: "${prefix}.metrics.tsv"
    publishDir "${params.outdir}/coeffs", pattern: "${prefix}.coeff.tsv"

    input:
        tuple val(motif_id), path(indicator_file), path(matrix)
    
    output:
        tuple val(motif_id), path("${prefix}.metrics.tsv"), path("${prefix}.coeff.tsv")
    
    script:
    prefix = "${motif_id}"
    """
    Rscript $moduleDir/bin/motif_enrichment.R \
        ${matrix} \
        ${indicator_file} \
        ${prefix}
    """
}

workflow logisticRegression {
    params.r_conda = "/home/afathul/miniconda3/envs/r-kernel"
    params.masterlist_file = "/net/seq/data2/projects/afathul/motif_enhancement/masterlist.filtered.bed"
    params.outdir = "/net/seq/data2/projects/afathul/motif_enhancement"
    params.matrix = "/net/seq/data2/projects/afathul/motif_enhancement/bin_new_unweight_full.16.H.npy"
    
    coeffs = Channel.fromPath("${params.moods_scans_dir}/*")
        | map (it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it, params.masterlist_file))
        | motif_hits_intersect
        | combine(Channel.fromPath(params.matrix))
	    | logistic_regression
    
    coeffs | map(it -> it[1])
        | collectFile(name: 'all.metrics.tsv',
            storeDir: "${params.outdir}",
            skip: 1,
            sort: true,
            keepHeader: true)

    coeffs | map(it -> it[2])
        | collectFile(name: 'all.coeff.tsv',
            storeDir: "${params.outdir}",
            skip: 1,
            sort: true,
            keepHeader: true)
}