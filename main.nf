include { fromMatrix as ldscFromMatrix } from './ldsc' addParams(by_cell_type: true)
include { fromMatrix as motifEnrichmentFromMatrix } from './masterlist_enrichment'

process find_top_samples {

    conda params.conda
    tag "${prefix}"
    publishDir "${params.outdir}", pattern: "${name}"
    publishDir "${params.outdir}", pattern: "${res}"

    input:
        tuple val(prefix), path(W_matrix), path(samples_order)

    output:
        tuple path("*.*.component_${prefix}.bw"), path(name), path(res)

    script:
    name = "${prefix}.top_samples.tsv"
    res = "${prefix}.density_tracks_meta.tsv"
    """
    python3 $moduleDir/bin/find_top_samples.py \
        ${W_matrix} \
        ${samples_order} \
        ${params.samples_file} \
        ${prefix} \
        ${params.top_count} \
        ${params.outdir}/top_samples/${prefix} \
    """
}


process top_samples_track {

    scratch true
    conda params.conda
    tag "${prefix}:${component}"
    publishDir "${params.outdir}/top_samples/${prefix}"

    input:
        tuple val(component), val(prefix), path(density_bw, stageAs: "?/*")
    
    output:
        tuple val(prefix), path(name), path(bg)
    
    script:
    name = "${prefix}.${component}.top_samples.bw"
    bg = "${prefix}.${component}.top_samples.bg"
    """
    wiggletools write_bg ${bg} mean ${density_bw}
    bedGraphToBigWig "${bg}" "${params.chrom_sizes}" "${name}"
    """
}


process prepare_mixings_data {
    conda params.conda
    publishDir "${params.outdir}/mixing/${prefix}"
    tag "${prefix}"

    input:
        tuple val(prefix), path(H_matrix)

    output:
        tuple val(clean_prefix), path(clean_comps_matrix), path(clean_comp_order), emit: pure
        tuple val(mixing_prefix), path(mixing_matrix), path(mixing_comp_order), emit: mixing

    script:
    clean_prefix = "${prefix}.pure"
    clean_comps_matrix = "${clean_prefix}.50pr.npy"
    clean_comp_order = "${clean_prefix}.50pr.order.txt"
    
    
    mixing_prefix = "${prefix}.mixing"
    mixing_matrix = "${mixing_prefix}.80pr.npy"
    mixing_comp_order = "${mixing_prefix}.80pr.order.txt"
    """
    python3 $moduleDir/bin/prepare_mixings_data.py \
        ${H_matrix} \
        ${prefix}
    """
}

process craft_config {
    conda params.conda
    publishDir "${params.outdir}"
    tag "${prefix}"

    input:
        tuple val(prefix), path(samples_order)

    output:
        tuple val(prefix), path("${prefix}.ini")

    script:
    """
    python3 $moduleDir/bin/craft_config.py \
        ${samples_order} \
        ${prefix} \
        ${params.nmf_metadata}
    """

}

workflow {
    input_data = Channel.fromPath(params.nmf_metadata)
        | splitCsv(header: true, sep: "\t")
        | map(row -> tuple(row.prefix, file(row.W), file(row.H), file(row.samples_order)))
    
    // Top samples tracks
    input_data
        | map(it -> tuple(it[0], it[1], it[3]))
        | find_top_samples
        | map(it -> it[0])
        | flatten()
        | combine(input_data)
        | map(it -> tuple(it[0].simpleName, it[1], it[0]))
        | groupTuple(by: [0, 1])
        | top_samples_track

    // Mixings
    mixing_data = input_data 
        | map(it -> tuple(it[0], it[2]))
        | prepare_mixings_data
    
    if (!file("${params.template_run}/proportion_accessibility.tsv").exists()) {
        error "No accessibility file found at ${params.template_run}/proportion_accessibility.tsv; please run masterlist_enrichment:fromBinaryMatrix first or specify template_run folder. Once per binary matrix."
    }
    
    mixing_data.pure
        | mix(mixing_data.mixing)
        | (motifEnrichmentFromMatrix & ldscFromMatrix) // ldsc. ALWAYS uses by_cell_type version if run from here.
    
    //craft_config()
}
