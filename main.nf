include { fromMatrix as ldscFromMatrix } from './ldsc' addParams(by_cell_type: true)
include { fromMatrix as motifEnrichmentFromMatrix; extract_from_anndata } from './masterlist_enrichment'


process find_top_samples {

    conda params.conda
    tag "${prefix}"
    publishDir "${params.outdir}/top_samples", pattern: "${name}"
    publishDir "${params.outdir}/top_samples", pattern: "${res}"

    input:
        tuple val(prefix), path(W_matrix), path(samples_order), path(anndata)

    output:
        tuple path("*.*.component_${prefix}.bw"), path(name), path(res)

    script:
    name = "${prefix}.top_samples.tsv"
    res = "${prefix}.density_tracks_meta.tsv"
    """
    python3 $moduleDir/bin/find_top_samples.py \
        ${W_matrix} \
        ${samples_order} \
        ${anndata} \
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
        tuple val(prefix), path(H_matrix), path(dhs_coordinates)

    output:
        tuple val(clean_prefix), path(clean_comps_matrix), path(clean_comp_order), path(dhs_coordinates), emit: pure
        tuple val(mixing_prefix), path(mixing_matrix), path(mixing_comp_order), path(dhs_coordinates), emit: mixing

    script:
    clean_prefix = "${prefix}.pure.50pr"
    clean_comps_matrix = "${clean_prefix}.npy"
    clean_comp_order = "${clean_prefix}.order.txt"
    
    
    mixing_prefix = "${prefix}.mixing.80pr"
    mixing_matrix = "${mixing_prefix}.npy"
    mixing_comp_order = "${mixing_prefix}.order.txt"
    """
    python3 $moduleDir/bin/prepare_mixings_data.py \
        ${H_matrix} \
        ${prefix}
    """
}

process craft_configs {
    conda params.conda
    publishDir "${params.outdir}"

    output:
        tuple path("*.config.ini"), path("*.matrix_meta.tsv")

    script:
    """
    python3 $moduleDir/bin/craft_configs.py \
        ${params.nmf_metadata} \
        ${params.outdir}
    """

}

workflow {
    input_data = Channel.fromPath(params.nmf_metadata)
        | splitCsv(header: true, sep: "\t")
        | map(row -> tuple(
            row.prefix,
            file(row.anndata_path),
            file(row.peaks_mask),
            file(row.W),
            file(row.H),
            row?.peaks_weights,
            row?.samples_weights,
            )
        )

    nmf_data = input_data
        | map(it -> tuple(it[0], it[1], it[2]))
        | unique { it[0] }
        | extract_from_anndata // prefix, binary, samples_order, masterlist
        | combine(input_data, by: 0) // prefix, binary, samples_order, masterlist, anndata, peaks_mask, W, H, peaks_weights, samples_weights
        | map(it -> tuple(it[0], it[6], it[7], it[2], it[3], it[8], it[9]))  // prefix, W, H,  samples_order, masterlist, peaks_weights, samples_weights
    
    // Top samples tracks
    nmf_data
        | map(it -> tuple(it[0], it[1], it[3]))
        | join(input_data.map(it -> tuple(it[0], it[1])))
        | find_top_samples
        | map(it -> it[0])
        | flatten()
        | combine(nmf_data.map(it -> it[0]))
        | map(it -> tuple(it[0].simpleName, it[1], it[0]))
        | groupTuple(by: [0, 1])
        | top_samples_track

    // Mixings
    mixing_data = nmf_data 
        | map(it -> tuple(it[0], it[2], it[4]))
        | prepare_mixings_data // prefix, matrix, names, dhs
    
    if (!file("${params.template_run}/proportion_accessibility.tsv").exists()) {
        error "No accessibility file found at ${params.template_run}/proportion_accessibility.tsv; please run masterlist_enrichment:fromBinaryMatrix first or specify template_run folder. Once per binary matrix."
    }
    
    data = mixing_data.pure
        | mix(mixing_data.mixing) // prefix, matrix, names, dhs
    
    data | ldscFromMatrix
       // | (motifEnrichmentFromMatrix & ldscFromMatrix) // ldsc. ALWAYS uses by_cell_type version if run from here.
    
    craft_configs()
}
