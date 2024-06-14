include { fromMatrix as ldscFromMatrix } from './ldsc' addParams(by_cell_type: true)
include { fromMatrix as motifEnrichmentFromMatrix } from './masterlist_enrichment'

process find_top_samples {

    conda params.conda
    publishDir "${params.outdir}/top_samples"

    input:
        tuple val(prefix), path(W_matrix), path(samples_order)

    output:
        path "*.*.component_${prefix}.bw"

    script:
    """
    python3 $moduleDir/bin/find_top_samples.py \
        ${W_matrix} \
        ${samples_order} \
        ${params.samples_file} \
        ${prefix} \
        ${params.top_count}
    """
}


process top_samples_track {

    scratch true
    conda params.conda
    tag "${component}"
    publishDir "${params.outdir}/top_samples"

    input:
        tuple val(component), path(density_bw, stageAs: "?/*")
    
    output:
        tuple path(name), path("${component}.top_samples.bg")
    
    script:
    name = "${component}.top_samples.bw"
    bg = "${component}.top_samples.bg"
    """
    wiggletools write_bg ${bg} mean ${density_bw}
    bedGraphToBigWig "${bg}" "${params.chrom_sizes}" "${name}"
    """
}


process prepare_mixings_data {
    conda params.conda
    publishDir "${params.outdir}/mixings"

    input:
        tuple val(prefix), path(W_matrix), path(H_matrix), path(samples_order)

    output:
        tuple val("${prefix}.clean"), path(clean_comps_matrix), path(clean_comp_order), emit: clean
        tuple val("${prefix}.mixing"), path(mixing_matrix), path(mixing_comp_order), emit: mixing

    script:
    clean_prefix = "${prefix}.clean"
    mixing_prefix = "${prefix}.mixing"
    clean_comps_matrix = "${clean_prefix}.50pr.npy"
    mixings_matrix = "${mixing_prefix}.mixings_80pr.npy"
    clean_comp_order = "${clean_prefix}.clean_50pr.order.txt"
    mixing_comp_order = "${mixing_prefix}.mixings_80pr.order.txt"
    """
    python3 $moduleDir/bin/prepare_mixings_data.py \
        ${W_matrix} \
        ${H_matrix} \
        ${samples_order} \
        ${params.samples_file} \
        ${prefix}
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
        | map(it -> tuple(it.simpleName, it))
        | groupTuple()
        | top_samples_track

    // // Mixings
    // mixing_data = prepare_mixings_data(input_data)
    
    // if !file("${params.template_run}/proportion_accessibility.tsv").exists() {
    //     error "No accessibility file found at ${params.template_run}/proportion_accessibility.tsv; please run masterlist_enrichment:fromBinaryMatrix first. Once per binary matrix."
        
    // }
    
    // mixing_data.clean
    //     | mix(mixing_data.mixing)
    //     | (motifEnrichmentFromMatrix & ldscFromMatrix) // ldsc. ALWAYS uses by_cell_type version if run from here.
}
