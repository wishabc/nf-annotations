include { splitMatrices; matricesListFromMeta } from './helpers'


process overlap_annotation {
    conda params.conda
    publishDir "${params.outdir}/finemapping/${matrix_name}"
    tag "${prefix}"
    
    input:
        tuple val(matrix_name), val(prefix), path(dhs_mask), path(dhs_coordinates), path(variants)
    
    output:
        tuple val(matrix_name), path(name)
    
    script:
    name = "${prefix}.overlap.txt"
    """
    awk 'NR==FNR { mask[FNR]=\$1; next } mask[FNR]==1' \
        ${dhs_mask} \
        <(grep -v '#' ${dhs_coordinates}) \
        | bedmap --indicator \
            <(zcat ${variants}) - > ${name}
    """
}

process mock_indicator {
    conda params.conda
    
    input:
        path dhs_coordinates
    
    output:
        tuple path(name), path(dhs_coordinates)
    
    script:
    name = "masterlist.all_ones.indicator.bed"
    """
    grep -v '#' ${dhs_coordinates} \
        | awk '{print 1}'
        > ${name}
    """
}


workflow {
    data = matricesListFromMeta()
       | splitMatrices // matrix_name, prefix, annotation_bool, dhs_coordinates


    data.first()
        | map(it -> it[3])
        | mock_indicator // mock_indicator, dhs_coordinates
        | map(it -> tuple("", "all_dhs", it[0], it[1])) // matrix_name, prefix, annotation_bool, dhs_coordinates
        | mix(data)
        | combine(
            Channel.fromPath(params.finemapped_variants_file)
        ) // matrix_name, prefix, annotation_bool, dhs_coordinates, variants
        | overlap_annotation

    data
        | collectFile(
            storeDir: "${params.outdir}/finemapping/",
            skip: 1,
            keepHeader: true
        ) { it -> [ "${it[0]}.finemapping_metadata.tsv", "prefix\tindicator\n${it[1]}\t${params.outdir}/finemapping/${it[0]}/${it[1]}.overlap.bed\n"] } // metadata
}