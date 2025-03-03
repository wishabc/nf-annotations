include { splitMatrices; matricesListFromMeta } from './helpers'


process overlap_annotation {
    conda params.conda
    publishDir "${params.outdir}/${matrix_name}"
    tag "${prefix}"
    
    input:
        tuple val(matrix_name), val(prefix), path(dhs_mask), path(dhs_coordinates), path(variants)
    
    output:
        tuple val(matrix_name), path(name)
    
    script:
    name = "${prefix}.overlap.bed"
    """
    awk 'NR==FNR { mask[FNR]=\$1; next } mask[FNR]==1' \
        ${dhs_mask} \
        <(grep -v '#' ${dhs_coordinates}) \
        | bedmap --indicator \
            <(zcat ${variants}) - > ${name}
    """
}

workflow {
    data = matricesListFromMeta()
       | splitMatrices // matrix_name, prefix, annotation_bool, dhs_coordinates
       | combine(
            Channel.fromPath(params.finemapped_variants_file)
       ) // matrix_name, prefix, annotation_bool, dhs_coordinates, variants
    data // 
       | overlap_annotation
    
    data
        | collectFile(
            storeDir: params.outdir
        ) { it -> [ "${it[0]}.metadata.bed", "${it[1]}\t${params.outdir}/${it[0]}/${it[1]}.overlap.bed"] } // metadata
}