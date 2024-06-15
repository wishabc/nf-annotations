include { splitMatrices } from './ldsc'



process overlap_annotation {
    conda params.conda
    publishDir "${params.outdir}/${matrix_name}"
    tag "${prefix}"
    
    input:
        tuple val(matrix_name), val(prefix), path(dhs_mask), path(variants)
    
    output:
        tuple val(matrix_name), path(name)
    
    script:
    name = "${prefix}.overlap.bed"
    """
    awk 'NR==FNR { mask[FNR]=\$1; next } mask[FNR]==1' \
        ${dhs_mask} \
        ${params.masterlist_file} \
        | bedmap --indicator \
            <(zcat ${variants}) - > ${name}
    """
}

workflow {
    matrices = Channel.fromPath(params.matrices_list)
       | splitCsv(header:true, sep:'\t')
       | map(row -> tuple(row.matrix_name, file(row.matrix), file(row.sample_names)))
       | splitMatrices // matrix_name, annotation_name, annotation_bool
       | combine(
            Channel.fromPath("${params.finemapped_variants_file}")
       )
       | overlap_annotation
}