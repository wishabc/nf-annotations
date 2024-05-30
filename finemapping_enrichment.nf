
process split_matrices {
    conda params.conda
    //publishDir "${params.outdir}"
    tag "${matrix_name}"
    
    input:
        tuple val(matrix_name), path(matrix), path(sample_names)
    
    output:
        path "${matrix_name}.*.txt"
    
    script:
    """
    python3 $moduleDir/bin/split_matrix.py \
        ${matrix} \
        ${sample_names} \
        ${matrix_name}
    """
}

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
       | split_matrices
       | flatten()
       | map(it -> tuple(it.simpleName, it.baseName, it))
       | combine(
            Channel.fromPath("${params.finemapped_variants_file}")
       )
       | overlap_annotation
}