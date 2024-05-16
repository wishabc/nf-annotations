
process split_matrices {
    conda params.conda
    //publishDir "${params.outdir}"
    
    input:
        tuple val(matrix_name), path(matrix), path(sample_names)
    
    output:
        path "${matrix_name}.*.txt"
    
    script:
    """
    awk -F'\t' -v names_file=${sample_names} -v matrix_name=${matrix_name} \
        'BEGIN {
            while ((getline line < names_file) > 0) {
                names[NR] = line
            }
            close(names_file)
        }
        {
            # For each column, print the value to the corresponding file
            for (i = 1; i <= NF; i++) {
                file = matrix_name"."names[i]".txt"
                print \$i > file
            }
        }' "${matrix}"
    """
}

process overlap_annotation {
    conda params.conda
    publishDir "${params.outdir}"
    
    input:
        tuple val(prefix), path(dhs_mask), path(variants)
    
    output:
        tuple val(component_id), path(name)
    
    script:
    name = "comp${component_id}.overlap.bed"
    """
    awk 'NR==FNR { mask[FNR]=\$1; next } mask[FNR]==1' \
        ${dhs_mask} \
        ${params.masterlist_file} \
        | bedmap --echo --indicator \
            ${variants} - > ${name}
    """
}

workflow {
    matrices = Channel.fromPath(params.matrices_list)
       | splitCsv(header:true, sep:'\t')
       | map(row -> tuple(row.matrix_name, file(row.matrix), file(row.sample_names)))
       | split_matrices
       | flatten()
       | map(it -> tuple(it.baseName, it))
       | combine(
            Channel.fromPath("${params.finemapped_variants_file}")
       )
       | overlap_annotation
       

}