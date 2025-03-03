process extract_from_anndata {
    conda params.conda
    label "high_mem"

    input:
        tuple val(prefix), path(index_anndata), path(peaks_mask)
    
    output:
        tuple val(prefix), path(name), path(sample_names), path(masterlist_file)
    
    script:
    name = "binary_matrix.npz"
    sample_names = "sample_names.txt"
    masterlist_file = "masterlist.bed"
    """
    python $moduleDir/bin/extract_from_anndata.py \
        ${index_anndata} \
        ${peaks_mask} \
        ${name} \
        ${sample_names} \
        ${masterlist_file}
    """
}

process split_matrices {
    conda params.conda
    //publishDir "${params.outdir}"
    tag "${matrix_name}"
    
    input:
        tuple val(matrix_name), path(matrix), path(sample_names), path(dhs_coordinates)
    
    output:
        tuple val(matrix_name), path(dhs_coordinates), path("${matrix_name}.*.txt")
    
    script:
    """
    python3 $moduleDir/bin/split_matrix.py \
        ${matrix} \
        ${sample_names} \
        ${matrix_name}
    """
}


workflow splitMatrices {
    take:
        matrices_list
    main:
        out = matrices_list
            | split_matrices
            | transpose()
            | map(it -> tuple(it[0], it[2].baseName, it[2], it[1]))
    emit:
        out
}


workflow matricesListFromMeta {
    main:
        meta = Channel.fromPath(params.matrices_list)
            | splitCsv(header:true, sep:'\t')
            | map(
                row -> tuple(
                    row.matrix_name, 
                    file(row.anndata_path),
                    file(row.peaks_mask),
                    file(row.matrix), 
                    file(row.sample_names),
                )
            )
        out = meta 
            | map(it -> tuple(it[0], it[1], it[2])) // matrix_name, anndata, peaks_mask
            | extract_from_anndata // matrix_name, matrix, names, dhs_names
            | join(meta) // matrix_name, matrix, names, dhs_names, anndata, peaks_mask, matrix, sample_names
            | map(it -> tuple(it[0], it[6], it[7], it[3])) // matrix_name, matrix, names, dhs_names
    emit:
        out
}

