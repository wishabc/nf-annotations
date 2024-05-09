process gwas_logistic_regression {
    conda params.r_conda
    tag "${prefix}"
    publishDir "${params.outdir}/metrics", pattern: "${prefix}.metrics.tsv"
    publishDir "${params.outdir}/coeffs", pattern: "${prefix}.coeff.tsv"

    input:
        tuple val(phen_id), path(indicator_file), path(matrix)
    
    output:
        tuple val(motif_id), path("${prefix}.metrics.tsv"), path("${prefix}.coeff.tsv")
    
    script:
    prefix = "${motif_id}.${n_components}"
    """
    Rscript $moduleDir/bin/pheno_enrichment.R \
        ${matrix} \
        ${indicator_file} \
        ${motif_id} \
        ${n_components} 

    """
}

process geom_odd_ratio {
    conda params.pyconda
    tag "${motif_id}"
    publishDir "${params.outdir}/logodd"

    input:
        tuple val(motif_id), path(indicator_file)
    
    output:
        tuple val(motif_id), path(name)
    
    script:
    prefix = "${motif_id}"
    name = "${motif_id}.stats.tsv"
    """
    python $moduleDir/bin/log_oddratio.py \
        ${motif_id} \
        ${indicator_file} \
        ${name} \
        ${params.nmf_matrix} \
        ${params.metadata_file}

    """
}

workflow logisticRegression {
    params.r_conda = "/home/afathul/miniconda3/envs/r-kernel"
    params.pyconda = "/home/afathul/miniconda3/envs/motif_enrichment"
    params.masterlist_file = "/net/seq/data2/projects/afathul/motif_enhancement/masterlist.filtered.bed"
    params.matrix = "/net/seq/data2/projects/afathul/motif_enhancement/bin_new_unweight_full.16.H.npy"
    params.matrix_file = "/net/seq/data2/projects/afathul/motif_enhancement/modified_meta.tsv"
    params.gc_content_file = "/net/seq/data2/projects/afathul/motif_enhancement/regions_gc_annotated.bed.gz"

    matrices = Channel.fromPath(params.matrix_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.ncomponents, file(row.matrix))) // n_comp, matrix

    coeffs = Channel.fromPath("${params.moods_scans_dir}/*")
        | map (it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it, params.masterlist_file2))
        | motif_hits_intersect // motif_id, indicator
        | combine(matrices) // motif_id, indicator, n_comp, matrix
	    | logistic_regression

    coeffs | map(it -> it[1])
        | collectFile(name: 'all.metrics.tsv',
            storeDir: "${params.outdir}",
            skip: 1,
            sort: true,
            keepHeader: true)

    coeffs 
        | map(it -> it[2])
        | collectFile(name: 'all.coeff.tsv',
            storeDir: "${params.outdir}",
            skip: 1,
            sort: true,
            keepHeader: true)
        // tf_by_components
    
    // extract_gc_content
}

// Workflow 
process logistic_regression {
    conda params.r_conda
    tag "${prefix}"
    publishDir "${params.outdir}/metrics", pattern: "${prefix}.metrics.tsv"
    publishDir "${params.outdir}/coeffs", pattern: "${prefix}.coeff.tsv"

    input:
        tuple val(motif_id), path(indicator_file), val(n_components), path(matrix)
    
    output:
        tuple val(motif_id), path("${prefix}.metrics.tsv"), path("${prefix}.coeff.tsv")
    
    script:
    prefix = "${motif_id}.${n_components}"
    """
    Rscript $moduleDir/bin/motif_enrichment.R \
        ${matrix} \
        ${indicator_file} \
        ${motif_id} \
        ${n_components} 

    """
}

workflow gwasLogisticRegression {
    matrices = Channel.fromPath(params.matrix_file_pheno)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, file(row.phen_bed), params.masterlist_file2))
        | pheno_hits_intersect
    
}