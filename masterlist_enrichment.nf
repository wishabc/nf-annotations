include { readMotifsList } from "./motif_enrichment"


process cut_matrix {
    conda params.conda
    tag "samples ${interval}"
    scratch true

    input:
        tuple val(sample_id), path(binary_matrix)

    output:
        tuple val(interval), path(name)

    script:
    end = sample_id + params.step - 1
    interval = "${sample_id}-${end}"
    name = "${interval}.cut_matrix.npy"
    """
    zcat ${binary_matrix} | cut -f${interval} > tmp.txt
    python3 $moduleDir/bin/convert_to_numpy.py tmp.txt ${name}
    """
}

process motif_hits_intersect {
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/output_file", pattern: "${motif_id}.hits.bed"

    input:
        tuple val(motif_id), path(moods_file), path(masterlist_file)

    output:
        tuple val(motif_id), path(indicator_file)

    script:
    indicator_file = "${motif_id}.hits.bed"
    """
    zcat ${moods_file} \
        | bedmap --indicator --sweep-all \
        --fraction-map 1 ${masterlist_file} - > ${indicator_file}
    """
}

process calc_index_motif_enrichment {
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/motif_stats_chunks"
    scratch true

    input:
        tuple val(motif_id), path(indicator_file), val(chunk_n), path(binary_matrix_chunk)
    
    output:
        path name

    script:
    name = "${motif_id}.${chunk_n}.enrichment.tsv"
    """
    python3 $moduleDir/bin/index_motif_enrichment.py  \
        ${binary_matrix_chunk} ${indicator_file} ${motif_id} \
        ${params.sample_names} ${chunk_n} > ${name}
    """
}


workflow indexEnrichment {
    chunks_count = file(params.sample_names).countLines().intdiv(params.step)

    c_mat = Channel.of(0..chunks_count) // 0, 1, 2, 3 ,4 , 5.. 9
        | map(it -> it * params.step + 1) // 1 201 401 601 1801
        | combine(file(params.binary_matrix)) // [chunk_n, binary_matrix]
        | cut_matrix // [1, 1-200.cut_matrix.npy], [201, 201-400.cut_matrix.npy]

    moods_scans = readMotifsList() // motif_id, motif_path
        | combine(file(params.masterlist_file)) // motif_id, motif_path, masterlist
        | motif_hits_intersect // motif_id, indicator_file
        | combine(c_mat) // motif_id, indicator_file, chunk_n, binary_matrix_chunk
        | calc_index_motif_enrichment // enrichment_file
        | collectFile(storeDir: "$launchDir/${params.outdir}", name: "enrichments.tsv")
}
// 
workflow {
    indexEnrichment()
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

// New process for generating plot for coeffs
process tf_by_components {
    conda params.pyconda
    publishDir "${params.outdir}/plot"

    input:
        path(all_coefs)
    
    output:
        path("${prefix}*.pdf")
    
    script:
    prefix = "components"
    """
    python $moduleDir/bin/plot_tf_by_components.py \
        ${all_coefs} \
        ${params.metadf} \
        ${prefix}
    """
}

// params.index_file = '/net/seq/data2/projects/sabramov/SuperIndex/dnase-0209/output/masterlist.filtered.bed'
// params.genome_fasta = '/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa'
// params.mappable_file = '/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.K76.mappable_only.bed'
process extract_gc_content {
	conda params.conda
	publishDir "${params.outdir}/gc_content"

	output:
		path gc_content_path

	script:
	gc_content_path = 'regions_gc_annotated.bed.gz'
	"""
	faidx -i nucleotide -b ${params.index_file} ${params.genome_fasta} \
		| awk -v OFS="\t" 'NR>1 { total =\$4+\$5+\$6+\$7+\$8; cg=\$6+\$7; print \$1, \$2-1, \$3,total, cg, cg/total;  }' \
		| bedmap --delim "\t" --echo --bases-uniq - ${params.mappable_file} \
		| paste - <(cut -f4,9 ${params.index_file}) \
	| bgzip -c > ${gc_content_path}
	"""
}

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

workflow gwasLogisticRegression {
    matrices = Channel.fromPath(params.matrix_file_pheno)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, file(row.phen_bed), params.masterlist_file2))
        | pheno_hits_intersect
    
}

workflow hyperGeom {
    params.pyconda = "/home/afathul/miniconda3/envs/motif_enrichment"

    params.reference_sample_meta = "/home/afathul/data2seq/motif_enrichment/odd_ratio/oddratio_dnase_3501/reference_samples/nmf/metadata_reference_sample.bed"
    params.all_samples_meta = "/home/afathul/data2seq/motif_enrichment/odd_ratio/oddratio_dnase_3501/all_samples/nmf/metadata_all_sample.bed"


    coeffs = Channel.fromPath("${params.moods_scans_dir}/*")
        | map (it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it, params.all_samples_meta))
        | motif_hits_intersect // motif_id, indicator
	    | geom_odd_ratio

    coeffs | map(it -> it[1])
        | collectFile(name: 'all.logodd.stats.tsv',
            storeDir: "${params.outdir}",
            skip: 1,
            sort: true,
            keepHeader: true)

}

process match_gc_background {
    conda params.pyconda
    tag "${motif_id}:${matrix_type}"
    publishDir "${params.outdir}/${matrix_type}"

    input:
        tuple val(motif_id), path(indicator_file), val(matrix_type), path(binary_matrix)
    
    output:
        tuple val(motif_id), val(matrix_type), path(name)
    
    script:
    name = "${motif_id}.${matrix_type}.z_score.tsv"
    """
    python $moduleDir/bin/subsample_proportion.py \
        ${motif_id} \
        ${matrix_type} \
        ${indicator_file} \
        ${name} \
        ${params.dhs_index_masterlist} \
        ${binary_matrix} \
        ${params.all_samples_meta} \
        ${params.metadata_file} \
        ${params.acc_proportion}

    """
}

workflow matchingBackground {
    params.pyconda = "/home/afathul/miniconda3/envs/motif_enrichment"

    matrices = Channel.of(tuple("NMF", file(params.binary_nmf)), tuple("DHS_Binary", file(params.binary_dhs_agid)))

    Channel.fromPath("${params.moods_scans_dir}/*")
        | map (it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it, params.all_samples_meta))
        | motif_hits_intersect // motif_id, indicator
        | combine(matrices)
	    | match_gc_background
        | map(it -> it[2])
        | collectFile(name: 'all.z_score.stats.tsv',
            storeDir: "${params.outdir}",
            skip: 1,
            sort: true,
            keepHeader: true
        )

}

workflow motifIndicator {
    params.pyconda = "/home/afathul/miniconda3/envs/motif_enrichment"
    params.all_samples_meta = "/home/afathul/data2seq/motif_enrichment/odd_ratio/oddratio_dnase_3501/all_samples/nmf/metadata_all_sample.bed"


    coeffs = Channel.fromPath("${params.moods_scans_dir}/*")
        | map (it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it, params.all_samples_meta))
        | motif_hits_intersect 

}