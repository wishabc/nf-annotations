
process filter_significant_hits {
    publishDir "${params.outdir}/per_phenotype/${phen_id}"
    conda params.conda
    tag "${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(bed_file)
    
    output:
        tuple val(phen_id), path(name)
    
    script:
    name = "${phen_id}.significant_hits.bed.gz"
    """
    zcat ${bed_file} \
        | awk -v OFS='\t' \
            '((NR == 1) || (\$10 >= 7.301)) {print }' \
        | bgzip -c > ${name}
    """
}


process sort_and_index {
    conda params.conda
    scratch true

    publishDir params.outdir

    input:
        path significant_hits_paths
    
    output:
        path name
    
    script:
    name = "significant_hits.bed.gz"
    """

    zcat \$(head -n 1 ${significant_hits_paths}) | head -1 || true > header.txt


    # Concatenate, sort, and merge with the header from the first file

    xargs -a ${significant_hits_paths} zcat \
        | tail -n +2 \
        | sort-bed - \
        | cat header.txt - \
        | bgzip -c > merged_and_sorted.gz
    """
}


workflow {
    phenotypes = Channel.fromPath(params.phenotypes_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, file(row.bed_file)))
        | filter_significant_hits
        | map(it -> it.toString())
        | collectFile(name: 'significant_hits_paths.txt', newLine: true)
        | sort_and_index
}



process annotate_ref_pop_with_gwas {


    input:
        tuple val(gwas_name), path(gwas_file)
    
    output:
        tuple val(gwas_name), path(name)
    
    script:
    name = "${gwas_name}.pop_annotated.bed"
    """
    head -1 ${params.ref_pop_file} > ${name}
    tail -n+2 ${params.ref_pop_file} \
        | bedmap -e 1 - <(zcat ${gwas_file} | grep -v '#') >> ${name}
    """
}

process get_n_per_bin {

    conda params.conda

    input:
        tuple val(gwas_name), path(pop_annotated_file)
    
    output:
        tuple val(gwas_name), path(name)

    script:
    name = "${gwas_name}.n_per_bin.txt"
    """
    python3 $moduleDir/bin/gwas_enrichment/cut_in_bins.py \
        ${pop_annotated_file} \
        ${name}
    """
}


process sample_from_ref_pop {

    publishDir "${params.outdir}/gwas_enrichment/${gwas_name}"
    conda params.conda

    input:
        tuple val(gwas_name), path(per_bin_counts), val(seed)
    
    output:
        tuple val(gwas_name), val(seed), path(name)

    script:
    name = "${gwas_name}.${seed}.sampled.bed"
    """
    python3 $moduleDir/bin/gwas_enrichment/sample.py \
        ${params.ref_pop_file} \
        ${per_bin_counts} \
        ${name}
    """
}

process extend_by_ld {

    publishDir "${params.outdir}/gwas_enrichment/${gwas_name}"
    conda params.conda

    input:
        tuple val(gwas_name), val(seed), path(sampled_variants)
    
    output:
        tuple val(gwas_name), val(seed), path(sampled_variants), path(ld_extended)
    
    script:
    ld_extended = "${gwas_name}.${seed}.ld_extended.bed"
    """
    bedops -e 1 ${sampled_variants} ${params.perfect_ld_variants} \
        | awk -v OFS="\t" '\$5==1 { print \$1, \$6, \$6+1, ".", \$5, \$2; }' \
        | sort-bed -> ${ld_extended}
    """
}


workflow sampleMatchedVariantsForTraits {

    seeds = Channel.from(1..params.n_samples) 
        //| map(seed -> seed * 1000)

    data = Channel.fromPath(params.gwas_files)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.gwas_name, file(row.gwas_file)))
        | annotate_ref_pop_with_gwas
        | get_n_per_bin
        | combine(seeds)
        | sample_from_ref_pop
        | extend_by_ld
        | collectFile(
            storeDir: "${params.outdir}/gwas_enrichment/",
            skip: 1,
            keepHeader: true,
        ) {
            it -> [
                "${it[0]}.sampled_file.txt", 
                "seed\tsampled\tld_extended\n${it[1]}\t${params.outdir}/gwas_enrichment/${it[0]}.${it[1]}.sampled.bed\t${params.outdir}/gwas_enrichment/${it[0]}.${it[1]}.ld_extended.bed\n"
            ]
        }
}


process overlap_with_sampled {

}


workflow overlapSampled {
    
}

































// DEFUNC

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