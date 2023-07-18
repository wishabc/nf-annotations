#!/usr/bin/env nextflow


process ld_scores {
	conda params.conda
	tag "${prefix}"
    publishDir "${params.outdir}"

	input:
		tuple val(chrom), path(snps_positions)
	
	output:
		path name
	
	script:
    prefix = snps_positions.simpleName
	name = "${chrom}:${prefix}.geno.ld"
    additional_params = chrom == 'all' ? "" : "--chr ${chrom}"
 	"""
    echo "chrom chromStart  chromEnd" > variants.bed
    cat ${snps_positions} \
        | grep -v '^#' \
        | awk -v OFS='\t' '\$NF == "True" { print \$1,\$2,\$3 }'  \
        | uniq >> variants.bed

	vcftools --geno-r2 \
		--gzvcf ${params.genotype_file} \
		--minDP 10 \
        --bed variants.bed \
		--ld-window 1 \
		--out ${prefix} \
        ${additional_params} \
	"""
}


workflow byChromosome {
    Channel.of(1..22)
        | map(it -> "chr${it}")
        | combine(Channel.fromPath(params.pval_file))
        | ld_scores
}


workflow bySample {
    params.by_sample_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase0620/dnase.auto/output/by_sample/*.bed"
    Channel.fromPath(params.by_sample_file)
        | map(it -> tuple("all", it))
        | ld_scores
}