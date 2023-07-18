#!/usr/bin/env nextflow


process ld_scores {
	conda params.conda
	tag "${prefix}"
    publishDir "${params.outdir}"

	input:
		path snps_positions
	
	output:
		path name
	
	script:
    prefix = snps_positions.simpleName
	name = "${prefix}.geno.ld"
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
		--out ${prefix}
	"""
}

workflow {
    params.genotype_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-genotypes-round2/output/genotypes/all.filtered.snps.annotated.vcf.gz"

    params.pval_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase0620/dnase.auto/output/by_sample/*.bed"
    pval_file = Channel.fromPath(params.pval_file) 
        | ld_scores
}