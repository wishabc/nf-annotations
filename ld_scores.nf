#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"

process ld_scores {
	conda params.conda
	tag "${prefix}"

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
        | awk -v OFS='\t' '{ print \$1,\$2,\$3 }'  \
        | uniq >> variants.bed

	vcftools --geno-r2 \
		--gzvcf ${params.genotype_file} \
        --chr ${prefix} \
		--minDP 10 \
        --bed variants.bed \
		--ld-window 1 \
		--out ${prefix}
	"""
}

workflow ldScores {
    params.genotype_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-genotypes-round2/output/genotypes/all.filtered.snps.annotated.vcf.gz"


    pval_file = Channel.fromPath("/net/seq/data2/projects/sabramov/ENCODE4/dnase0620/dnase/cavs/output/by_sample/*.bed") 
        | ld_scores
}