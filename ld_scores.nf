#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"

process ld_scores {
	conda params.conda
	tag "${chromosome}"
	scratch true

	input:
		tuple val(chromosome), path(snps_positions)
	
	output:
		tuple val(chromosome), path(name)
	
	script:
	name = "${chromosome}.geno.ld"
	"""
    echo "chrom chromStart  chromEnd" > variants.bed
    awk -v OFS='\t' '{ print \$1,\$2,\$3 }' ${snps_positions} >> variants.bed

	vcftools --geno-r2 \
		--gzvcf ${params.genotype_file} \
        --chr ${chromosome} \
		--minDP ${params.min_DP} \
		--bed variants.bed \
		--ld-window 1 \
		--chr ${chromosome} \
		--out ${chromosome}
	"""
}

workflow ldScores {
    params.genotype_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-genotypes-round2/output/genotypes/all.filtered.snps.annotated.vcf.gz"

    chroms = Channel.of(1..22)
		| map(it -> "chr${it}")

    pval_file = Channel.fromPath(params.pval_file) 
        | filterUniqVariants
        | combine(chroms)
        | ld_scores
        | map(it -> it[1])
        | collectFile(name: 'ld_scores.geno.ld', storeDir: "$launchDir/${params.outdir}")
}