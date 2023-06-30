#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"

process ld_scores {
	conda params.conda
	tag "${snps_positions.simpleName}"

	input:
		tuple path(snps_positions)
	
	output:
		tuple path(name)
	
	script:
	name = "${snps_positions.simpleName}.geno.ld"
	"""
    echo "chrom chromStart  chromEnd" > variants.bed
    awk -v OFS='\t' '{ print \$1,\$2,\$3 }' ${snps_positions} \
        | uniq >> variants.bed

	vcftools --geno-r2 \
		--gzvcf ${params.genotype_file} \
        --chr ${chromosome} \
		--minDP 10 \
        --bed variants.bed \
		--ld-window 1 \
		--chr ${chromosome} \
		--out ${chromosome}
	"""
}

workflow ldScores {
    params.genotype_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-genotypes-round2/output/genotypes/all.filtered.snps.annotated.vcf.gz"


    pval_file = Channel.fromPath("/net/seq/data2/projects/sabramov/ENCODE4/dnase0620/dnase/cavs/output/by_sample/*.bed") 
        | filterUniqVariants
        | ld_scores
        | map(it -> it[1])
        | collectFile(
            sort: true,
            skip: 1,
            keepHeader: true,
            name: 'ld_scores.geno.ld',
            storeDir: "$launchDir/${params.outdir}")
}