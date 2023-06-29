#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"

process ld_scores {
	conda params.conda
	tag "${chromosome}"

	input:
		tuple path(snps_positions), val(chromosome)
	
	output:
		tuple val(chromosome), path(name)
	
	script:
	name = "${chromosome}.ld.bed"
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

    chroms = Channel.of(1..22)
		| map(it -> "chr${it}")

    pval_file = Channel.fromPath(params.pval_file) 
        | filterUniqVariants
        | combine(chroms)
        | ld_scores
        | map(it -> it[1])
        | collectFile(
            sort: true
            name: 'ld_scores.geno.ld',
            storeDir: "$launchDir/${params.outdir}")
}