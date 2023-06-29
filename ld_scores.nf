#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"

process ld_scores {
	conda params.conda
	tag "${chromosome}"
	scratch true

	input:
		tuple path(snps_positions), val(chromosome)
	
	output:
		tuple val(chromosome), path(name)
	
	script:
	name = "${chromosome}.geno.ld"
	"""
    awk -v OFS='\t' '{ print \$1,\$2,\$3 }' ${snps_positions} \
        | uniq \
        | python3 reformat_to_r2.py > variants.bed

	vcftools --geno-r2 \
		--gzvcf ${params.genotype_file} \
        --chr ${chromosome} \
		--minDP 10 \
		--ld-window-bp 500000 \
		--chr ${chromosome} \
		--out ${chromosome}
    
    cat ${chromosome}.geno.ld \
        | awk -v OFS='\t' '{ print \$1,\$2-1,\$2,\$3,\$4,\$5 }' \
        | bedtools intersect -wa -a stdin -b variants.bed 
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