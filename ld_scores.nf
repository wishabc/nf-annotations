#!/usr/bin/env nextflow


process ld_scores {
	conda params.conda
	tag "${prefix}"

	input:
		tuple val(chrom), path(snps_positions)
	
	output:
		path name
	
	script:
    prefix = "${chrom}:${snps_positions.simpleName}"
	name = "${prefix}.geno.ld"
    additional_params = chrom == 'all' ? "" : "--chr ${chrom}"
 	"""
    echo "chrom chromStart  chromEnd" > variants.bed
    cat ${snps_positions} \
        | grep -v '^#' \
        | awk -v OFS='\t' '{ print \$1,\$2,\$3 }'  \
        | uniq >> variants.bed

	vcftools --geno-r2 \
		--gzvcf ${params.genotype_file} \
		--minDP 10 \
        --bed variants.bed \
		--out ${prefix} \
        ${additional_params} \
	"""
}


workflow byChromosome {
    Channel.of(1..22)
        | map(it -> "chr${it}")
        | combine(Channel.fromPath(params.pval_file))
        | ld_scores
        | collectFile(
            storeDir: params.outdir,
            keepHeader: true,
            sort: true,
            skip: 1,
            name: "ld_scores.geno.ld"
        )
}


workflow bySample {
    params.by_sample_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase0620/dnase.auto/output/by_sample/*.bed"
    Channel.fromPath(params.by_sample_file)
        | map(it -> tuple("all", it))
        | ld_scores
}