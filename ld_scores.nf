#!/usr/bin/env nextflow

params.by_sample_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase0620/dnase.auto/output/by_sample/*.bed"

// sort-bed - \
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
        --ld-window-bp 50000 \
        --bed variants.bed \
		--out ${prefix} \
        ${additional_params} \
	"""
}
process sort {
    conda params.conda
    input:
        path ld_scores
    
    output:
        path name
    
    script:
    name = "${ld_scores.baseName}.sorted.bed"
    """
    tail -n+2 ${ld_scores} \
        | awk -v OFS='\t' '{print \$1,\$2-1,\$2,\$3,\$4,\5}' \
        | sort-bed - > ${name}
    """
}
process intersect_with_variants {
    conda params.conda
	tag "${variants_file.simpleName}"

    input:
        tuple path(ld_scores), path(variants_file)
    
    output:
        path name
    
    script:
    name = "${variants_file.simpleName}.bed"
    """
    # Format: chr, start, end, end2, distance, start_es, end_es, sample, chr, start, end, end3 ...
    python3 $moduleDir/bin/find_neighbors.py ${variants_file} \
        | sort-bed - \
        | bedtools intersect -a stdin -b ${ld_scores} -wa -wb -sorted \
        | awk -v OFS='\t' '\$4==\$12 { print; }' > ${name}
    """
}


workflow byChromosome {
    samples = Channel.fromPath(params.by_sample_file)
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
        | sort
        | combine(samples)
        | collectFile(
            storeDir: params.outdir,
            name: "ld_scores.annotated_samples.geno.ld"
        )
}


workflow bySample {
    samples = Channel.fromPath(params.by_sample_file)
        | map(it -> tuple("all", it))
        | ld_scores
}