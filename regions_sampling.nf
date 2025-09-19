#!/usr/bin/env nextflow
params.conda = "$moduleDir/environment.yml"


process split_masterlist_in_chunks {
    conda params.conda

    input:
        tuple val(prefix), path(masterlist_file)

    output:
        path "${prefix}*.bed"

    script:
    """
    grep -v '#' ${masterlist_file} \
        | cut -f 1-3 \
        | split \
            -l 100000 \
            --suffix-length=4 \
            -d \
            --additional-suffix=.bed \
            - \
            ${prefix}
    """
}
process sample_matching_bg {

    conda '/home/afathul/miniconda3/envs/r-kernel'
    tag "${prefix}:${iter}"

    input:
        tuple val(iter), val(prefix), path(bed_file)
    
    output:
        tuple val(iter), path(name)

    script:
    name = "${prefix}.${iter}.bed"
    """
    Rscript $moduleDir/bin/motif_enrichment/delta_svm_match_bg.R \
        ${bed_file} \
        tmp.bed
    
    grep -w -F -f ${params.nuclear_chroms} tmp.bed > ${name}
    """
}

process annotate_regions {

    conda '/home/sabramov/miniconda3/envs/super-index'
    publishDir "${params.outdir}/motif_enrichment/", pattern: "index.annotated.bed"
    tag "${prefix}"
    label "high_mem"
    scratch true

    input:
        tuple val(prefix), path(bed_file)

    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.annotated.bed"
    """
    cat ${bed_file} \
        | cut -f 1-3 \
        | sort-bed - > tmp.bed

    faidx -i nucleotide \
        -b tmp.bed \
        ${params.genome_fasta} \
        | awk -v OFS="\t" \
            'NR>1 {  \
                total=\$4+\$5+\$6+\$7+\$8; \
                cg=\$6+\$7; \
                print \$1, \$2-1, \$3, \$3-\$2, cg, cg/total; }' > annotated.bed
    
    python3 $moduleDir/bin/motif_enrichment/assign_bins.py annotated.bed ${params.length_bins_bounds} ${name}
    """
}


process to_parquet {
    conda params.conda
    tag "${prefix}"
    publishDir "${params.outdir}/motif_enrichment/"
    label "high_mem"

    input:
        tuple val(prefix), path(bed_file)

    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.parquet"
    """
    python3 $moduleDir/bin/motif_enrichment/convert_to_parquet.py ${bed_file} ${name}
    """
}

workflow getRegionsPool {
    masterlist = Channel.fromPath(params.masterlist_file)
        | map(it -> tuple("index", it))

    chunks = masterlist
        | split_masterlist_in_chunks
        | flatten()
        | map(it -> tuple(it.baseName, it)) // prefix, bed_file
    
    sampled_bg = Channel.of(1..100)
        | combine(chunks)
        | sample_matching_bg
        | mix(chunks)
        | mix(masterlist)
        | annotate_regions
        | filter { it[0] != 'index' }
        | map(it -> it[1])
        | collectFile(
            sort: true,
            name: 'sampled_regions_pool.bed',
        )
        | map(it -> tuple('sampled_regions_pool', it))
        | to_parquet
}
