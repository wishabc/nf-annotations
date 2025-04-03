#!/usr/bin/env nextflow
include { matricesListFromMeta; splitMatrices } from './helpers'
params.conda = "$moduleDir/environment.yml"


workflow fromMatrix {
    take:
        matrices // prefix, matrix, names
    main:
        accessibility = Channel.fromPath("${params.template_run}/proportion_accessibility.tsv", checkIfExists: true)
        Channel.fromPath("${params.template_run}/motif_hits/*.hits.bed")
            | map(it -> tuple(it.name.replaceAll('.hits.bed', ''), it)) // motif_id, motif_hits
            | combine(matrices) // motif_id, motif_hits, prefix, matrix, names
            | combine(accessibility) // motif_id, motif_hits, prefix, matrix, names, accessibility
            | motifEnrichment
}

workflow fromMatricesList {
    Channel.fromPath(params.matrices_list)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.matrix_name, file(row.matrix), file(row.sample_names)))
        | fromMatrix
}


workflow fromBinaryMatrix {
    matrices = Channel.of(tuple("DHS_binary", file(params.index_anndata), file(params.peaks_mask))) // kinda hotfixes
        | extract_from_anndata // prefix, matrix, sample_names, dhs_names
    
    prop_accessibility = matrices
        | calc_prop_accessibility
    
    Channel.fromPath("${params.moods_scans_dir}/*") // result of nf-genotyping scan_motifs pipeline
        | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))
        | combine(matrices.map(it -> it[3]))
        | motif_hits_intersect // motif_id, indicator
        | combine(matrices)
        | combine(prop_accessibility)
        | motifEnrichment
}

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
    
    python3 $moduleDir/bin/motif_enrichment/assign_bins.py annotated.bed ${name}
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

workflow getRegionsSamplingPool {
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



process overlap_and_sample {
    conda params.conda
    tag "${annotation_name}"
    publishDir "${params.outdir}/motif_enrichment/per_category_samples/"
    label "med_mem"

    input:
        tuple val(annotation_name), path(annotation), path(annotation_coordinates), path(sampled_regions_pool), path(masterlist)

    output:
        tuple val(annotation_name), val('sampled'), path(name), emit: sampled
        tuple val(annotation_name), val('reference'), path(reference_dhs), emit: reference

    script:
    reference_dhs = "${motif_id}.${annotation_name}.reference_regions.bed"
    name = "${motif_id}.${annotation_name}.sampled_regions.bed"
    """
    python3 $moduleDir/bin/motif_enrichment/sample_regions.py \
        ${masterlist} \
        ${sampled_regions_pool} \
        ${motif_indicator} \
        ${annotation} \
        ${annotation_coordinates} \
        ${name} \
        ${reference_dhs}
    """
}

process motif_hits_intersect {
    tag "${motif_id}"
    conda params.conda
    // publishDir "${params.outdir}/motif_hits/${prefix}"

    input:
        tuple val(prefix), val(sampling_type), path(bed_file), val(motif_id), path(moods_file), 

    output:
        tuple val(prefix), val(sampling_type), path(name)

    script:
    name = "${prefix}.${motif_id}.${sampling_type}.stats.tsv"
    """
    zcat ${moods_file} \
        | bedmap --indicator --sweep-all \
        --fraction-map 1 <(grep -v '#' ${bed_file}) - > mask.txt
    
    echo "motif_id\tprefix\tsampling_type\toverlaps\ttotal" > ${name}
    awk -v OFS='\t' \
        '{ ones += \$1; total++ } \
        END { print "${motif_id}", "${prefix}", "${sampling_type}", ones, total }' \
        mask.txt >> ${name}

    """
}

workflow randomFromMatricesList {
    matricesListFromMeta()
        | splitMatrices
        | randomRegionsEnrichment
}


workflow randomRegionsEnrichment {
    take:
        annotations // prefix, annotation_mask, dhs_coordinates
    main:
        motif_hits = Channel.fromPath("${params.template_run}/motif_hits/index/*.hits.bed")
            | map(it -> tuple(it.name.replaceAll('.hits.bed', ''), it))
            | filter { it[0] == "M02739_2.00" }
        
        motifs_meta = Channel.fromPath("${params.moods_scans_dir}/*") // result of nf-genotyping scan_motifs pipeline
            | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))

        sampled = Channel.fromPath("${params.template_run}/motif_enrichment/sampled_regions_pool.parquet")
         
        annotated_masterlist = Channel.fromPath("${params.template_run}/motif_enrichment/index.annotated.bed")

        sampled_regions = annotations
            | combine(sampled)
            | combine(annotated_masterlist)
            | overlap_and_sample
        
        result = sampled_regions.sampled
            | mix(sampled_regions.reference)
            | combine(motifs_meta)
            | motif_hits_intersect
            | collectFile(
                storeDir: "${params.outdir}/motif_enrichment/",
                skip: 1,
                name: "enrichment_stats.tsv",
                keepHeader: true
            )
    emit:
        result
}
