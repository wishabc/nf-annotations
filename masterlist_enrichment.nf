#!/usr/bin/env nextflow
include { matricesListFromMeta; splitMatrices } from './helpers'
params.conda = "$moduleDir/environment.yml"


process motif_hits_intersect {
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/motif_hits/${prefix}"

    input:
        tuple val(motif_id), path(moods_file), val(prefix), path(bed_file)

    output:
        tuple val(motif_id), path(indicator_file)

    script:
    indicator_file = "${motif_id}.hits.bed"
    """
    zcat ${moods_file} \
        | bedmap --indicator --sweep-all \
        --fraction-map 1 <(grep -v '#' ${bed_file}) - > ${indicator_file}
    """
}

process motif_enrichment_z_score {
    conda params.conda
    tag "${motif_id}:${matrix_type}"
    publishDir "${params.outdir}/motif_enrichment/${matrix_type}"
    label "med_mem"

    input:
        tuple val(motif_id), path(indicator_file), val(matrix_type), path(binary_matrix), path(sample_names), path(accessibility_proportion)
    
    output:
        tuple val(matrix_type), val(motif_id), path(name)
    
    script:
    name = "${motif_id}.${matrix_type}.z_score.tsv"
    """
    python $moduleDir/bin/subsample_proportion.py \
        ${motif_id} \
        ${indicator_file} \
        ${name} \
        ${binary_matrix} \
        ${accessibility_proportion} \
        --sample_names ${sample_names} \
        --n_bins ${params.matching_bins}
    """
}

workflow motifEnrichment {
    take:
        data
    main: 
        out = data
            | motif_enrichment_z_score
            | collectFile(
                storeDir: params.outdir,
                skip: 1,
                sort: true,
                keepHeader: true
            ) { it -> [ "${it[0]}.z_score_stats.tsv", it[2].text ] }
    emit:
        out
}

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
                print \$1, \$2, \$3, \$3-\$2, cg, cg/total; }' > annotated.bed
    
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

    motifs_meta = Channel.fromPath("${params.moods_scans_dir}/*") // result of nf-genotyping scan_motifs pipeline
        | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))

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
        | branch {
            masterlist: it[0] == 'index'
            sampled: true
        }
    
    sampled_bg.sampled
        | map(it -> it[1])
        | collectFile(
            sort: true,
            name: 'sampled_regions_pool.bed',
        )
        | map(it -> tuple('sampled_regions_pool', it))
        | to_parquet

    sampled_bg.masterlist
        | combine(motifs_meta)
        | map(it -> tuple(it[2], it[3], it[0], it[1]))
        | motif_hits_intersect
}



process overlap_and_sample {
    conda params.conda
    tag "${motif_id}:${annotation_name}"
    publishDir "${params.outdir}/motif_enrichment/per_motif_samples/${motif_id}"
    label "med_mem"

    input:
        tuple val(motif_id), path(motif_indicator), val(annotation_name), path(annotation), path(annotation_coordinates), path(sampled_regions_pool), path(masterlist)

    output:
        tuple val(motif_id), val(annotation_name), path(name)

    script:
    name = "${motif_id}.${annotation_name}.sampled_regions.bed"
    """
    python3 $moduleDir/bin/motif_enrichment/sample_regions.py \
        ${sampled_regions_pool} \
        ${masterlist} \
        ${motif_indicator} \
        ${annotation} \
        ${annotation_coordinates} \
        ${name}
    """
}

workflow randomFromMatricesList {
    matricesListFromMeta()
        | splitMatrices
        | map(it -> tuple(it[1], it[2], it[3]))
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

        sampled = Channel.fromPath("${params.template_run}/motif_enrichment/sampled_regions_pool.parquet.bed")
         
        annotated_masterlist = Channel.fromPath("${params.template_run}/motif_enrichment/index.annotated.bed")

        sampled_regions = motif_hits
            | combine(annotations)
            | combine(sampled)
            | combine(annotated_masterlist)
            | overlap_and_sample
            | combine(motifs_meta, by: 0)
            | map(it -> tuple(it[0], it[2], it[1].baseName, it[1]))
            | motif_hits_intersect
            // | count_number_of_hits
            // | collectFile(
            //     storeDir: "${params.outdir}/motif_enrichment/${it[0]}/",
            //     skip: 1,
            //     sort: true,
            //     keepHeader: true
            // )
    emit:
        motif_hits_intersect.out
}

// DEFUNC
// process calc_prop_accessibility {
//     conda params.conda
//     publishDir "${params.outdir}"
//     label "high_mem"
    
//     input:
//         tuple val(id), path(binary_matrix), path(sample_names), path(masterlist_file)

//     output:
//         path name
    
//     script:
//     name = "proportion_accessibility.tsv"
//     """
//     python $moduleDir/bin/calc_prop_accessibility.py \
//         ${binary_matrix} \
//         ${masterlist_file} \
//         ${name} \
//         --samples_weights ${params.sample_weights}
//     """
// }