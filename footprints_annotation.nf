
process annotate_motif_hits {
    conda params.conda
    tag "${prefix}"

    input:
        tuple val(prefix), path(motif)
    output:
        tuple val(prefix), path(name)
    
    script:
    name = "${prefix}.footprints_hits.bed"
    """
    bedtools intersect \
        -a ${params.footprints_index} \
        -b ${motifs} \
        -F 0.5 \
        -wa -wb > ${name}
    """
}

process sort_and_index {
    conda params.conda
    publishDir params.outdir
    label "highmem"

    input:
        path overlaps
    
    output:
        tuple path(name), path("${name}.tbi")
    
    script:
    name = "motif_overlaps.bed.gz"
    """
    cat ${overlaps} | sort-bed - \
        | bgzip -c > ${name}
    tabix ${name}
    """

}

workflow {
    params.footprints_index = "/net/seq/data2/projects/ENCODE4Plus/footprints/4078_Index/footprint_index_0521/output/unfiltered_masterlists/masterlist_DHSs_Altius_all_chunkIDs.bed"

    Channel.fromPath("${params.template_run}/motif_hits/*")
        | map(it -> tuple(it.name.replaceAll('.hits.bed', ''), it))
        | annotate_motif_hits
        | collectFile(
            name: "motifs_annotated.bed"
            sort: true,
        )
        | sort_and_index
}


// DEFUNC
process tmp_reformat {
    conda params.conda
    tag "${prefix}"
    publishDir params.outdir 

    input:
        tuple val(prefix), path(wrong_format_files)
    
    output:
        tuple val(prefix), path(name)

    script:
    name = "${prefix}.hg38.bed.gz"
    """
    zcat "${wrong_format_files}" \
        | sed 's/chrchr/chr/g' \
        | bgzip > ${name}
    """
}

workflow tmp {
    Channel.fromPath("/net/seq/data2/projects/sabramov/UKBB_2023/harmonized_n_samples/output/per_phenotype/**/*.hg38.bed.gz")
        | map(it -> tuple(it.name.replaceAll('.hg38.bed.gz', ''), it))
        | tmp_reformat
}