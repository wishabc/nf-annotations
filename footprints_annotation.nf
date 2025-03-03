
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
    zcat ${motif} \
        | bedtools intersect \
        -a ${params.footprints_index} \
        -b stdin \
        -F 0.9 -f 0.9 \
        -wa -wb \
        | cut -f1-4,12- > ${name}
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

    Channel.fromPath("${params.moods_scans_dir}/*") // result of nf-genotyping scan_motifs pipeline
        | map(it -> tuple(it.name.replaceAll('.moods.log.bed.gz', ''), it))
        | annotate_motif_hits
        | map(it -> it[1])
        | collectFile(
            name: "motifs_annotated.bed",
            sort: true
        )
        | sort_and_index
}
