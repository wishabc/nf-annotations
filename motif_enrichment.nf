#!/usr/bin/env nextflow

include { filterTestedVariants } from "./main"
// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"


process scan_with_moods {
    conda params.conda
    tag "${motif_id}"
    scratch true
    publishDir "${params.moods_scans_dir}", pattern: "${name}", mode: "move"

    input:
        tuple val(motif_id), path(pwm_path)

    output:
        tuple val(motif_id), path(pwm_path), path(name)
    
    script:
    name = "${motif_id}.moods.log.bed.gz"
    moods_params = file(params.bg_file).exists() ? "--lo-bg `cat ${params.bg_file}`" : ""
    """
    moods-dna.py --sep ";" -s ${params.alt_fasta_file} \
        --p-value ${params.motif_pval_tr} \
        ${moods_params} \
        -m "${pwm_path}" \
        -o moods.log

    
    cat moods.log | awk '{print \$1}' > chroms.txt

    cat moods.log \
        | cut -d";" -f2- \
        | sed 's/;\$//g' \
        | awk -v FS=";" -v OFS="\t" \
            '{ print \$2, \$2+length(\$5), \$1, \$4, \$3, \$5; }' \
        | sed 's/".pfm"/""/g' \
        | paste chroms.txt - \
        | sort-bed - \
        | bgzip -c \
        > ${name}
    """
}

process motif_counts {
    // scratch true
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/counts", pattern: "${counts_file}"

    input:
        tuple val(motif_id), path(pwm_path), path(moods_file), path(pval_file)

    output:
        tuple val(motif_id), path(counts_file)

    script:
    counts_file = "${motif_id}.counts.bed"
    """
    echo -e "#chr\tstart\tend\tID\tref\talt\tmotif\toffset\twithin\tstrand\tref_score\talt_score\tseq" > ${counts_file}
    zcat ${moods_file} \
        | bedmap \
            --skip-unmapped \
            --sweep-all \
            --range ${params.flank_width} \
            --delim "|" \
            --multidelim ";" \
            --echo \
            --echo-map <(cat ${pval_file}) \
            - \
        | python3 $projectDir/bin/parse_variants_motifs.py \
            ${params.genome_fasta_file} \
            ${pwm_path} \
        >> ${counts_file}
    """
}

process tabix_index {
    conda params.conda
    publishDir "${params.outdir}"
    label "high_mem"


    input:
        path counts

    output:
        tuple path(name), path("${name}.tbi")
    script:
    name = "all_counts.merged.bed.gz"
    """
    head -1 ${counts} > ${name}
    sort-bed ${counts} >> ${name}
    bgzip ${name}
    tabix ${name}
    """
}

process calc_enrichment {
    tag "${motif_id}"
    conda params.conda

    input:
        tuple val(motif_id), path(counts_file), path(pval_file)

    output:
        tuple val(motif_id), path(pval_file), path(name)
    
    script:
    name = "${motif_id}.stats.tsv"
    """
    # Counts file
    python3 ${projectDir}/bin/motif_stats.py  \
        ${pval_file} ${counts_file} ${name} \
        --flank_width ${params.flank_width} \
        --fdr_tr ${params.fdr_tr}
    """
}

// process split_by_sample {

//     input:
//         path motif_file

//     output:
//         path "*.bed"
    
//     script:
//     """
//     cat ${pval_file} | grep -v '#' | awk '{print> \$19".bed"}'
//     """
// }

workflow motifCounts {
    take:
        data
    main:
        out = readMotifsList()
            | map(it -> tuple(it[0], it[1], "${params.moods_scans_dir}/${it[0]}.moods.log.bed.gz"))
            | combine(data)
            | motif_counts
        out 
            | map(it -> it[1])
            | collectFile(
                name: "all.counts.bed",
                keepHeader: true,
                skip: 1
            ) 
            | tabix_index
    emit:
        out
}

workflow readMotifsList {
    main:
        scans = Channel.fromPath(params.motifs_list)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.motif, file(row.motif_file)))
    emit:
        scans
}

// ------------ Entry workflows -------------------
workflow scanWithMoods {
    readMotifsList() | scan_with_moods
}

workflow {
    Channel.fromPath("${params.by_sample_pval_files}/*.bed")
        | collect(sort: true)
        | filterTestedVariants
        | motifCounts // motif_hits, motif_hits_index
        | combine(Channel.fromPath(params.result_pval_file))
        | calc_enrichment
        | map(it -> it[2])
        | collectFile(
            storeDir: params.outdir,
            name: "motif_enrichmnent.tsv",
            keepHeader: true,
            skip: 1
        )
}