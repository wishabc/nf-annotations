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
    scratch true
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/counts"

    input:
        tuple val(motif_id), path(pwm_path), path(moods_file), path(pval_file)

    output:
        tuple val(motif_id), path(counts_file)

    script:
    counts_file = "${motif_id}.counts.bed"
    """
    zcat ${moods_file} | bedmap \
        --skip-unmapped \
        --sweep-all \
        --range 20 \
        --delim "|" \
        --multidelim ";" \
        --echo \
        --echo-map <(cat ${pval_file}) \
        -    \
        | python $projectDir/bin/parse_variants_motifs.py \
            ${params.genome_fasta_file} \
            ${pwm_path} \
        > ${counts_file}
    """
}

process tabix_index {
    conda params.conda
    publishDir "${params.outdir}"

    input:
        path counts

    output:
        tuple path(name), path("${name}.tbi")
    script:
    name = "all_counts.merged.bed.gz"
    """
    sort-bed ${counts} | bgzip -c > ${name}
    tabix ${name}
    """
}

process get_motif_stats {
    tag "${pval_file.simpleName}:${prefix}"
    conda params.conda

    input:
        tuple path(pval_file), val(prefix), path(counts_file)

    output:
        tuple path(pval_file), path(motif_stats)
    
    script:
    motif_stats = "${pval_file.baseName}.${prefix}.stats.tsv"
    """
    # Counts file
    python3 ${projectDir}/bin/motif_stats.py  \
        ${pval_file} ${counts_file} > ${motif_stats}
    """
}

// Development workflows

workflow calcEnrichment {
    take:
        counts
        pvals_files
    main:
        motif_ann = pvals_files
            | combine(counts)
            | get_motif_stats
            | collectFile(
                    storeDir: "${params.outdir}/stats"
                ) { it -> ["${it[0].simpleName}.motif_stats.txt", it[1].text] }
    emit:
        motif_ann
}

workflow motifCounts {
    take:
        pval_file
    main:
        uniq_vars = filterTestedVariants(pvals_files)

        counts =  readMotifsList()
            | map(it -> tuple(it[0], it[1], "${params.moods_scans_dir}/${it[0]}.moods.log.bed.gz"))
            | combine(uniq_vars)
            | motif_counts
            | map(it -> it[1])
            | collectFile(name: "all.counts.bed") 
            | tabix_index
    emit:
        counts
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
    Channel.fromPath(params.pval_file) | motifCounts
}

workflow cavsEnrichment {
    pvals_files = Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))

    calcEnrichment("${params.outdir}/all_counts.merged.bed.gz", pvals_files)
}
