#!/usr/bin/env nextflow
include { filterUniqVariants; readMotifsList } from "./motif_enrichment"
// Put in the Apptainer
params.conda = "$moduleDir/environment.yml"

params.window = 20

params.sarus = "/home/sabramov/projects/ENCODE4/sarus/sarus/sarus-latest.jar"
params.ape = "/home/sabramov/projects/ENCODE4/sarus/macro-perfectos-ape/ape.jar"

process cut_sequence {
    conda params.conda

    input:
        path pval_file
    
    output:
        path fasta_file
    
    script:
    fasta_file = "snps_context.fa"
    """
    cat ${pval_file} | awk -v OFS='\t' \
     '{print \$1,\$2-${params.window},\$3+${params.window}," "\$1"@"\$3"@"\$5"@"\$6"@ref"}' > regions.bed
    
    bedtools getfasta -fi ${params.genome_fasta_file} \
        -bed regions.bed -nameOnly | sed 's/[YMSRDVBHWK]/N/g' > ${fasta_file}
    
    python3 $moduleDir/bin/change_context_to_alt.py ${fasta_file} ${params.window} out.fa
    
    cat out.fa >> ${fasta_file}
    """
}

process precalc_thresholds {
    conda params.conda

    input:
        path "./pwms/*"
    
    output:
        path "./motif_thresholds/*"
    
    script:
    """
    java -cp ${params.ape} ru.autosome.ape.PrecalculateThresholds ./pwms ./motif_thresholds
    """
}

process scan_with_sarus {
    conda params.conda
    tag "${motif_id}"

    input:
        tuple val(motif_id), path(pwm_path), path(pval_mapping), path(fasta_file)
    
    output:
        tuple val(motif_id), path(name)
    
    script:
    name = "${motif_id}.sarus.tsv"
    """
    java -cp ${params.sarus} ru.autosome.SARUS \
        ${fasta_file} \
        ${pwm_path} \
        -10000000 \
        --transpose \
        --threshold-mode score \
        --pvalues-file ${pval_mapping} \
        --output-scoring-mode logpvalue > sarus.log
    python3 $moduleDir/bin/parse_sarus_log.py sarus.log \
            `cat ${pwm_path} | wc -l` ${params.window} ${name} ${pwm_path}
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


workflow runSarus {
    take:
        pval_files
        motifs
    main:
        unique_variants = filterUniqVariants(pval_files) | cut_sequence
        thresholds = motifs
            | map(it -> it[1])
            | collect(sort: true)
            | precalc_thresholds
            | flatten()
            | map(it -> tuple(it.baseName, it))
            | view()

        out = motifs 
            | join(thresholds)
            | combine(unique_variants)
            | take(3)
            | scan_with_sarus
    emit:
        out
}

workflow {
    pvals_files = Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
    motifs_list = readMotifsList()
    runSarus(pvals_files, motifs_list)
}
