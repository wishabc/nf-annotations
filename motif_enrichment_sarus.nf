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
        '{print \$1,\$2-${params.window},\$3+${params.window}," "\$1"@"\$3"@"\$4"@"\$5"@"\$6"@ref"}' > regions.bed
    
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
        path "motif_thresholds/*.thr"
    
    script:
    """
    java -cp ${params.ape} ru.autosome.ape.PrecalculateThresholds \
        ./pwms ./motif_thresholds --transpose
    """
}

process scan_with_sarus {
    conda params.conda
    tag "${motif_id}"
    publishDir "${params.outdir}/sarus_logs"

    input:
        tuple val(motif_id), path(pwm_path), path(pval_mapping), path(fasta_file)
    
    output:
        tuple val(motif_id), path(name), path(pwm_path)
    
    script:
    name = "${motif_id}.sarus.log"
    """
    java -cp ${params.sarus} ru.autosome.SARUS \
        ${fasta_file} \
        ${pwm_path} \
        -10000000 \
        --transpose \
        --threshold-mode score \
        --pvalues-file ${pval_mapping} \
        --output-scoring-mode logpvalue > ${name}
    """
}

process parse_log {
    conda params.conda
    tag "${motif_id}"
    publishDir "${params.outdir}/sarus"
    label "med_mem"

    input:
        tuple val(motif_id), path(sarus_log), path(pwm_path)
    
    output:
        tuple val(motif_id), path(name)

    script:
    name = "${motif_id}.sarus.tsv"
    """
    python3 $moduleDir/bin/parse_sarus_log.py ${sarus_log} \
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
    name = "all_counts.signif.bed.gz"
    """
    cat ${counts} | awk '((NR==1) || (\$8==1)) {print;}' \
        | sort-bed - \
        | bgzip -c > ${name}
    tabix ${name}
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

        out = motifs 
            | join(thresholds)
            | combine(unique_variants)
            | scan_with_sarus
            | parse_log
        
        out.map(it -> it[1])
            | collectFile(name: 'all.sarus.tsv') 
            | tabix_index
    emit:
        out
}

workflow {
    pvals_files = Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
    motifs_list = readMotifsList()
    runSarus(pvals_files, motifs_list)
}
