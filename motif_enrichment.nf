#!/usr/bin/env nextflow

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

process filter_uniq_variants {
    conda params.conda
    scratch true

    input:
        path(pval_files)
    output:
        path name
    script:
    name = "merged.snps.sorted.uniq.bed"
    """
    echo "${pval_files}" | tr " " "\n" > filelist.txt
    while read file; do
        # extract chr, start, end, ID, ref, alt
        cat \$file | awk -v OFS='\t' '\$1 ~ /^[^;#]/ {print \$1,\$2,\$3,\$4,\$5,\$6}' >> merged_files.bed
    done < filelist.txt
    sort-bed merged_files.bed | uniq > ${name}
    """
}

process motif_counts {
    scratch true
    tag "${motif_id}"
    conda params.conda

    input:
        tuple val(motif_id), path(pwm_path), path(moods_file)
        path pval_file

    output:
        path counts_file

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

process get_motif_stats {
    tag "${pval_file.simpleName}"
    publishDir "${params.outdir}/motif_stats", pattern: "${motif_stats}"
    conda params.conda
    errorStrategy 'terminate'

    input:
        tuple path(pval_file), path(counts_file)

    output:
        path motif_stats
    
    script:
    motif_stats = "${pval_file.simpleName}.motif_stats.tsv"
    """
    # Counts file
    python3 ${projectDir}/bin/motif_stats.py  \
        ${pval_file} ${counts_file} > ${motif_stats}
    """
}
process collect_files {
    conda params.conda
    
    input:
        path motif_counts_file
    output:
        path name

    script:
    name = "merged_motif_chunks.bed"
    """
    for file in ${motif_counts_file}; do
        cat \$file >> merged_motifs.bed
    done
    sort-bed merged_motifs.bed > ${name}
    """
}

workflow calcEnrichment {
    take:
        moods_scans
        pvals_files
    main:
        pval_file = filter_uniq_variants(pvals_files.collect(sort: true))
        counts = motif_counts(moods_scans, pval_file) 
            | collate(30)
            | collect_files
        motif_ann = pvals_files
            | combine(counts)
            | get_motif_stats
            | collectFile(storeDir: "${params.outdir}") { it.text }
    emit:
        motif_ann
}

params.redo_moods = false
workflow readMoods {
    // Check if moods_scans_dir exists, if not run motifEnrichment pipeline
    main:
        motifs = Channel.fromPath(params.motifs_list)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.motif, file(row.motif_file)))
        if (file(params.moods_scans_dir).exists() && !params.redo_moods) {
            moods_scans = motifs.map(it -> tuple(it[0], it[1], "${params.moods_scans_dir}/${it[0]}.moods.log.bed.gz"))
        } else {
            moods_scans = scan_with_moods(motifs)
        }
    emit:
        moods_scans
}

workflow {
    pvals = Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
        | collect(sort: true)
        | flatten()
    moods_scans = readMoods()
    calcEnrichment(moods_scans, pvals)
}


process motif_hits_intersect {
    publishDir "${params.outdir}/counts", pattern: "${counts_file}"
    tag "${motif_id}"
    conda params.conda

    input:
        tuple val(motif_id), path(moods_file), path(index_file)

    output:
        tuple val(motif_id), path(counts_file)

    script:
    counts_file = "${motif_id}.hits.bed"
    """
    zcat ${moods_file} | bedmap --indicator --sweep-all --fraction-map 1 ${index_file} - > ${counts_file}
    """
}

process cut_matrix {
    conda params.conda
    tag "samples ${interval}"
    scratch true

    input:
        val sample_id

    output:
        tuple val(interval), path(name)

    script:
    end = sample_id + params.step - 1
    interval = "${sample_id}-${end}"
    name = "${interval}.cut_matrix.npy"
    """
    zcat ${params.binary_matrix} | cut -f${interval} > tmp.txt
    python3 $moduleDir/bin/convert_to_numpy.py tmp.txt ${name}
    """
}


process calc_index_motif_enrichment {
    tag "${motif_id}"
    conda params.conda
    publishDir "${params.outdir}/motif_stats_chunks"
    scratch true

    input:
        tuple val(motif_id), path(counts_file), val(sample_id),  path(matrix)
    
    output:
        tuple val(sample_id), path(name)

    script:
    name = "${motif_id}.${sample_id}.enrichment.tsv"
    """
    python3 $moduleDir/bin/index_motif_enrichment.py  \
        ${matrix} ${counts_file} ${motif_id} ${params.sample_names} ${sample_id} > ${name}
    """

}

workflow indexEnrichment {
    samples_count = file(params.sample_names).countLines().intdiv(params.step)
    sample_names = Channel.of(0..samples_count)
        | map(it -> it * params.step + 1)
        | toInteger()

    index = Channel.fromPath(file(params.index_file))

    moods_scans = readMoods().map(it -> tuple(it[0], it[2]))

    c_mat = cut_matrix(sample_names)
    out = motif_hits_intersect(moods_scans.combine(index)).combine(c_mat)
        | calc_index_motif_enrichment
        | collectFile(name: 'motif_enrichment.tsv', 
                      storeDir: "$launchDir/${params.outdir}")
}

workflow debug_counts {
    samples_count = file(params.sample_names).countLines().intdiv(params.step)
    sample_names = Channel.of(0..samples_count)
        | map(it -> it * params.step + 1)
        | toInteger()
    index = Channel.fromPath(file(params.index_file))

    moods_scans = readMoods().map(it -> tuple(it[0], it[2]))

    c_mat = cut_matrix(sample_names)
    out = motif_hits_intersect(moods_scans.combine(index)
            .filter { !file("$launchDir/${params.outdir}/counts/${it[0]}.hits.bed").exists() })
}

workflow debug2 {
    samples_count = file(params.sample_names).countLines().intdiv(params.step)
    sample_names = Channel.of(0..samples_count).map(it -> it * params.step + 1)
    index = Channel.fromPath(file(params.index_file))

    c_mat = cut_matrix(sample_names)
    motif_hits = readMoods()
        | map(it -> tuple(it[0], file("$launchDir/${params.outdir}/counts/${it[0]}.hits.bed")))
        | combine(c_mat)
        | filter { !file("$launchDir/${params.outdir}/motif_stats_chunks/${it[0]}.${it[2]}.enrichment.tsv").exists() }
        | calc_index_motif_enrichment
}