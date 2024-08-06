
process download_file {
    tag "${phen_id}"
    publishDir "${params.outdir}/per_phenotype/${phen_id}"
    maxForks 20

    input:
        tuple val(phen_id), val(sumstats_link), val(metadata_link)
    
    output:
        tuple val(phen_id), path(name), path(metadata_name)
    
    script:
    name = "${phen_id}.tsv.gz"
    metadata_name = "${phen_id}.meta.yaml"
    """
    wget '${sumstats_link}' -O ${name}
    wget '${metadata_link}' -O ${metadata_name}
    """
}

process download_meta {
    tag "${phen_id}"
    publishDir "${params.outdir}/per_phenotype/${phen_id}"
    maxForks 4

    input:
        tuple val(phen_id), val(link)
    
    output:
        tuple val(phen_id), path(name)
    
    script:
    name = "${phen_id}.meta.yaml"
    """
    wget '${link}' -O ${name}
    """
}

process munge_sumstats {
    conda params.ldsc_conda
    tag "${phen_id}"
    label "ldsc"
    scratch true
    errorStrategy 'ignore'
    publishDir "${params.outdir}/per_phenotype/${phen_id}", pattern: "${prefix}.log"

    input:
        tuple val(phen_id), path(sumstats_file), val(n_samples), val(ref_is_effect)

    output:
        tuple val(phen_id), path("${prefix}.sumstats.gz"), path("${prefix}.log")
    
    script:
    prefix = "${phen_id}.${ref_is_effect}.munge"
    """
    python ${moduleDir}/bin/preprocess_munge.py \
        ${sumstats_file} \
        ${params.ldsc_scripts_path}/munge_sumstats.py \
        ${params.tested_snps} \
        ${n_samples} \
        ${prefix} \

    """
}

process choose_correct_orientation {
    conda params.conda
    tag "${phen_id}"
    publishDir "${params.outdir}/per_phenotype/${phen_id}"
    errorStrategy 'ignore'

    input:
        tuple val(phen_id), path(sumstats_files)

    output:
        tuple val(phen_id), path(name)
    
    script:
    name = "${phen_id}.munge.sumstats.gz"
    """
    max_count=0
    max_file=""

    for file in ${sumstats_files}; do
        count=\$(zcat "${file}" \
            | awk -F'\t' \
                '\$3 != "" && \$3 != "NA" {count++} END {print count}')

        if (( count > max_count )); then
            max_count=\$count
            max_file=\$file
        fi
    done
    cp \$max_file ${name}
    """
}

workflow {
    meta = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, row.ebi_link, row.ebi_meta_link))
        | filter { it[0] != "null" }
        | download_file
}

workflow tmp {
    meta = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, file(row.sumstats_file), row.n_samples.toInteger(), file("/net/seq/data2/projects/sabramov/EBI.Summary_statistics/output/per_phenotype/${row.phen_id}/${row.phen_id}.munge.sumstats.gz")))
        | filter{ !it[3].exists() & it[2] > 0 }
        | map(it -> tuple(it[0], it[1], it[2]))
        | munge_sumstats
        | map(it -> tuple(it[0], it[1]))
        | groupTuple(by: 0, size: 2, remainder: true)
        | choose_correct_orientation

}