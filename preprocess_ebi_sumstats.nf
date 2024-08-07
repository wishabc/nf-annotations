
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
    tag "${phen_id}:${ref_is_effect}"
    label "ldsc"
    scratch true
    errorStrategy 'ignore'
    publishDir "${params.outdir}/per_phenotype/${phen_id}"

    input:
        tuple val(phen_id), path(sumstats_file), val(n_samples)

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
        ${prefix}
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
        //| filter{ !it[3].exists() }
        | map(it -> tuple(it[0], it[1], it[2]))
        | filter{ it[2] > 0 }
        | munge_sumstats
        | map(it -> tuple(it[0], it[1]))
        | groupTuple(by: 0, size: 2, remainder: true)
        | choose_correct_orientation

}