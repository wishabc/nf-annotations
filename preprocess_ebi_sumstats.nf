
process download_file {
    tag "${phen_id}"
    publishDir "${params.outdir}/per_phenotype/${phen_id}"
    maxForks 4

    input:
        tuple val(phen_id), val(link)
    
    output:
        tuple val(phen_id), path(name)
    
    script:
    name = "${phen_id}.tsv.gz"
    """
    wget '${link}' -O ${name}
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

workflow {
    meta = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, row.ebi_link, row.sumstats_exists))
        | filter{ it[2] == "TRUE" }
        | map(it -> tuple(it[0], it[1]))
        | download_file
}

workflow tmp {
    meta = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, row.ebi_meta_link, row.sumstats_exists))
        | filter{ it[2] == "TRUE" }
        | map(it -> tuple(it[0], it[1]))
        | download_meta
}