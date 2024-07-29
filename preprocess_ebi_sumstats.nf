
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
    publishDir "${params.outdir}/per_phenotype/${phen_id}"
    scratch true
    errorStrategy 'ignore'

    input:
        tuple val(phen_id), path(sumstats_file), val(n_samples)

    output:
        tuple val(phen_id), path("${prefix}.sumstats.gz"), path("${prefix}.log")
    
    script:
    prefix = "${phen_id}.munge"
    """
    check_column_presence() {
        local column_name=\$1
        local file_name=\$2

        awk -v col_name="\$column_name" '
            BEGIN {FS="\t"}
            NR==1 {for (i=1; i<=NF; i++) if (\$i == col_name) col=i}
            NR > 1 && col {print "present"; exit}
        ' "\$file_name"
    }

    effect_allele_frequency_flag=\$(check_column_presence "effect_allele_frequency" "$sumstats_file" && echo " --frq effect_allele_frequency" || echo "")

    # Check for rs_id column
    snp_column=\$(check_column_presence "rs_id" "$sumstats_file" && echo "rs_id" || echo "variant_id")


    
    python ${params.ldsc_scripts_path}/munge_sumstats.py \
        --sumstats ${sumstats_file} \
        --merge-alleles ${params.tested_snps} \
        --a1 effect_allele \
        --a2 other_allele \
        --snp \${snp_column} \
        --N ${n_samples} \
        \${effect_allele_frequency_flag} \
        --out ${prefix}
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
        | map(row -> tuple(row.phen_id, file(row.sumstats_file), row.n_samples, file("/net/seq/data2/projects/sabramov/EBI.Summary_statistics/output/per_phenotype/" + row.phen_id + ".munge.sumstats.gz")))
        | filter{ !it[3].exists() }
        | munge_sumstats
}