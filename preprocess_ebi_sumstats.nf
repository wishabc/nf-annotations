
process download_file {
    tag "${phen_id}"
    publishDir "${params.outdir}/per_phenotype/${phen_id}"
    maxForks 4

    input:
        tuple val(phen_id), val(sumstats_link), val(metadata_link)
    
    output:
        tuple val(phen_id), path(name), path(metadata_name)
    
    script:
    name = "${phen_id}.tsv.gz"
    metadata_name = "${phen_id}.meta.yaml"
    """
    wget '${link}' -O ${name}
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

    input:
        tuple val(phen_id), path(sumstats_file), val(n_samples)

    output:
        tuple val(phen_id), path("${prefix}.sumstats.gz"), path("${prefix}.log")
    
    script:
    prefix = "${phen_id}.munge"
    """
    check_effect_allele_frequency=\$( \
        awk 'BEGIN {FS="\t"} \
            NR==1 {for (i=1; i<=NF; i++) \
                if (\$i == "effect_allele_frequency") col=i} 
                    NR > 1 && col && \$col != "NA" {found=1} 
            END {if (found) print "present"}'\
            "${sumstats_file}")

    effect_allele_frequency_flag=\$([ "\$check_effect_allele_frequency" == "present" ] && echo " --frq effect_allele_frequency" || echo "")

    
    python ${params.ldsc_scripts_path}/munge_sumstats.py \
        --sumstats ${sumstats_file} \
        --merge-alleles ${params.tested_snps} \
        --a1 effect_allele \
        --a2 other_allele \
        --snp variant_id \
        --N ${n_samples} \
        \${effect_allele_frequency_flag} \
        --out ${prefix}
    """
}

workflow {
    meta = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, row.ebi_link, row.ebi_meta_link))
        | download_file
}

workflow tmp {
    meta = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, file(row.sumstats_file), row.n_samples))
        | filter{ it[2] }
        | munge_sumstats
}