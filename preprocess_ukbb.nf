#!/usr/bin/env nextflow

process check_file {
    tag "${phen_id}"

    input:
        tuple val(phen_id), path(sumstats)
    
    output:
        tuple val(phen_id), stdout
    
    script:
    """
    if [ -f '${sumstats}' ]; then
        md5sum '${sumstats}' | awk '{ print \$1 }' | tr -d "\n"
    else
        echo 'no file'
    fi
    """
}

process download_file {
    tag "${phen_id}"
    publishDir "${params.outdir}/${phen_id}"
    maxForks 4

    input:
        tuple val(phen_id), val(aws_link), val(fname), val(md5)
    
    output:
        tuple val(phen_id), path(fname)
    
    script:
    """
    wget '${aws_link}'
    if [ "\$(md5sum '${fname}' | awk '{ print \$1 }' | tr -d '\n')" != "${md5}" ]; then
        echo "md5 are not matching!!!"
        exit 1
    fi
    """
}


process convert_to_hg38 {
    tag "${phen_id}"
    publishDir "${params.outdir}/${phen_id}", pattern: "${hg38_bed}"
    conda params.conda
    scratch true
    label "med_mem"

    input:
        tuple val(phen_id), path(sumstats), val(n_samples)
    
    output:
        tuple val(phen_id), path(hg38_bed)
    
    script:
    hg38_bed = "${phen_id}.hg38.bed.gz"
    """
    zcat ${params.variants_manifest} \
        | cut -f-5,11 > variants.txt
    
    # returns file with columns: 
    # chr, start, end, ref, alt, Beta, Beta_se, P, neglog10_p
    zcat ${sumstats} \
        | paste - variants.txt \
        | python3 $moduleDir/bin/reformat_sumstats.py \
            hg19.bed \
            ${params.population}

    echo -e "#chr\tstart\tend\tSNP\tref\talt\tBeta\tBeta_se\tP\tneglog10_p\tINFO\tphen_id\tN" > tmp.bed

    # don't do all the operations if file is empty
    if [ -s hg19.bed ]; then
        liftOver -bedPlus=3 \
            hg19.bed \
            ${params.chain_file} \
            unsorted \
            unMapped
    
        sort-bed unsorted \
            | awk -v OFS='\t' \
                '{print \$0, "${phen_id}", "${n_samples}"}' >> tmp.bed
    fi

    bgzip -c tmp.bed> ${hg38_bed}
    """
}

process filter_significant_hits {
    publishDir "${params.outdir}/${phen_id}"
    conda params.conda
    tag "${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(bed_file)
    
    output:
        tuple val(phen_id), path(name)
    
    script:
    name = "${phen_id}.significant_hits.bed"
    """
    zcat ${bed_file} \
        | awk -v OFS='\t' \
            '((NR == 1) || (\$10 >= 7.301)) {print }' \
        > ${name}
    """
}

process munge_sumstats {
    conda params.ldsc_conda
    tag "${phen_id}"
    publishDir "${params.outdir}/${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file)

    output:
        tuple val(phen_id), path("${prefix}.sumstats.gz"), path("${prefix}.log")
    
    script:
    prefix = "${phen_id}.munge"
    """
    python ${params.ldsc_scripts_path}/munge_sumstats.py \
        --sumstats ${sumstats_file} \
        --merge-alleles ${params.tested_snps} \
        --a1 alt \
        --a2 ref \
        --out ${prefix}
    """
}

process sort_and_index {
    conda params.conda
    scratch true

    publishDir params.outdir

    input:
        path siginficant_gwas_hits
    
    output:
        path name
    
    script:
    name = "significant_hits.bed.gz"
    """
    head -1 ${siginficant_gwas_hits} > result.bed
    tail -n +2 ${siginficant_gwas_hits} \
        | sort-bed - >> result.bed
    
    bgzip -c result.bed > ${name}
    """
}

workflow {
    // Code needs some adjustment for meta study!
    params.population = "EUR"
    params.variants_manifest = "/net/seq/data2/projects/GWAS/UKBB_2023/sabramov/full_variant_qc_metrics.txt.bgz"
    params.chain_file = "/home/ehaugen/refseq/liftOver/hg19ToHg38.over.chain.gz"
    data = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id,
            file(row.sumstats_file),
            row["n_cases_${params.population}"],
            row.pops))
        | filter { it[3] =~ /${params.population}/ }
        | map(it -> tuple(*it[0..2]))
        | convert_to_hg38

    data
        | filter_significant_hits
        | map(it -> it[1])
        | collectFile(
            sort: true,
            keepHeader: true,
            skip: 1,
            name: 'significant_hits.bed',
        )
        | sort_and_index
    

    munge_sumstats(data)
}

workflow checkData {
    meta = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
                row.phen_id,
                file(row.sumstats_file),
                row.aws_link,
                row.filename,
                row.md5_hex)
            )
    
    meta.map(it -> tuple(it[0], it[1]))
        | check_file // phen_id, md5
        | join(
            meta.map(
                it -> tuple(it[0], *it[2..4])
            )
        ) // phen_id, md5_file, aws_link, fname, md5_meta
        | filter { it[1] != it[4] }
        | map(it -> tuple(it[0], *it[2..4]))
        | download_file
}