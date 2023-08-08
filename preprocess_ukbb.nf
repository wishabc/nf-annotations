#!/usr/bin/env nextflow

// Intended to work only on Altius server
process check_file {
    tag "${phen_id}"

    input:
        tuple val(phen_id), path(sumstats)
    
    output:
        tuple val(phen_id), stdout
    
    script:
    """
    if [ -f "${sumstats}" ]; then
        md5sum ${sumstats} | awk '{ print \$1 }'
    else
        echo 'no file'
    fi
    """
}

process download_file {
    tag "${phen_id}"
    publishDir "${params.outdir}/downloads/${phen_id}"

    input:
        tuple val(phen_id), val(md5_file), val(aws_link), val(fname)
    
    output:
        tuple val(phen_id), path(fname)
    
    script:
    """
    wget ${aws_link}
    md5=\$(md5sum ${fname} | awk '{ print \$1 }')
    if [ "\$(md5sum ${fname} | awk '{ print \$1 }')" != ];
        echo "md5 are not matching!!!"
        exit 1
    fi
    """
}

process convert_to_hg38 {
    tag "${phen_id}"
    publishDir "${params.outdir}/${phen_id}"
    module "kentutil:bedops"

    input:
        tuple val(phen_id), path(sumstats)
    
    output:
        tuple val(phen_id), path(hg38_bed), emit: bed
        tuple val(phen_id), path(hg38_sumstats), emit: sumstats
    
    script:
    hg38_bed = "${phen_id}.hg38.bed"
    hg38_sumstats = "${phen_id}.hg38.sumstats.gz"
    """
    zcat ${sumstats} \
        | awk  -F'\t' -v OFS='\t' \
        'NR > 1 {print "chr"\$1,\$2-1,\$2,\$9,NR-1}' \
        > hg19.bed

    liftOver -bedPlus=3 \
        hg19.bed \
        ${params.chain_file} \
        .unsorted \
        .unMapped
    
    echo -e "#chr\tstart\tend\tpval" > ${hg38_bed}
    sort-bed .unsorted >> ${hg38_bed}

    zcat ${sumstats} | head -1 
    """
}

// process collect_significant_hits {


// }

process munge_sumstats {
    conda params.ldsc_conda
    tag "${phen_id}"
    publishDir "${params.outdir}/${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file)

    output:
        tuple val(phen_id), path("${prefix}.sumstats.gz")
    
    script:
    prefix = "${phen_id}.munge"
    """
    python ${params.ldsc_scripts_path}/munge_sumstats.py \
        --sumstats ${sumstats_file} \
        --merge-alleles ${params.tested_snps} \
        --out ${prefix}
    """
}

params.ukbb_meta = "/net/seq/data2/projects/GWAS/UKBB_2023/sabramov/UKBB.metadata+sumstats.080823.tsv"

workflow checkData {
    meta = Channel.fromPath(params.ukbb_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(
                row.phenotype_id,
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
        | map(it -> tuple(*it[0..4]))
        | download_file
}

workflow {
    params.chain_file = "/home/ehaugen/refseq/liftOver/hg19ToHg38.over.chain.gz"
    Channel.fromPath(params.ukbb_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phenotype_id, file(row.sumstats_file)))
        
    convert_to_hg38
}