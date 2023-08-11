#!/usr/bin/env nextflow
// include { filterTestedVariants } from "./main"

params.conda = "$moduleDir/environment.yml"

// TODO wrap in apptainer
process calc_ld {
    publishDir "${params.outdir}/ldsc/l2_logs", pattern: "${prefix}.log"
    publishDir "${params.outdir}/ldsc/l2", pattern: "${prefix}.l2.*"
    publishDir "${params.outdir}/ldsc/l2", pattern: "${annotation_file}"
    
    tag "${prefix}"
    scratch true
    conda params.ldsc_conda

    input:
        tuple val(prefix), path(annotation_file), val(is_baseline)
    
    output:
        tuple val(prefix), path("${prefix}.l2.*"), path("${prefix}.log")
    
    script:
    annot_type = is_baseline ? "" : "--thin-annot"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    
    awk 'NR>1 {print \$1}' \
        ${params.tested_snps} > tested_snps.txt

    ${params.ldsc_scripts_path}/ldsc.py \
        --print-snps tested_snps.txt \
        --ld-wind-cm 1.0 \
        --out ${prefix} \
        --bfile ${params.gtfiles}${chrom} \
        --annot ${annotation_file} \
        ${annot_type} \
        --l2
    """
}

// this version doesn't return h^2 estimates!
process run_ldsc_cell_types {
    conda params.ldsc_conda
    publishDir "${params.outdir}/ldsc/ldsc_coefs", pattern: "${name}"
    publishDir "${params.outdir}/ldsc/ldsc_logs", pattern: "${phen_id}.log"
    tag "${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file), path("data_files/*")
    
    output:
        tuple val(phen_id), path(name), path("${phen_id}.log")

    script:
    name = "${phen_id}.results.tsv"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    ls -1 data_files \
        | cut -d"." -f 1 \
        | sort \
        | uniq \
        | awk -v OFS='\\t' \
            '{ print \$1,"data_files/"\$1"." }' > per_sample.ldcts
    
    ${params.ldsc_scripts_path}/ldsc.py \
        --h2-cts ${sumstats_file} \
        --ref-ld-chr ${params.baseline_ld} \
        --ref-ld-chr-cts per_sample.ldcts \
        --frqfile-chr ${params.frqfiles} \
        --w-ld-chr ${params.weights} \
        --overlap-annot \
        --print-coefficients \
        --out ${phen_id}

    awk -v OFS='\t' \
        '{print (NR==1 ? "Phenotype_ID" : "${phen_id}"), \$0}' \
        ${phen_id}.cell_type_results.txt > ${name}
    """
}


process run_ldsc_single_sample {
    conda params.ldsc_conda
    publishDir "${params.outdir}/ldsc/ldsc_coefs${prefix}", pattern: "${name}.results"
    publishDir "${params.outdir}/ldsc/ldsc_logs/${prefix}", pattern: "${name}.log"
    tag "${prefix}:${phen_id}"
    //scratch true

    input:
        tuple val(phen_id), path(sumstats_file), val(prefix), path(ld_files)
    
    output:
        tuple val(prefix), val(phen_id), path("${name}.results"), path("${name}.log")

    script:
    name = "${prefix}.${phen_id}"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    ${params.ldsc_scripts_path}/ldsc.py \
        --h2 ${sumstats_file} \
        --ref-ld-chr ${params.baseline_ld},${prefix}. \
        --frqfile-chr ${params.frqfiles} \
        --w-ld-chr ${params.weights} \
        --overlap-annot \
        --print-coefficients \
        --out ${name}
    """
}


process collect_ldsc_results {
    scratch true
    publishDir "${params.outdir}"

    input:
        tuple path(ldsc_files), path(ldsc_logs)
    
    output:
        path name

    script:
    name = "ldsc_enrichments_results.tsv"
    """
    # copy header from the first file
    head -1 ${ldsc_files[0]} \
        | xargs -I % echo "group_name\tphenotype_id\t%" > result.txt
    
    # Aggregate the data
    echo ${ldsc_files} | tr ' ' '\n' > filelist.txt
    echo "h^2" > h2.stats
    while read line; do
        fname=\$(basename "\$line" .results)
        echo \$fname
        grep "Total Observed scale h2" \${fname}.log \
            | sed 's/[:(]/\t/g' \
            | sed 's/)//g' \
            | cut -f 3 >> h2.stats

        echo \$fname \
            | tr "." "\t" \
            | xargs -I % echo "%\t`tail -1 "\$line"`" >> result.txt
    done < filelist.txt
    paste result.txt h2.stats > ${name}
    """
}

// Make annotation workflows
process filter_cavs {

    input:
        path pval_file

    output:
        path "*.bed"
    
    script:
    """
    cat ${pval_file} \
        | grep -v '#' \
        | awk '\$NF <= ${params.fdr_tr} \
            {print> \$19".bed"}'
    """
}

process make_ldsc_annotation {
    conda params.ldsc_conda
    tag "${prefix}"
    scratch true

    input:
        tuple val(chrom), path(custom_annotation)

    output:
        tuple val(prefix), path(name)
    
    script:
    prefix = "${custom_annotation.simpleName}.${chrom}"
    baseannotation = "${params.base_ann_path}${chrom}.annot.gz"
    name = "${prefix}.annot.gz"
    """
    python ${params.ldsc_scripts_path}/make_annot.py \
        --bimfile ${params.gtfiles}${chrom}.bim \
        --bed-file ${custom_annotation} \
        --annot-file ${name}
    """
}

workflow LDSCcellTypes {
    take:
        ld_data
        sumstats_files
    main:
        out = sumstats_files
            | combine(
                ld_data.map(it -> it[1]).collect(sort: true)
            )
            | run_ldsc_cell_types
            | map(it -> it[1])
            | collectFile(
                storeDir: params.outdir,
                skip: 1,
                keepHeader: true,
                sort: true,
                name: 'ldsc_ct_results.tsv'
            )
    emit:
        out
}


workflow LDSC {
    take:
        ld_data
        sumstats_files
    main:
        out = sumstats_files
            | combine(ld_data)
            | run_ldsc_single_sample
            | map(it -> tuple(it[2], it[3]))
            | collect(sort: true, flat: false)
            | collect_ldsc_results
    emit:
        out
}

workflow fromAnnotations {
    take:
        annotations
    main:
        sumstats_files = Channel.fromPath(params.phenotypes_meta)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.phen_id, file(row.munge_sumstats_file)))
            | filter { it[1].exists() }

        ldsc_annotations = Channel.of(1..22)
            | combine(annotations)
            | make_ldsc_annotation

        ld_data = ldsc_annotations
            | combine(
                Channel.of(false)
            ) // prefix, annotation, is_baseline
            | calc_ld // prefix, ld, ld_log
            | join(ldsc_annotations) // prefix, ld, ld_log, annotation
            | map(it -> tuple(it[0], [it[1], it[3]].flatten()))
            | groupTuple(size: 22)
            | map(
                it -> tuple(it[0], it[1].flatten())
            )
        if (params.by_cell_type) {
            out = LDSCcellTypes(ld_data, sumstats_files)
        } else {
            out = LDSC(ld_data, sumstats_files) 
        }   
    emit:
        out
}

// Entry workflows
workflow calcBaseline {
    Channel.of(1..22)
        | map(it -> tuple(it, file("${params.base_ann_path}${it}.annot.gz", checkIfExists: true), true))
        | calc_ld
}

workflow fromPvalFiles {
    Channel.fromPath(params.result_pval_file) 
        | filter_cavs
        | flatten()
        | filter { it.countLines() >= 8000 }
        | fromAnnotations
   
}

// workflow {
//     custom_annotations = Channel.fromPath("${params.annotations_dir}/*.bed") 
//         | filterTestedVariants
//         | fromAnnotations
// }
