#!/usr/bin/env nextflow
include { splitMatrices; matricesListFromMeta; convert_to_bed } from './helpers'
params.conda = "$moduleDir/environment.yml"


// TODO wrap in apptainer
process calc_ld {
    publishDir "${params.outdir}/ldsc/l2_logs", pattern: "${prefix}.log"
    publishDir "${params.outdir}/ldsc/l2", pattern: "${prefix}.l2.*"
    publishDir "${params.outdir}/ldsc/l2", pattern: "${new_annot}"
    
    tag "chr${chrom}:${group_id}"
    scratch true
    conda params.ldsc_conda
    label "ldsc"

    input:
        tuple val(matrix_prefix), val(group_id), val(chrom), path(annotation)
    
    output:
        tuple val(matrix_prefix), val(group_id), val(chrom), path("${prefix}.l2.*"), path("${prefix}.log"), path(new_annot)
    
    script:
    prefix = "${group_id}.${chrom}"
    annot_type = group_id == "baseline" ? "" : "--thin-annot"
    new_annot = "${prefix}.annot.gz"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    
    if [ "${annotation}" != "${new_annot}" ]; then
        cp ${annotation} ${new_annot}
    fi
    
    awk 'NR>1 {print \$1}' \
        ${params.tested_snps} > tested_snps.txt

    ${params.ldsc_scripts_path}/ldsc.py \
        --print-snps tested_snps.txt \
        --ld-wind-cm 1.0 \
        --out ${prefix} \
        --bfile ${params.gtfiles}${chrom} \
        --annot ${new_annot} \
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
    label "ldsc"

    input:
        tuple val(matrix_prefix), path("data_files/*"), val(phen_id), path(sumstats_file), val(baseline_ld)
    
    output:
        tuple val(matrix_prefix), val(phen_id), path(name), path("${phen_id}.log")

    script:
    name = "${phen_id}.results.tsv"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    python $moduleDir/bin/metadata_cl_ldsc.py \
        ${matrix_prefix} \
        data_files/ \
        per_sample.ldcts
    
    ${params.ldsc_scripts_path}/ldsc.py \
        --h2-cts ${sumstats_file} \
        --ref-ld-chr ${baseline_ld} \
        --ref-ld-chr-cts per_sample.ldcts \
        --frqfile-chr ${params.frqfiles} \
        --w-ld-chr ${params.weights} \
        --overlap-annot \
        --print-coefficients \
        --out ${phen_id}

    awk -v OFS='\t' \
        '{print (NR==1 ? "phen_id" : "${phen_id}"), \$0}' \
        ${phen_id}.cell_type_results.txt > ${name}
    """
}


process run_ldsc_single_sample {
    conda params.ldsc_conda
    publishDir "${params.outdir}/ldsc/ldsc_coefs_${prefix}", pattern: "${name}.results"
    publishDir "${params.outdir}/ldsc/ldsc_logs/${prefix}", pattern: "${name}.log"
    tag "${prefix}:${phen_id}"
    label "ldsc"
    scratch true

    // phen_id, sumstats_file, baseline_ld, matrix_prefix, group_id, ld_files
    input:
        tuple val(phen_id), path(sumstats_file), val(baseline_ld), val(matrix_prefix), val(prefix), path(ld_files)
    
    output:
        tuple val(matrix_prefix), val(prefix), val(phen_id), path("${name}.results"), path("${name}.log"), path("${name}.summary_result")

    script:
    name = "${prefix}.${phen_id}"
    additional_ld_chr = prefix != "baseline" ? ",${prefix}." : ""
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    ${params.ldsc_scripts_path}/ldsc.py \
        --h2 ${sumstats_file} \
        --ref-ld-chr ${baseline_ld}${additional_ld_chr} \
        --frqfile-chr ${params.frqfiles} \
        --w-ld-chr ${params.weights} \
        --overlap-annot \
        --print-coefficients \
        --out ${name}
    
    head -1 ${name}.results \
        | xargs -I % echo -e "group_id\tphen_id\t%\th^2\th^2_err\tldsc_path" > ${name}.summary_result

    grep "Total Observed scale h2" ${name}.log \
        | awk -F'[:()]' -v OFS='\t' \
        '{gsub(/[[:space:]]/, "", \$2); \
         gsub(/[[:space:]]/, "", \$3); \
            print \$2, \$3}' > h2_data.txt

    # Combine and format the output
    output_path="${params.outdir}/ldsc/ldsc_coefs_${prefix}/${name}.results"
    ldsc_annotation_output=\$(tail -n 1 ${name}.results)
    echo -e "${prefix}\t${phen_id}\t\$ldsc_annotation_output\t\$(cat h2_data.txt)\t\$output_path" >> ${name}.summary_result
    """
}

process make_ldsc_annotation {
    conda params.conda
    tag "chr${chrom}:${group_id}"
    publishDir "${params.outdir}/pipeline_annotations"
    scratch true

    input:
        tuple val(chrom), val(matrix_prefix), val(group_id), path(custom_annotation)

    output:
        tuple val(matrix_prefix), val(group_id), val(chrom), path(name)
    
    script:
    name = "${group_id}.${chrom}.annot.gz"
    oper = custom_annotation.name.endsWith('.bed') ? 'cat' : 'zcat'
    """
    echo ANNOT | gzip > ${name}

    cat ${params.gtfiles}${chrom}.bim \
        | awk -v OFS='\t' '{ print "chr"\$1,\$4-1,\$4,NR }' \
        | sort-bed - > sorted_bim.bed

    ${oper} ${custom_annotation} \
        | grep -v '#' \
        | cut -f1-3 \
        | sort-bed - \
        | bedmap --indicator --sweep-all sorted_bim.bed - \
        | paste sorted_bim.bed - \
        | sort -k4,4n \
        | cut -f 5 \
        | gzip >> ${name}
    """
}


workflow LDSCcellTypes {
    take:
        ld_data
        sumstats_files
    main:
        out = ld_data // matrix_prefix, group_id, ld_files
            | map(it -> tuple(it[0], it[2])) // matrix_prefix, ld_files
            | groupTuple()
            | map(it -> tuple(it[0], it[1].flatten())) // matrix_prefix, ld_files_flatten
            | combine(sumstats_files) // matrix_prefix, ld_files, phen_id, sumstats_file, baseline_ld
            | run_ldsc_cell_types // matrix_prefix, phen_id, result, log
            | collectFile(
                storeDir: params.outdir,
                skip: 1,
                keepHeader: true,
                sort: true,
            ) { it -> [ "${it[0]}.ldsc_cell_types_results.tsv", it[2].text ] }
    emit:
        out
}


workflow LDSC {
    take:
        ld_data
        sumstats_files
    main:
        out = sumstats_files // phen_id, sumstats_file, baseline_ld
            | combine(ld_data) // phen_id, sumstats_file, baseline_ld, matrix_prefix, group_id, ld_files
            | run_ldsc_single_sample // matrix_prefix, group_id, phen_id, result, log, annotation_result
            | collectFile(
                storeDir: params.outdir,
                skip: 1,
                keepHeader: true
            ) { it -> [ "${it[0]}.ldsc_enrichments_results.tsv", it[5].text ] }
    emit:
        out
}

workflow fromAnnotations {
    take:
        annotations
    main:
        sumstats_files = Channel.fromPath(params.phenotypes_meta)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.phen_id, file(row.munge_sumstats_file), params.baseline_ld)) 

        ld_data = Channel.of(1..22) // chroms, move to a separate file
            | combine(annotations) // matrix_prefix, annotation_name, annotation
            | make_ldsc_annotation // matrix_prefix, annotation_name, chrom, annotation
            | calc_ld //  matrix_prefix, annotation_name, chrom, ld, ld_log, annotation
            | map(it -> tuple(it[0], it[1], [it[3], it[5]].flatten())) // matrix_prefix, group_id, ld_files
            | groupTuple(size: 22, by: [0, 1]) // also count # of lines in the file
            | map(it -> tuple(it[0], it[1], it[2].flatten())) // matrix_prefix, group_id, ld_files_flatten

        if (params.by_cell_type) {
            out = LDSCcellTypes(ld_data, sumstats_files)
        } else {
            out = LDSC(ld_data, sumstats_files)
        }   
    emit:
        out
}

workflow fromMatrix {
    take:
        matrices // matrix_name, matrix, names, dhs_names
    main:
        out = matrices
            | splitMatrices // matrix_name, dhs_names, annotation_bool
            | convert_to_bed // matrix_name, group_id, annotation_bed
            | fromAnnotations
    emit:
        out
}


// Entry workflows
workflow {
    println "Looking for custom annotations in ${params.custom_annotations_meta}"
    custom_annotations = Channel.fromPath(params.custom_annotations_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.prefix, row.annotation_name, file(row.annotation)))
        | fromAnnotations
}

// Start from matrices list
workflow fromMatricesList {
    matricesListFromMeta()
        | fromMatrix
}


// Start from CAV calling pval file
process split_cell_specific_aggregation {
    conda params.conda
    publishDir "${params.outdir}/ldsc/annotations", pattern: "*${suffix}"

    input:
        path aggregated_pval_file

    output:
        path "*${suffix}"
    
    script:
    suffix = ".annotation.bed"
    """
    python3 $moduleDir/bin/split_cell_specific_aggregation.py \
        ${aggregated_pval_file} \
        --fdr_tr ${params.fdr_tr} \
        --suffix ${suffix}
    """
}

workflow fromAggregatedCavs {
    params.result_pval_file = "${params.outdir}/aggregated.${params.aggregation_key}.bed"
    pvals = Channel.fromPath(params.result_pval_file) 
        | split_cell_specific_aggregation
        | flatten()
        | filter { it.countLines() >= params.min_snps }
        | map(it -> tuple(params.aggregation_key, it.baseName, it)) // matrix_name, group_id, aggregated_pval_file
        | fromAnnotations
   
}



// workflow fromBinaryMatrix {
//     // Doesn't work just yet
//     Channel.fromPath(params.binary_matrix)
//         | map(it -> tuple("DHS_Binary", it, file(params.sample_names), file(params.masterlist_file)))
//         | fromMatrix
// }
