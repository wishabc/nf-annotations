#!/usr/bin/env nextflow
include { split_matrices } from './finemapping_enrichment'

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
        tuple val(group_id), val(chrom), path(annotation)
    
    output:
        tuple val(group_id), val(chrom), path("${prefix}.l2.*"), path("${prefix}.log"), path(new_annot)
    
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
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file), val(baseline_ld), path(filepaths)
    
    output:
        tuple val(phen_id), path(name), path("${phen_id}.log")

    script:
    name = "${phen_id}.results.tsv"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    cat ${filepaths} \
        | cut -d"." -f 1 \
        | sort \
        | uniq \
        | awk -v OFS='\\t' \
            '{ print \$1,"data_files/"\$1"." }' > per_sample.ldcts
    
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

    input:
        tuple val(phen_id), path(sumstats_file), val(baseline_ld), val(prefix), path(ld_files)
    
    output:
        tuple val(prefix), val(phen_id), path("${name}.results"), path("${name}.log"), path("${name}.summary_result")

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

// Make annotation workflows
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

process make_ldsc_annotation {
    conda params.conda
    tag "chr${chrom}:${group_id}"
    publishDir "${params.outdir}/pipeline_annotations"
    scratch true

    input:
        tuple val(chrom), val(group_id), path(custom_annotation)

    output:
        tuple val(group_id), val(chrom), path(name)
    
    script:
    name = "${group_id}.${chrom}.annot.gz"
    """
    echo ANNOT | gzip > ${name}

    cat ${params.gtfiles}${chrom}.bim \
        | awk -v OFS='\t' '{ print "chr"\$1,\$4-1,\$4,NR }' \
        | sort-bed - > sorted_bim.bed

    cut -f1-3 ${custom_annotation} \
        | sort-bed - \
        | bedmap --indicator --sweep-all sorted_bim.bed - \
        | paste sorted_bim.bed - \
        | sort -k4,4n \
        | cut -f 5 \
        | gzip >> ${name}
    """
}

process convert_to_bed {

    conda params.conda
    tag "${mask_name}:${prefix}"

    input:
        tuple val(mask_name), path(mask)

    output:
        tuple val(prefix), path(name)
    
    script:
    prefix = "${mask_name.replaceAll(/\./, '__')}"
    name = "${prefix}.bed"
    """
    awk -v OFS='\t' \
        'NR==FNR {mask[NR]=\$1; next} \
            mask[FNR] == 1' \
            ${mask} ${params.masterlist_file} \
        | cut -f 1-3 > ${name}

    """

}

workflow LDSCcellTypes {
    take:
        ld_data
        sumstats_files
        out_prefix
    main:
        dat = ld_data
            | map(it -> it[1].toString())
            | collectFile(name: 'all.paths.txt', newLine: true)

        out = sumstats_files
            | combine(dat)
            | run_ldsc_cell_types
            | map(it -> it[1])
            | collectFile(
                storeDir: params.outdir,
                skip: 1,
                keepHeader: true,
                sort: true,
                name: "${out_prefix}.ldsc_cell_types_results.tsv"
            )
    emit:
        out
}


workflow LDSC {
    take:
        ld_data
        sumstats_files
        out_prefix
    main:
        out = sumstats_files
            | combine(ld_data)
            | run_ldsc_single_sample
            | map(it -> it[4])
            | collectFile(
                name: "${out_prefix}.ldsc_enrichments_results.tsv",
                storeDir: params.outdir,
                skip: 1,
                keepHeader: true
            )
    emit:
        out
}

workflow fromAnnotations {
    take:
        annotations
        out_prefix
    main:
        sumstats_files = Channel.fromPath(params.phenotypes_meta)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.phen_id, file(row.munge_sumstats_file), params.baseline_ld))
            | filter { it[1].exists() }

        ld_data = Channel.of(1..22)
            | combine(annotations)
            | make_ldsc_annotation // group_id, chrom, annotation
            | calc_ld //  group_id, chrom, ld, ld_log, annotation
            | map(it -> tuple(it[0], [it[2], it[4]].flatten()))
            | groupTuple(size: 22)
            | map(it -> tuple(it[0], it[1].flatten()))

        if (params.by_cell_type) {
            out = LDSCcellTypes(ld_data, sumstats_files, out_prefix)
        } else {
            out = LDSC(ld_data, sumstats_files, out_prefix) 
        }   
    emit:
        out
}

workflow fromMatrix {
    take:
        matrices
        out_prefix
    main:
        data = matrices
            | split_matrices
            | flatten()
            | map(it -> tuple(it.baseName, it)) // mask_name, mask
            | convert_to_bed // group_id, annotation
        out = fromAnnotations(data, out_prefix)
    emit:
        out
}

// Entry workflows
workflow {
    custom_annotations = Channel.fromPath("${params.annotations_dir}/*.bed") 
        | map(it -> tuple(it.baseName, it)) // group_id, custom_annotation
    
    params.custom_annotation_name = params.custom_annotation_name ?: "custom_annotations"

    fromAnnotations(custom_annotations, params.custom_annotation_name)
}


workflow fromMatricesList {
    matrices = Channel.fromPath(params.matrices_list)
       | splitCsv(header:true, sep:'\t')
       | map(row -> tuple(row.matrix_name, file(row.matrix), file(row.sample_names)))
    fromMatrix(matrices, "${file(params.matrices_list).baseName}")
}

workflow fromBinaryMatrix {
    matrices = Channel.fromPath(params.binary_matrix)
       | map(it -> tuple("DHS_Binary", it, file(params.sample_names)))
    fromMatrix(matrices, "DHS_Binary")
}

// Start from CAV calling pval file
workflow fromAggregatedCavs {
    params.result_pval_file = "${params.outdir}/aggregated.${params.aggregation_key}.bed"
    pvals = Channel.fromPath(params.result_pval_file) 
        | split_cell_specific_aggregation
        | flatten()
        | filter { it.countLines() >= params.min_snps }
        | map(it -> tuple(it.baseName, it)) // group_id, aggregated_pval_file
    
    fromAnnotations(pvals, "CAVs.${params.aggregation_key}")
   
}