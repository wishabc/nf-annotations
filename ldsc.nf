#!/usr/bin/env nextflow
// include { filterTestedVariants } from "./main"

params.conda = "$moduleDir/environment.yml"

// TODO wrap in apptainer
process calc_ld {
    publishDir "${params.outdir}/ldsc/l2_logs", pattern: "${prefix}.log"
    publishDir "${params.outdir}/ldsc/l2", pattern: "${prefix}.l2.*"
    publishDir "${params.outdir}/ldsc/l2", pattern: "${annotation_file}"
    
    tag "chr${chrom}:${group_id}"
    scratch true
    conda params.ldsc_conda

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
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file), val(baseline_ld), path("data_files/*")
    
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
        --ref-ld-chr ${baseline_ld} \
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
    publishDir "${params.outdir}/ldsc/ldsc_coefs_${prefix}", pattern: "${name}.results"
    publishDir "${params.outdir}/ldsc/ldsc_logs/${prefix}", pattern: "${name}.log"
    tag "${prefix}:${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file), val(baseline_ld), val(prefix), path(ld_files)
    
    output:
        tuple val(prefix), val(phen_id), path("${name}.results"), path("${name}.log")

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
    """
}


process collect_ldsc_results {
    scratch true
    publishDir "${params.outdir}"

    input:
        val ldsc_files
    
    output:
        path name

    script:
    name = "ldsc_enrichments_results.tsv"
    """
    echo "${ldsc_files.join('\n')}" \
        | grep '.results' \
         > filelist.txt

    # copy header from the first file
    head -1 filelist.txt \
        | xargs head -1 \
        | xargs -I % echo "group_name\tphenotype_id\t%" > result.txt
    
    # Aggregate the data
    echo -e "h^2\th^2_err" > h2.stats
    while read line; do
        fname="\${line%.*}"
        grep "Total Observed scale h2" \${fname}.log \
            | sed 's/[:(]/\t/g' \
            | sed 's/)//g' \
            | awk -F'\t' -v OFS='\t' '{print \$2,\$3}' >> h2.stats

        basename \$fname \
            | sed "s/\\./\t/" \
            | xargs -I % echo "%\t`tail -1 "\$line"`" >> result.txt
    done < filelist.txt
    paste result.txt h2.stats > ${name}
    """
}

// Make annotation workflows
process filter_cavs {
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
        tuple val(chrom), path(custom_annotation)

    output:
        tuple val(group_id), val(chrom), path(name)
    
    script:
    group_id = "${custom_annotation.simpleName}"
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
                name: "ldsc_${params.aggregation_key}_results.tsv"
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
            | collect(sort: true, flat: true)
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
            | map(row -> tuple(row.phen_id, file(row.munge_sumstats_file), params.baseline_ld))
            | filter { it[1].exists() }

        ldsc_annotations = Channel.of(1..22)
            | combine(annotations)
            | make_ldsc_annotation // group_id, chrom, annotation

        ld_data = ldsc_annotations
            | calc_ld //  group_id, chrom, ld, ld_log, annotation
            | map(it -> tuple(it[0], [it[2], it[4]].flatten()))
            | groupTuple(size: 22)
            | map(it -> tuple(it[0], it[1].flatten()))

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
    sumstats_files = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, file(row.munge_sumstats_file), "baseline."))
        | filter { it[1].exists() } // phen_id, sumstats, baseline

    ld_data = Channel.of(1..22)
        | map(it -> tuple('baseline', it, file("${params.base_ann_path}${it}.annot.gz", checkIfExists: true)))
        | calc_ld //  group_id, chrom, ld, ld_log, annotation
        | flatMap(it -> [it[2], it[4]])
        | collect(sort: true)
        | map(it -> tuple('baseline', it))
    
    // phen_id, sumstats_file, baseline_ld, val(prefix), path(ld_files)
    LDSC(ld_data, sumstats_files)
        
}

workflow fromPvalFiles {
    params.result_pval_file = "${params.outdir}/aggregated.${params.aggregation_key}.bed"
    Channel.fromPath(params.result_pval_file) 
        | filter_cavs
        | flatten()
        | filter { it.countLines() >= params.min_snps }
        | fromAnnotations
   
}

workflow {
    custom_annotations = Channel.fromPath("${params.annotations_dir}/*.bed") 
        | fromAnnotations
}
