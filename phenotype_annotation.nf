#!/usr/bin/env nextflow
include { filterTestedVariants } from "./motif_enrichment"

params.conda = "$moduleDir/environment.yml"

is_baseline = false

process filter_cavs {

    input:
        path pval_file

    output:
        path "*.bed"
    
    script:
    """
    cat ${pval_file} | grep -v '#' | awk '\$NF <= ${params.fdr_tr} {print> \$19".bed"}'
    """
}

process make_ldsc_annotation {
    conda params.conda
    tag "chr${chrom}:${annotation.simpleName}"
    scratch true

    input:
        tuple val(chrom), path(annotation)

    output:
        tuple val(chrom), path(name)
    
    script:
    suffix = "${chrom}.annot.gz"
    baseannotation = "${params.base_ann_path}${suffix}"
    name = "${annotation.simpleName}.${suffix}"
    """
    cat ${annotation} |  sed -e "s/^chr//" | sort-bed -  > annot_numchr.bed
    echo "CHR\tBP\tSNP\tCM\t${annotation.simpleName}" | gzip > ${name}
    zcat ${baseannotation} \
        | awk -v OFS='\t' -F'\t' '(NR > 1) { print \$1,\$2-1,\$2,\$3,\$4 }'\
        | bedtools intersect -wa -c -a stdin -b annot_numchr.bed \
        | awk -v OFS='\t' '{ print \$1,\$3,\$4,\$5,\$6}' | gzip >> ${name}
    """
}

process munge_sumstats {
    conda params.ldsc_conda
    tag "${phen_id}"
    publishDir "${params.outdir}/munge_sumstats", pattern: "${prefix}.sumstats.gz"
    publishDir "${params.outdir}/munge_sumstats_logs", pattern: "${prefix}.log"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file)

    output:
        tuple val(phen_id), path("${prefix}.sumstats.gz"), emit: result
        path("${prefix}*"), emit: all
    
    script:
    prefix = "UKBB_${sumstats_file.simpleName}"
    """
    python ${params.ldsc_scripts_path}/munge_sumstats.py \
        --sumstats ${sumstats_file} \
        --merge-alleles ${params.tested_snps} \
        --out ${prefix}
    """
}

// TODO wrap in apptainer
process calc_ld {
    publishDir "${outdir}/logs", pattern: "${name}.log", enabled: !is_baseline
    publishDir "${outdir}", pattern: "${name}.l2.ldscore.gz"
    publishDir "${outdir}", pattern: "${name}.l2.M*"
    publishDir "${outdir}", pattern: "${annotation_file}", enabled: !is_baseline
    tag "chr${chrom}:${annotation_file.simpleName}"
    scratch true
    conda params.ldsc_conda

    input:
        tuple val(chrom), path(annotation_file)
    
    output:
        tuple val(annotation_file.simpleName), path("${name}.l2.*"), path(annotation_file), emit: result
        path("${name}.log"), emit: logs
    
    script:
    if (is_baseline) {
        outdir = file(params.base_ann_path).parent
    } else {
        outdir = "${params.outdir}/l2/${annotation_file.simpleName}"
    }
    name = "${annotation_file.simpleName}.${chrom}"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    
    # TODO: Check if --print-snps parameter is needed
    ${params.ldsc_scripts_path}/ldsc.py \
        --print-snps ${params.tested_snps} \
        --ld-wind-cm 1.0 \
        --out ${name} \
        --bfile ${params.gtfiles}${chrom} \
        --annot ${annotation_file} \
        --l2
    """
}
// TODO wrap in apptainer
process run_ldsc_cell_types {
    conda params.ldsc_conda
    publishDir "${params.outdir}/ldsc", pattern: "${name}.cell_type_results.txt"
    publishDir "${params.outdir}/ldsc_logs", pattern: "${name}.log"
    tag "${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file)
        path "data_files/*"
    
    output:
        path "${name}.cell_type_results.txt", emit: results
        path "${name}.log", emit: logs
        tuple val(phen_id), path("${name}.*"), emit: all_data


    script:
    name = "${phen_id}"
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}

    ls -1 data_files | cut -d"." -f 1 | sort \
        | uniq | awk -v OFS='\\t' '{ print \$1,"data_files/"\$1"." }' > per_sample.ldcts
    
    ${params.ldsc_scripts_path}/ldsc.py \
        --h2-cts ${sumstats_file} \
        --ref-ld-chr ${params.base_ann_path} \
        --ref-ld-chr-cts per_sample.ldcts \
        --frqfile-chr ${params.frqfiles} \
        --w-ld-chr ${params.weights} \
        --overlap-annot \
        --print-coefficients \
        --print-delete-vals \
        --out ${name}
    """
}

process run_ldsc_single_sample {
    conda params.ldsc_conda
    publishDir "${params.outdir}/${prefix}/ldsc", pattern: "${name}.results"
    publishDir "${params.outdir}/${prefix}/ldsc_logs", pattern: "${name}.log"
    tag "${prefix}:${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file), val(prefix), path(ld_files)
    
    output:
        path "${name}.results", emit: results
        path "${name}.log", emit: logs
        tuple val(prefix), val(phen_id), path("${name}.*"), emit: all_data


    script:
    name = "${prefix}.${phen_id}"
    pref = "${prefix}."
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    ${params.ldsc_scripts_path}/ldsc.py \
        --h2 ${sumstats_file} \
        --ref-ld-chr ${params.base_ann_path},${pref} \
        --frqfile-chr ${params.frqfiles} \
        --w-ld-chr ${params.weights} \
        --overlap-annot \
        --print-coefficients \
        --print-delete-vals \
        --out ${name}
    """
}


process collect_ldsc_results {
    scratch true
    publishDir "${params.outdir}"

    input:
        // expected to be more than one file
        val ldsc_files
    
    output:
        path name
        

    script:
    name = "ldsc_result.tsv"
    """
    # copy header from the first file
    head -1 ${ldsc_files[0]} | xargs -I % echo "group_name\tphenotype_id\t%" > ${name}
    
    # Aggregate the data
    echo '${ldsc_files.join("\n")}' > filelist.txt
    while read line; do
        basename "\$line" .results | tr "." "\t" | xargs -I % echo "%\t`tail -1 "\$line"`" >> ${name}
    done < filelist.txt
    """
}

workflow LDSCcellTypes {
    take:
        ld_data
    main:
        sumstats = Channel.fromPath(params.phenotypes_meta)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.phen_id, file(row.sumstats_file)))
            | filter { it[1].exists() }
            // | munge_sumstats
        
        ldsc_res = run_ldsc_cell_types(
            sumstats, 
            ld_data.map(it -> it[1]).collect(sort: true))

        out = ldsc_res.results.collect(sort: true)
        // out = collect_ldsc_results(l) FIXME
    emit:
        out
}


workflow LDSC {
    take:
        ld_data
    main:
        ldsc_res = Channel.fromPath(params.phenotypes_meta)
            | splitCsv(header:true, sep:'\t')
            | map(row -> tuple(row.phen_id, file(row.sumstats_file)))
            | filter { it[1].exists() }
            | combine(ld_data)
            | run_ldsc_single_sample

        out = ldsc_res.results
            | collect(sort: true)
            | collect_ldsc_results
    emit:
        out
}

workflow calcBaseline {
    is_baseline = true
    data = Channel.of(1..22)
        | map(it -> tuple(it, file("${params.base_ann_path}${it}.annot.gz", checkIfExists: true)))
        | calc_ld
}

workflow fromAnnotations {
    take:
        annotations
    main:
        ldsc_data = Channel.of(1..22)
            | combine(annotations)
            | make_ldsc_annotation
            | calc_ld
            | results
            | map(it -> tuple(it[0], [it[1], it[2]].flatten()))
            | groupTuple(size: 22)
            | map(
                it -> tuple(it[0], it[1].flatten())
            )
        if (params.by_cell_type) {
            out = LDSCcellTypes(ldsc_data)
        } else {
            out = LDSC(ldsc_data)    
        }
        
    emit:
        out
}

workflow fromPvalFiles {
    params.fdr_tr = 0.05
    Channel.fromPath(params.pval_file) 
        | filter_cavs
        | flatten()
       // | fromAnnotations
   
}
workflow {
    custom_annotations = Channel.fromPath("${params.annotations_dir}/*.bed") 
        | map(it -> file(it))
        | filterTestedVariants
        | fromAnnotations
}

// Defunc
workflow mergeResults {
    Channel.fromPath(
        "/net/seq/data2/projects/sabramov/ENCODE4/dnase-annotations/LDSC.clusters/output/*/ldsc/*.results"
    )   
        | collect(sort: true, flat: true) 
        | collect_ldsc_results
}
