#!/usr/bin/env nextflow
include { filterUniqVariants } from "./motif_enrichment"

params.conda = "$moduleDir/environment.yml"

params.phenotypes_data = "/home/sabramov/phenotypes_data"

is_baseline = false
// Annotates with pheWAS, clinvar, finemapping, grasp, ebi-gwas phenotypes
process annotate_with_phenotypes {
    conda params.conda
    publishDir "${params.outdir}/phenotypes"

    input:
        path pval_file

    output:
        path name

    script:
    name = "phenotypes_ann.bed"
    """
    python3 $moduleDir/bin/annotate_with_phenotypes.py ${params.phenotypes_data} ${pval_file} ${name}
    """
}

process filter_cavs {
    tag "${prefix}"

    input:
        path pval_file

    output:
        path name
    
    script:
    prefix = pval_file.simpleName
    name = "${prefix}.fdr${params.fdr_tr}.bed"
    """
    head -1 ${pval_file} > ${name}
    cat ${pval_file} | awk '((\$NF <= ${params.fdr_tr}) && (NR>1)) {print}' >> ${name}
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
    conda params.conda
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
    baseannotation = "${params.base_ann_path}${suffix}"
    prefix = "UKBB_${sumstats_file.simpleName}"
    """
    python ${params.ldsc_scripts_path}/munge_sumstats.py \
        --sumstats ${sumstats_file} \
        --merge-alleles ${params.tested_snps} \
        --out ${name}
    """
}

// TODO wrap in apptainer
process calc_ld {
    publishDir "${outdir}/l2_logs", pattern: "${name}.log", enabled: !is_baseline
    publishDir "${outdir}", pattern: "${name}.l2.ldscore.gz"
    publishDir "${outdir}", pattern: "${name}.l2.M*"
    publishDir "${outdir}/l2", pattern: "${annotation_file}", enabled: !is_baseline
    tag "chr${chrom}:${annotation_file.simpleName}"
    scratch true
    conda params.ldsc_conda

    input:
        tuple val(chrom), path(annotation_file)
    
    output:
        tuple val(annotation_file.simpleName), path("${name}*"), path(annotation_file)
    
    script:
    if (is_baseline) {
        outdir = file(params.base_ann_path).parent
        name = "${annotation_file.simpleName}.${chrom}"
    } else {
        outdir = "${params.outdir}/${annotation_file.simpleName}"
        name = "l2/${annotation_file.simpleName}.${chrom}"
    }
    """
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export GOTO_NUM_THREADS=${task.cpus}
    export OMP_NUM_THREADS=${task.cpus}
    
    mkdir l2
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
    export OPENBLAS_NUM_THREADS="${task.cpus}"
    export GOTO_NUM_THREADS="${task.cpus}"
    export OMP_NUM_THREADS="${task.cpus}"


    
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
        path ldsc_files
    
    output:
        tuple path(name), path(pr_h2)
        

    script:
    name = "ldsc_result.tsv"
    """
    head -1 ${ldsc_files[0]} | xargs -I % echo "group_name\tphenotype_id\t%" > ${name}
    echo '${ldsc_files}' | tr ' ' '\n' > filelist.txt
    while read line; do
        echo "\$line `basename "\$line" .results | tr "." "\t"`"
        tail -1 "\$line" > ann.txt
        basename "\$line" .results | tr "." "\t" | xargs -I % echo "%\t`cat ann.txt`" >> ${name}
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
            | munge_sumstats
        
        ldsc_res = run_ldsc_cell_types(sumstats.result, ld_data.map(it -> it[1]).collect(sort: true))

        l = ldsc_res.results.collect(sort: true)
        out = collect_ldsc_results(l)
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

        l = ldsc_res.results.collect(sort: true)
        out = collect_ldsc_results(l)
    emit:
        out
}

workflow calcBaseline {
    data = Channel.of(1..22).map(
        it -> tuple(it, file("${params.base_ann_path}${it}.annot.gz", checkIfExists: true))
    )
    is_baseline = true
    calc_ld(data)
}

workflow fromAnnotations {
    take:
        annotations
    main:
        data = Channel.of(1..22).combine(annotations)
        anns = make_ldsc_annotation(data) 
        lds = calc_ld(anns)
        ldsc_data = lds.map(it -> tuple(it[0], [it[1], it[2]].flatten()))
            .groupTuple(size: 22)
            .map(
                it -> tuple(it[0], it[1].flatten())
            )
        //out = LDSC(ldsc_data)
        out = LDSCcellTypes(ldsc_data)
    emit:
        out
}

workflow fromPvalFiles {
    params.fdr_tr = 0.05
    Channel.fromPath("${params.pval_file_dir}/*.bed") 
        | map(it -> file(it))
        | filter_cavs
        | fromAnnotations
   
}
workflow {
    custom_annotations = Channel.fromPath("${params.annotations_dir}/*.bed") 
        | map(it -> file(it))
        | filterUniqVariants
        | fromAnnotations
}

workflow annotateWithPheno {
    out =  Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
        | collect(sort: true)
        | filterUniqVariants
        | annotate_with_phenotypes
}

