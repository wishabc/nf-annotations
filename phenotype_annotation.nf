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
process run_ldsc {
    conda params.ldsc_conda
    publishDir "${params.outdir}/${prefix}/ldsc", pattern: "${name}.results"
    publishDir "${params.outdir}/${prefix}/ldsc_logs", pattern: "${name}.logs"
    publishDir "${params.outdir}/${prefix}/ldsc_logs", pattern: "${name}.part_delete"
    tag "${prefix}:${phen_id}"
    scratch true

    input:
        tuple val(phen_id), path(sumstats_file), val(prefix), path(ld_files)
    
    output:
        tuple val(phen_id), path("${name}.results"), emit: results
        tuple val(phen_id), path("${name}.*"), emit: all_data


    script:
    name = "${prefix}.${phen_id}"
    pref = "${prefix}."
    """
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

    input:
        paths ldsc_files
    
    output:
        path name

    script:
    name = "ldsc_result.tsv"
    """
    echo "group_name\tphenotype_id\t`head -1 ${ldsc_files[0]}`"" > ${name}
    echo "${ldsc_files}" | tr " " "\n" > filelist.txt
    while read line; do
        echo "`basename "\$line" .results | tr "." "\t"`\t`tail -1 \$line`" >> ${name}
    done < filelist.txt
    """
}


workflow LDSC {
    take:
        ld_data
    main:
        ldsc_res = Channel.fromPath(params.phenotypes_meta)
            .splitCsv(header:true, sep:'\t')
            .map(row -> tuple(row.phen_id, file(row.sumstats_file)))
            | combine(ld_data)
            | run_ldsc
        out = ldsc_res.results
            | map(it -> it[1])
            | collect()
            | collect_ldsc_results
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
        out = LDSC(ldsc_data)
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

// defunc
workflow test {
        Channel.fromPath("/net/seq/data2/projects/sabramov/ENCODE4/dnase-annotations/LDSC.clusters/output/*/ldsc/*.results")
            | map(it -> file(it))
            | collect()
            | collect_ldsc_results
}

workflow annotateWithPheno {
    out =  Channel.fromPath("${params.pval_file_dir}/*.bed")
        | map(it -> file(it))
        | collect(sort: true)
        | filterUniqVariants
        | annotate_with_phenotypes
}
