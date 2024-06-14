include { LDSC; calc_ld } from './ldsc'


workflow {
    sumstats_files = Channel.fromPath(params.phenotypes_meta)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.phen_id, file(row.munge_sumstats_file), "baseline."))
        | filter { it[1].exists() } // phen_id, sumstats, baseline

    ld_data = Channel.of(1..22)
        | map(it -> tuple('baseline', 'baseline', it, file("${params.base_ann_path}${it}.annot.gz", checkIfExists: true)))
        | calc_ld //  matrix_name, group_id, chrom, ld, ld_log, annotation
        | flatMap(it -> [it[3], it[5]])
        | collect(sort: true)
        | map(it -> tuple('baseline', 'baseline', it))
    
    // phen_id, sumstats_file, baseline_ld, val(prefix), path(ld_files)
    LDSC(ld_data, sumstats_files)
        
}