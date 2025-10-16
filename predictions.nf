
params.fasta_file = "/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa"
params.sample_genotype_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-wasp.v5/output/meta+sample_ids.tsv"
params.genotype_file = "/net/seq/data2/projects/sabramov/ENCODE4/dnase-wasp.v5/output/all_variants_stats.bed.gz"


process predict {
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir "${params.outdir}/predictions"
    label "gpu"

    input:
        tuple val(prefix), path(embeddings_file), path(dhs_dataset), path(checkpoint), val(model_type)
    
    output:
        tuple val(prefix), path(name)
    
    script:
    name = "${prefix}.npy"
    """
    python3 $moduleDir/bin/dhs_prediction/predict_DHS_model.py \
        ${dhs_dataset} \
        --embeddings_file ${embeddings_file} \
        --checkpoint ${checkpoint} \
        --fasta_file ${params.fasta_file} \
        --sample_genotype_file ${params.sample_genotype_file} \
        --genotype_file ${params.genotype_file} \
        --num_workers ${task.cpus} \
        --model_type ${model_type} \
        --output ${name} 
    """
}


process annotate_with_predictions {
    conda "/home/sabramov/miniconda3/envs/pytorch"
    publishDir "${params.outdir}/"

    output:
        path name

    script:
    name = "${file(params.samples_file).baseName}.annotated_with_predictions.tsv"
    """
    python3 $moduleDir/bin/dhs_prediction/annotate_meta.py \
        ${params.samples_file} \
        ${params.outdir}/predictions \
        ${name}
    """
}

workflow {
    Channel.fromPath(params.samples_file)
        | splitCsv(header:true, sep:'\t')
        | map(row -> tuple(row.prefix, file(row.embeddings_file), file(row.dhs_dataset), file(row.checkpoint), row.model_type))
        | predict
    
    annotate_with_predictions()
}