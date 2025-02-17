import configparser
import sys
import pandas as pd


def main(metadata, outdir):
    for _, row in metadata.iterrows():
        config = configparser.ConfigParser()
        prefix = row['prefix']
        base_path = f"{outdir}/mixing/{prefix}/{prefix}"
        # Add sections and settings
        config['METADATA'] = {
            'INDEX_ANNDATA': row['anndata_path']
        }

        config['NMF'] = {
            'PREFIX': prefix,
            'N_COMPONENTS': row['n_components'],
            'W': row['W'],
            'H': row['H'],
            'ANNDATA': row['anndata_path'],
            'PEAKS_MASK': row['peaks_mask'],
            'PEAK_WEIGHTS': row.get("peaks_weights", ""),
            'SAMPLE_WEIGHTS': row.get("samples_weights", ""),

            'PURE.50PR_ANNOTATION': f"{base_path}.pure.50pr.npy",
            'PURE.50PR_ORDER': f"{base_path}.pure.50pr.order.txt",
            'MIXING.80PR_ANNOTATION': f"{base_path}.mixing.80pr.npy",
            'MIXING.80PR_ORDER': f"{base_path}.mixing.80pr.order.txt",

            'TOP_SAMPLES': f"{outdir}/top_samples/{prefix}.top_samples.tsv",
            'DENSITY_TRACKS_META': f"{outdir}/top_samples/{prefix}.density_tracks_meta.tsv"
        }


        config['LDSC'] = {
            'Z_SCORE_SUMMARY.PURE.50pr': f"{outdir}/{prefix}.pure.50pr.ldsc_cell_types_results.tsv",
            'Z_SCORE_SUMMARY.MIXING.80pr': f"{outdir}/{prefix}.mixing.80pr.ldsc_cell_types_results.tsv"
        }

        config['MOTIF.ENRICHMENT'] = {
            'Z_SCORE_SUMMARY.PURE.50pr': f"{outdir}/{prefix}.pure.50pr.z_score_stats.tsv",
            'Z_SCORE_SUMMARY.MIXING.80pr': f"{outdir}/{prefix}.mixing.80pr.z_score_stats.tsv"
        }

        with open(f'{prefix}.matrix_meta.tsv', 'w') as f:
            f.write('\t'.join(
                [
                    'matrix_name',
                    'matrix',
                    'sample_names',
                    'anndata_path'
                ]
            ) + '\n')
            f.write('\t'.join(
                [
                    config['NMF']['PREFIX'] + '.pure.50pr', 
                    config['NMF']['PURE.50PR_ANNOTATION'],
                    config['NMF']['PURE.50PR_ORDER'],
                    config['NMF']['ANNDATA'],
                    config['NMF']['PEAKS_MASK']
                ]
            ) + '\n')
            f.write('\t'.join(
                [
                    config['NMF']['PREFIX'] + '.mixing.80pr', 
                    config['NMF']['MIXING.80PR_ANNOTATION'],
                    config['NMF']['MIXING.80PR_ORDER'],
                    config['NMF']['ANNDATA'],
                    config['NMF']['PEAKS_MASK']
                ]
            ) + '\n')

        # Write the configuration file
        with open(f'{prefix}.config.ini', 'w') as configfile:
            config.write(configfile)


if __name__ == '__main__':
    nmf_meta = pd.read_table(sys.argv[1])

    main(nmf_meta, sys.argv[2])
