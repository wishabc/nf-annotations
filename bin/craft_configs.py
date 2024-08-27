import configparser
import sys
import pandas as pd


def main(metadata, samples_meta_path, outdir):
    for _, row in metadata.iterrows():
        config = configparser.ConfigParser()
        prefix = f"{row['prefix']}.{row['n_components']}"
        base_path = f"{outdir}/mixing/{prefix}/{prefix}"
        # Add sections and settings
        config['METADATA'] = {
            'SAMPLES_META': samples_meta_path
        }

        config['NMF'] = {
            'PREFIX': prefix,
            'N_COMPONENTS': row['n_components'],
            'W': row['W'],
            'H': row['H'],
            'PEAK_WEIGHTS': row['peaks_weights'],
            'SAMPLE_WEIGHTS': row['samples_weights'],
            'INDEX': row['dhs_meta'],
            'SAMPLES_ORDER': row['sample_names'],

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
                    'dhs_coordinates'
                ]
            ) + '\n')
            f.write('\t'.join(
                [
                    config['NMF']['PREFIX'] + '.pure.50pr', 
                    config['NMF']['PURE.50PR_ANNOTATION'],
                    config['NMF']['PURE.50PR_ORDER'],
                    config['NMF']['INDEX']
                ]
            ) + '\n')
            f.write('\t'.join(
                [
                    config['NMF']['PREFIX'] + '.mixing.80pr', 
                    config['NMF']['MIXING.80PR_ANNOTATION'],
                    config['NMF']['MIXING.80PR_ORDER'],
                    config['NMF']['INDEX']
                ]
            ) + '\n')

        # Write the configuration file
        with open(f'{prefix}.config.ini', 'w') as configfile:
            config.write(configfile)


if __name__ == '__main__':
    nmf_meta = pd.read_table(sys.argv[1])

    main(nmf_meta, sys.argv[2], sys.argv[3])
