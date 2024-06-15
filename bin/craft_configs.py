import configparser
import sys
import pandas as pd

def main(metadata, outdir):
    for _, row in metadata.iterrows():
        config = configparser.ConfigParser()
        prefix = row['prefix']
        base_path = f"{outdir}/mixing_annotations/{prefix}/{prefix}"
        # Add sections and settings
        config['NMF'] = {
            'PREFIX': prefix,
            'N_COMPONENTS': prefix.split('.')[-1],
            'W': row['W'],
            'H': row['H'],
            'PEAK_WEIGHTS': row['peak_weights'],
            'SAMPLE_WEIGHTS': row['sample_weights'],
            'INDEX': row['masterlist'],
            'SAMPLES_ORDER': row['samples_order'],

            'PURE.50PR_ANNOTATION': f"{base_path}.pure.50pr.npy",
            'PURE.50PR_ORDER': f"{base_path}.pure.order.txt",
            'MIXING.80PR_ANNOTATION': f"{base_path}.mixing.80pr.npy",
            'MIXING.80PR_ORDER': f"{base_path}.mixing.80pr.order.txt",

            'TOP_SAMPLES': f"{outdir}/{prefix}.top_samples.tsv",
            'DENSITY_TRACKS_META': f"{outdir}/{prefix}.density_tracks_meta.tsv"
        }

        config['LDSC'] = {
            'Z_SCORE_SUMMARY.PURE': f"{outdir}/{prefix}.pure.ldsc_cell_types_results.tsv",
            'Z_SCORE_SUMMARY.MIXING': f"{outdir}/{prefix}.mixing.ldsc_cell_types_results.tsv"
        }

        config['MOTIF.ENRICHMENT'] = {
            'Z_SCORE_SUMMARY.PURE': f"{outdir}/{prefix}.pure.z_score_stats.tsv",
            'Z_SCORE_SUMMARY.MIXING': f"{outdir}/{prefix}.mixing.z_score_stats.tsv"
        }

        # Write the configuration file
        with open(f'{prefix}.config.ini', 'w') as configfile:
            config.write(configfile)


if __name__ == '__main__':
    nmf_meta = pd.read_table(sys.argv[1])
