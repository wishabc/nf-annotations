import subprocess
import pandas as pd
import sys

def check_column_flag(columns, column_name, flag_value, default_value):
    
    if column_name in columns:
        return flag_value
    else:
        return default_value


def main(sumstats_file, script_path, tested_snps, n_samples, prefix):
    df = pd.read_table(sumstats_file, nrows=1)
    if 'beta' in df.columns and df['beta'].dropna().shape[0] > 0:
        sumstats_flag = '--signed-sumstats beta'
        ignore = 'odds_ratio '
    else:
        sumstats_flag = '--signed-sumstats odds_ratio'
        ignore = 'beta'
    effect_allele_frequency_flag = check_column_flag(sumstats_file, "effect_allele_frequency", "--frq effect_allele_frequency", "")
    snp_flag = check_column_flag(sumstats_file, "rs_id", "--snp rs_id", "--snp variant_id")

    cmd = [
        'python2', script_path,
        '--sumstats', sumstats_file,
        '--merge-alleles', tested_snps,
        '--a1', 'effect_allele',
        '--a2', 'other_allele',
        '--N', n_samples,
        sumstats_flag,
        snp_flag,
        effect_allele_frequency_flag,
        '--out', prefix,
        '--ignore', ignore
    ]

    # Remove empty strings from cmd list
    cmd = [arg for arg in cmd if arg]

    subprocess.call(cmd)


if __name__ == "__main__":
    sumstats_file = sys.argv[1]
    script_path = sys.argv[2]
    tested_snps = sys.argv[3]
    n_samples = sys.argv[4]
    prefix = sys.argv[5]
    
    main(sumstats_file, script_path, tested_snps, n_samples, prefix)
