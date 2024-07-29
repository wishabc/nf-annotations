import subprocess
import pandas as pd
import sys

def check_column_flag(columns, column_names, flag_value, default_value):
    for column_name in column_names:
        if column_name in columns:
            return [flag_value[0], column_name], ','.join([x for x in column_names if x != column_name] + [default_value[1]])
    return default_value, ""


def main(sumstats_file, script_path, tested_snps, n_samples, prefix):
    df = pd.read_csv(sumstats_file, sep='\t')
    if 'beta' in df.columns and df['beta'].dropna().shape[0] > 0:
        sumstats_flag = ['--signed-sumstats', 'beta,0']
        ignore = 'odds_ratio'
    else:
        sumstats_flag = ['--signed-sumstats', 'odds_ratio,1']
        ignore = 'beta'
    if not 'other_allele' in df.columns or len(df['other_allele'].dropna()) == 0:
        print('Other allele is not present! Exiting...')
        exit(5)
    
    if df['p_value'].notna().sum() != len(df):
        print('P-values are not present for all variants! Exiting...')
        exit(5)

    effect_allele_frequency_flag, _ = check_column_flag(df.columns, ["effect_allele_frequency"], ["--frq", "effect_allele_frequency"], ["", ""])
    snp_flag, to_ignore = check_column_flag(df.columns, ["rs_id", 'rsid'], ["--snp", ""], ["--snp", "variant_id"])

    ignore = ','.join([ignore, to_ignore])
    cmd = [
        'python2', script_path,
        '--sumstats', sumstats_file,
        '--merge-alleles', tested_snps,
        '--a1', 'effect_allele',
        '--a2', 'other_allele',
        '--N', n_samples,
        '--out', prefix,
        '--ignore', ignore
    ]
    cmd += sumstats_flag
    cmd += effect_allele_frequency_flag
    cmd += snp_flag

    # Remove empty strings from cmd list
    cmd = [arg for arg in cmd if arg]

    exit_code = subprocess.call(cmd)
    if exit_code != 0:
        raise ValueError("Failed to run command:", cmd, exit_code)


if __name__ == "__main__":
    sumstats_file = sys.argv[1]
    script_path = sys.argv[2]
    tested_snps = sys.argv[3]
    n_samples = sys.argv[4]
    prefix = sys.argv[5]
    
    main(sumstats_file, script_path, tested_snps, n_samples, prefix)
