'''
simulate_null_phenotypes.py

This script will simulate null traits for testing SAIGE-GPU and estimating 
genotype correlations from beta correlations. 

***This needs to be run in the simulate_pheno conda environment ***
'''

import argparse
import pandas as pd
import numpy as np
from joblib import Parallel, delayed, cpu_count

def simulate_phenotypes_chunk(chunk, num_pheno, pheno_type, param, seed):
    rng = np.random.default_rng(seed)
    n = len(chunk)
    if pheno_type == 'binary':
        pheno_matrix = rng.binomial(1, param, (n, num_pheno))
    else:
        pheno_matrix = rng.normal(0, param, (n, num_pheno))
    pheno_df = pd.DataFrame(pheno_matrix, columns=[f'Pheno{i+1}' for i in range(num_pheno)])
    return pd.concat([chunk.reset_index(drop=True), pheno_df], axis=1)

def main():
    parser = argparse.ArgumentParser(description='Simulate phenotypes for genomic data (parallelized).')
    parser.add_argument('-i', '--input', required=True, help='Input file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    parser.add_argument('-n', '--num-pheno', type=int, default=1, help='Number of phenotypes to simulate (default: 1)')
    parser.add_argument('-t', '--type', choices=['binary', 'continuous'], default='continuous',
                        help='Phenotype type (default: continuous)')
    parser.add_argument('-p', '--param', type=float, help='Simulation parameter: case probability (binary) or SD (continuous)')
    parser.add_argument('--n-jobs', type=int, default=-1, help='Number of parallel jobs (-1: all cores)')
    args = parser.parse_args()

    if args.param is None:
        args.param = 0.5 if args.type == 'binary' else 1.0

    if args.type == 'binary' and not (0 <= args.param <= 1):
        raise ValueError("Case probability must be between 0 and 1")
    if args.type == 'continuous' and args.param <= 0:
        raise ValueError("Standard deviation must be positive")

    # Read input data
    df = pd.read_csv(args.input, sep=r'\s+')
    n_jobs = args.n_jobs if args.n_jobs != 0 else 1
    if n_jobs == -1:
        n_jobs = cpu_count()

    # Split data into chunks for each core
    chunks = np.array_split(df, n_jobs)
    # Use a unique seed for each chunk for reproducibility
    seeds = [42 + i for i in range(n_jobs)]

    # Parallel phenotype simulation
    results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(simulate_phenotypes_chunk)(chunk, args.num_pheno, args.type, args.param, seed)
        for chunk, seed in zip(chunks, seeds)
    )

    # Combine results and write to output
    final_df = pd.concat(results, ignore_index=True)
    final_df.to_csv(args.output, sep='\t', index=False, float_format='%.6f')

if __name__ == '__main__':
    main()
