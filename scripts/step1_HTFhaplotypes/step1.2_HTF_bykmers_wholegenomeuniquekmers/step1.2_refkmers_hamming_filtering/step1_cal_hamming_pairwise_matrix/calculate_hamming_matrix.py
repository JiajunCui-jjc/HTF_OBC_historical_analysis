#!/usr/bin/env python3

import os
import sys
import pandas as pd
from itertools import product
from multiprocessing import Pool, cpu_count

def hamming(seq1, seq2):
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))

def process_pair(args):
    query_ref, query_kmer, target_ref, target_kmer = args
    dist = hamming(query_kmer, target_kmer)
    return (query_ref, target_ref, query_kmer, target_kmer, dist)

def load_kmers(file):
    with open(file) as f:
        return [line.strip() for line in f if line.strip()]

def main(input_dir, output_file, nthreads):
    files = [f for f in os.listdir(input_dir) if f.endswith("_unique.txt")]
    refs = {os.path.splitext(f)[0].replace("_unique",""): load_kmers(os.path.join(input_dir,f)) for f in files}

    print(f"Loaded {len(refs)} references:")
    for ref, kmers in refs.items():
        print(f"  {ref}: {len(kmers)} kmers")

    all_tasks = []
    for query_ref, query_kmers in refs.items():
        for target_ref, target_kmers in refs.items():
            if query_ref != target_ref:
                for qkmer, tkmer in product(query_kmers, target_kmers):
                    all_tasks.append( (query_ref, qkmer, target_ref, tkmer) )

    print(f"Total pairwise comparisons: {len(all_tasks):,}")

    with Pool(processes=nthreads) as pool:
        results = pool.map(process_pair, all_tasks)

    df = pd.DataFrame(results, columns=["QueryRef","TargetRef","QueryKmer","TargetKmer","HammingDist"])
    df.to_csv(output_file, sep='\t', index=False)

    print(f"✅ Saved output to {output_file}")

if __name__ == "__main__":
    input_dir = sys.argv[1] if len(sys.argv)>1 else "."
    output_file = sys.argv[2] if len(sys.argv)>2 else "pairwise_hamming_distance_longformat.tsv"
    nthreads = int(sys.argv[3]) if len(sys.argv)>3 else cpu_count()
    main(input_dir, output_file, nthreads)
