from collections import Counter
from Bio import SeqIO
import pandas as pd
from pathlib import Path

def get_kmers(sequence, k=3):
    return [sequence[i:i+k] for i in range(len(sequence) - k + 1)]


def extract_kmer_features(fastq_file, k=3):
    kmers_counter = Counter()
    total_sequences = 0

    for record in SeqIO.parse(fastq_file, "fastq"):
        seq = str(record.seq)
        kmers = get_kmers(seq, k)
        kmers_counter.update(kmers)
        total_sequences += 1

    # normalize frequencies
    total_kmers = sum(kmers_counter.values())

    features = {kmer: count / total_kmers for kmer, count in kmers_counter.items()}

    features["sample"] = Path(fastq_file).stem

    return features


def run_kmer_extraction(input_folder, output_csv, k=3):
    results = []

    for file in Path(input_folder).glob("*.fastq"):
        features = extract_kmer_features(file, k)
        results.append(features)

    df = pd.DataFrame(results).fillna(0)
    df.to_csv(output_csv, index=False)

    print(f"K-mer features saved to {output_csv}")


if __name__ == "__main__":
    run_kmer_extraction(
        input_folder="D:/project_fastq/data/raw",
        output_csv="D:/project_fastq/kmer_features.csv",
        k=3
    )