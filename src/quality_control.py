"""
Quality Control Module
======================
Computes QC statistics for FASTQ files
and generates human-readable reports.
"""

from pathlib import Path
from datetime import datetime
import statistics

from Bio.SeqUtils import gc_fraction

from fastq_parser import parse_fastq


def get_gc_content(sequence):
    """
    Calculate GC content percentage.
    """

    return gc_fraction(sequence) * 100


def check_quality_thresholds(qualities):
    """
    Count Q20 and Q30 bases.
    """

    q20_count = sum(1 for q in qualities if q >= 20)
    q30_count = sum(1 for q in qualities if q >= 30)

    return q20_count, q30_count


def update_per_base_quality(qualities, sums, counts):
    """
    Accumulate per-base quality scores.
    """

    for i, q in enumerate(qualities):

        if i >= len(sums):
            sums.append(0.0)
            counts.append(0)

        sums[i] += q
        counts[i] += 1

    return sums, counts


def compute_qc_stats(filepath):
    """
    Compute QC statistics for one FASTQ file.
    """

    filepath = Path(filepath)

    read_lengths = []
    all_qualities = []

    total_gc_bases = 0.0
    total_bases = 0

    total_q20 = 0
    total_q30 = 0

    pos_sums = []
    pos_counts = []

    for read in parse_fastq(filepath):

        seq = read.sequence
        quals = read.quality_scores

        read_len = len(seq)

        read_lengths.append(read_len)

        total_bases += read_len

        all_qualities.extend(quals)

        total_gc_bases += (
            get_gc_content(seq) / 100
        ) * read_len

        q20, q30 = check_quality_thresholds(quals)

        total_q20 += q20
        total_q30 += q30

        pos_sums, pos_counts = update_per_base_quality(
            quals,
            pos_sums,
            pos_counts,
        )

    if total_bases == 0:
        return None

    per_base_quality = [

        round(s / c, 2)

        for s, c in zip(pos_sums, pos_counts)
    ]

    return {

        "filename": filepath.name,

        "total_reads": len(read_lengths),

        "avg_read_length": round(
            statistics.mean(read_lengths), 2
        ),

        "min_read_length": min(read_lengths),

        "max_read_length": max(read_lengths),

        "avg_quality_score": round(
            statistics.mean(all_qualities), 2
        ),

        "gc_content_percent": round(
            (total_gc_bases / total_bases) * 100,
            2
        ),

        "q20_percent": round(
            (total_q20 / total_bases) * 100,
            2
        ),

        "q30_percent": round(
            (total_q30 / total_bases) * 100,
            2
        ),

        "per_base_quality": per_base_quality,
    }


def generate_qc_report(stats, output_path):
    """
    Generate a human-readable QC report.
    """

    if stats is None:
        print("No reads found.")
        return

    output_path = Path(output_path)

    output_path.parent.mkdir(
        parents=True,
        exist_ok=True
    )

    lines = [

        "=" * 60,

        "             FASTQ Quality Control Report",

        f"Generated: "
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",

        "=" * 60,

        "\n[File Information]",

        f"File Name:          {stats['filename']}",

        f"Total Reads:        {stats['total_reads']:,}",

        f"Average Read Length:{stats['avg_read_length']} bp",

        f"Length Range:       "
        f"{stats['min_read_length']} - "
        f"{stats['max_read_length']} bp",

        "\n[Quality Metrics]",

        f"Average Quality:    "
        f"{stats['avg_quality_score']} (Phred)",

        f"GC Content:         "
        f"{stats['gc_content_percent']}%",

        f"Q20 Bases:          "
        f"{stats['q20_percent']}%",

        f"Q30 Bases:          "
        f"{stats['q30_percent']}%",

        "\n[Per-Base Quality (First 15 Positions)]",
    ]

    pbq = stats["per_base_quality"]

    for i, q in enumerate(pbq[:15]):

        lines.append(
            f"Position {i+1:>3}: {q}"
        )

    lines.append("\n" + "=" * 60)

    with open(output_path, "w", encoding="utf-8") as f:

        f.write("\n".join(lines) + "\n")

    print(f" Report saved: {output_path}")


def process_directory(input_dir, output_dir):
    """
    Run QC on all FASTQ files in a directory.
    """

    input_dir = Path(input_dir)
    output_dir = Path(output_dir)

    fastq_files = list(input_dir.glob("*.fastq"))

    if not fastq_files:
        print("No FASTQ files found.")
        return

    print(f"\nFound {len(fastq_files)} FASTQ files.\n")

    for fq_file in fastq_files:

        print(f"Processing: {fq_file.name}")

        stats = compute_qc_stats(fq_file)

        report_path = (
            output_dir /
            f"{fq_file.stem}_report.txt"
        )

        generate_qc_report(
            stats,
            report_path
        )

    print("\n QC completed for all files.")


if __name__ == "__main__":

    import sys

    if len(sys.argv) < 2:

        print(
            "Usage:\n"
            "python quality_control.py <fastq_file>\n"
            "OR\n"
            "python quality_control.py <directory>"
        )

        sys.exit(1)

    target = Path(sys.argv[1])

    results_dir = Path("../results/reports")

    if target.is_file():

        stats = compute_qc_stats(target)

        report_path = (
            results_dir /
            f"{target.stem}_report.txt"
        )

        generate_qc_report(
            stats,
            report_path
        )

    # Directory mode (automation)
    elif target.is_dir():

        process_directory(
            target,
            results_dir
        )

    else:

        print(f"Error: {target} not found")