"""
FASTQ Parser Module
===================
Parses FASTQ files using BioPython and exposes structured read data.
"""

from collections import namedtuple
from pathlib import Path
from Bio import SeqIO
import sys

# Structured representation of a single FASTQ read
FastqRead = namedtuple(
    "FastqRead",
    ["id", "sequence", "quality_scores"]
)


def parse_fastq(filepath):
    """
    Generator that yields FastqRead objects from a FASTQ file.
    """

    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"FASTQ file not found: {filepath}")

    for record in SeqIO.parse(str(filepath), "fastq"):

        yield FastqRead(
            id=record.id,
            sequence=str(record.seq),
            quality_scores=record.letter_annotations["phred_quality"],
        )


def count_reads(filepath):
    """
    Count total reads in a FASTQ file.
    """

    return sum(1 for _ in parse_fastq(filepath))


def get_read_lengths(filepath):
    """
    Return a list of read lengths.
    """

    return [len(read.sequence) for read in parse_fastq(filepath)]


if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Usage: python fastq_parser.py <fastq_file>")
        sys.exit(1)

    fq_path = sys.argv[1]

    print(f"\nParsing: {Path(fq_path).name}")

    total = 0

    for read in parse_fastq(fq_path):

        total += 1

        if total <= 3:

            avg_quality = (
                sum(read.quality_scores)
                / len(read.quality_scores)
            )

            print(
                f"Read {total}: "
                f"{read.id} | "
                f"len={len(read.sequence)} | "
                f"avg_qual={avg_quality:.1f}"
            )

    print(f"\nTotal reads: {total}")