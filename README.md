# Stroke miRNA Classification Dataset Preparation

## Project Purpose

This project analyzes miRNA sequencing (miRNA-Seq) data from peripheral blood plasma to study microRNA expression changes associated with ischemic stroke.

The dataset compares two groups:

* Stroke patients
* Healthy controls

The goal is to explore circulating microRNAs that may serve as potential biomarkers for stroke diagnosis or biological insight into cerebrovascular injury.

---

## Why This Dataset Is Useful

Circulating miRNAs found in blood plasma are promising non-invasive biomarkers.
They can reflect disease-specific gene regulation changes and may help distinguish stroke-related molecular patterns using accessible blood samples.

---

## Quick Start Guide

Create a Python virtual environment and install the required packages.
python -m venv bio_env
bio_env\Scripts\activate
pip install -r requirements.txt

---

## Dataset Sources

Database: NCBI Sequence Read Archive (SRA)

Study: Stroke-induced miRNA expression in peripheral blood plasma

SRA Study Link:
https://www.ncbi.nlm.nih.gov/sra/?term=SRP435442

### Selected Runs
SRR24389252
SRR24389265
SRR24389268
SRR24389279
SRR24389283
SRR24389284

---

## Dataset Description

The dataset contains miRNA-Seq data generated using Illumina HiSeq 2500 from plasma samples of stroke patients and healthy individuals.
It enables analysis of miRNA expression patterns associated with ischemic stroke.
