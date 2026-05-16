import json
import pandas as pd
from pathlib import Path

def extract_fastp_features(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)

    summary = data.get("summary", {})
    before_filtering = summary.get("before_filtering", {})
    after_filtering = summary.get("after_filtering", {})

    features = {
        "total_reads_before": before_filtering.get("total_reads", 0),
        "total_reads_after": after_filtering.get("total_reads", 0),

        "q20_rate_before": before_filtering.get("q20_rate", 0),
        "q20_rate_after": after_filtering.get("q20_rate", 0),

        "q30_rate_before": before_filtering.get("q30_rate", 0),
        "q30_rate_after": after_filtering.get("q30_rate", 0),

        "gc_content_before": before_filtering.get("gc_content", 0),
        "gc_content_after": after_filtering.get("gc_content", 0),
    }

    return features


def run_feature_extraction(input_folder, output_csv):
    results = []

    for json_file in Path(input_folder).glob("*.json"):
        features = extract_fastp_features(json_file)
        features["sample"] = json_file.stem
        results.append(features)

    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"Saved features to {output_csv}")


if __name__ == "__main__":
    run_feature_extraction(
        input_folder="D:/project_fastq/results/fastp",
        output_csv="D:/project_fastq/features.csv"
    )