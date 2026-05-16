import pandas as pd

def build_dataset(features_path, kmer_path, metadata_path, output_path):

    features_df = pd.read_csv(features_path)
    kmer_df = pd.read_csv(kmer_path)
    meta_df = pd.read_csv(metadata_path)

    # نخلي sample موحد
    meta_df = meta_df[["Run", "disease"]]
    meta_df = meta_df.rename(columns={"Run": "sample"})

    # تنظيف اسم sample في features لو فيه _fastp
    features_df["sample"] = features_df["sample"].str.replace("_fastp", "")

    # merge 1: features + kmer
    merged = pd.merge(features_df, kmer_df, on="sample")

    # merge 2: add labels
    merged = pd.merge(merged, meta_df, on="sample")

    merged = merged.fillna(0)

    merged.to_csv(output_path, index=False)

    print("ML dataset created successfully!")
    print(f"Saved at: {output_path}")


if __name__ == "__main__":
    build_dataset(
        features_path="D:/project_fastq/features.csv",
        kmer_path="D:/project_fastq/kmer_features.csv",
        metadata_path="D:/project_fastq/SraRunTable.csv",
        output_path="D:/project_fastq/ml_dataset.csv"
    )