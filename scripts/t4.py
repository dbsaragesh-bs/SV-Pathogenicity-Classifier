import pandas as pd
b = pd.read_parquet("data/processed/batches.parquet")
print(b.groupby("batch").size().describe())
print(b.groupby("batch").size().values)
print(b.chrom.value_counts().head(8))
print(b.groupby("chrom").size().describe())
