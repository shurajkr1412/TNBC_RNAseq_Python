import pandas as pd
from pathlib import Path

# ============================================================
# 1. Define input file
# ============================================================

count_file = Path("data/counts/gene_counts.txt")


# ============================================================
# 2. Check that the file exists
# ============================================================

if not count_file.exists():
    raise FileNotFoundError(
        f"Count file not found: {count_file}\n"
        "Please check that gene_counts.txt is inside data/counts/"
    )


# ============================================================
# 3. Read featureCounts output
# ============================================================

counts = pd.read_csv(
    count_file,
    sep="\t",
    comment="#"
)

print("Count table loaded successfully!")
print("Original shape:", counts.shape)


# ============================================================
# 4. Display original columns
# ============================================================

print("\nOriginal columns:")
for column in counts.columns:
    print(column)


# ============================================================
# 5. Remove featureCounts annotation columns
# ============================================================

annotation_columns = [
    "Chr",
    "Start",
    "End",
    "Strand",
    "Length"
]

count_matrix = counts.drop(columns=annotation_columns)


# ============================================================
# 6. Use Geneid as the row identifier
# ============================================================

count_matrix = count_matrix.set_index("Geneid")
count_matrix = count_matrix.apply(pd.to_numeric, errors="coerce")


# ============================================================
# 7. Clean sample names
# ============================================================

count_matrix.columns = (
    count_matrix.columns
    .str.split("/")
    .str[-1]
    .str.replace(".sorted.bam", "", regex=False)
)


# ============================================================
# 8. Check the cleaned sample names
# ============================================================

print("\nClean sample names:")
for sample in count_matrix.columns:
    print(sample)


# ============================================================
# 9. Check the cleaned count matrix
# ============================================================

print("\nClean count matrix shape:", count_matrix.shape)

print("\nFirst 5 genes:")
print(count_matrix.head())


# ============================================================
# 10. Check for missing values
# ============================================================

print("\nMissing values:", count_matrix.isna().sum().sum())


# ============================================================
# 11. Check whether all counts are numeric
# ============================================================

non_numeric = count_matrix.apply(
    pd.to_numeric,
    errors="coerce"
).isna().sum().sum()

print("Non-numeric values:", non_numeric)


# ============================================================
# 12. Save cleaned count matrix
# ============================================================

output_file = Path("results/clean_count_matrix.csv")

count_matrix.to_csv(output_file)

print(f"\nClean count matrix saved to: {output_file}")

print("\n========== STEP 1 COMPLETED ==========")