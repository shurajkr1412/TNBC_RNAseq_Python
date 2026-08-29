import pandas as pd
from pathlib import Path


# ============================================================
# 1. Define input and output files
# ============================================================

metadata_file = Path("data/metadata/metadata.csv")
count_file = Path("results/clean_count_matrix.csv")
output_file = Path("results/metadata_for_deseq2.csv")


# ============================================================
# 2. Check that the required files exist
# ============================================================

if not metadata_file.exists():
    raise FileNotFoundError(
        f"Metadata file not found: {metadata_file}"
    )

if not count_file.exists():
    raise FileNotFoundError(
        f"Clean count matrix not found: {count_file}\n"
        "Run 01_load_counts.py first."
    )


# ============================================================
# 3. Read metadata
# ============================================================

metadata = pd.read_csv(metadata_file)

print("Metadata loaded successfully!")
print("Metadata shape:", metadata.shape)


# ============================================================
# 4. Display metadata columns
# ============================================================

print("\nMetadata columns:")
print(metadata.columns.tolist())


# ============================================================
# 5. Keep only the information needed for RNA-seq analysis
# ============================================================

metadata = metadata[["Run", "source_name"]].copy()


# ============================================================
# 6. Create the experimental condition
# ============================================================

def assign_condition(source_name):

    source_name = str(source_name).lower()

    if "tnbc subtype" in source_name:
        return "TNBC"

    elif "paracancerous" in source_name:
        return "Paracancerous"

    else:
        return "Other"


metadata["condition"] = metadata["source_name"].apply(assign_condition)


# ============================================================
# 7. Keep only the 8 samples present in our count matrix
# ============================================================

count_matrix = pd.read_csv(
    count_file,
    index_col=0
)

count_samples = count_matrix.columns.tolist()

metadata = metadata[
    metadata["Run"].isin(count_samples)
].copy()


# ============================================================
# 8. Check that all count samples are present in metadata
# ============================================================

missing_metadata = set(count_samples) - set(metadata["Run"])

if missing_metadata:
    raise ValueError(
        f"These count samples are missing from metadata: "
        f"{sorted(missing_metadata)}"
    )


# ============================================================
# 9. Check for unexpected conditions
# ============================================================

unexpected_conditions = set(metadata["condition"]) - {
    "TNBC",
    "Paracancerous"
}

if unexpected_conditions:
    raise ValueError(
        f"Unexpected conditions found: {unexpected_conditions}"
    )


# ============================================================
# 10. Reorder metadata to exactly match count matrix
# ============================================================

metadata = metadata.set_index("Run")

metadata = metadata.loc[count_samples]


# ============================================================
# 11. Keep only the columns needed by PyDESeq2
# ============================================================

metadata_for_deseq2 = metadata[["condition"]].copy()


# ============================================================
# 12. Display final metadata
# ============================================================

print("\nFinal metadata:")
print(metadata_for_deseq2)


# ============================================================
# 13. Display number of samples in each condition
# ============================================================

print("\nSamples per condition:")
print(metadata_for_deseq2["condition"].value_counts())


# ============================================================
# 14. Check that sample order matches
# ============================================================

print("\nCount matrix sample order:")
print(count_samples)

print("\nMetadata sample order:")
print(metadata_for_deseq2.index.tolist())


if count_samples == metadata_for_deseq2.index.tolist():
    print("\nSample order check: PASSED")
else:
    raise ValueError(
        "Sample order does not match between count matrix and metadata."
    )


# ============================================================
# 15. Save final metadata
# ============================================================

metadata_for_deseq2.to_csv(output_file)

print(f"\nMetadata saved to: {output_file}")

print("\n========== STEP 2 COMPLETED ==========")