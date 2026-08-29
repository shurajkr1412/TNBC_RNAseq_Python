import pandas as pd
from pathlib import Path

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


# ============================================================
# 1. Define input files
# ============================================================

count_file = Path("results/clean_count_matrix.csv")
metadata_file = Path("results/metadata_for_deseq2.csv")


# ============================================================
# 2. Check that the files exist
# ============================================================

if not count_file.exists():
    raise FileNotFoundError(
        f"Count matrix not found: {count_file}"
    )

if not metadata_file.exists():
    raise FileNotFoundError(
        f"Metadata not found: {metadata_file}\n"
        "Run 02_prepare_metadata.py first."
    )


# ============================================================
# 3. Read the count matrix
# ============================================================

counts = pd.read_csv(
    count_file,
    index_col=0
)


# ============================================================
# 4. Read the metadata
# ============================================================

metadata = pd.read_csv(
    metadata_file,
    index_col=0
)


# ============================================================
# 5. Display the original dimensions
# ============================================================

print("Original count matrix shape:", counts.shape)
print("Metadata shape:", metadata.shape)


# ============================================================
# 6. Transpose the count matrix
# ============================================================

counts = counts.T


# ============================================================
# 7. Make sure sample order matches metadata
# ============================================================

metadata = metadata.loc[counts.index]


# ============================================================
# 8. Convert counts to integers
# ============================================================

counts = counts.round().astype(int)


# ============================================================
# 9. Remove genes with zero counts in every sample
# ============================================================

counts = counts.loc[
    :, counts.sum(axis=0) > 0
]


# ============================================================
# 10. Check the final dimensions
# ============================================================

print("\nAfter transposing:")
print("Count matrix shape:", counts.shape)
print("Metadata shape:", metadata.shape)


# ============================================================
# 11. Check that samples match
# ============================================================

if list(counts.index) != list(metadata.index):
    raise ValueError(
        "Sample names or sample order do not match "
        "between counts and metadata."
    )

print("\nSample matching: PASSED")


# ============================================================
# 12. Check the experimental conditions
# ============================================================

print("\nExperimental conditions:")
print(metadata["condition"].value_counts())


# ============================================================
# 13. Create the PyDESeq2 dataset
# ============================================================

dds = DeseqDataSet(
    counts=counts,
    metadata=metadata,
    design="~condition",
    refit_cooks=True,
    n_cpus=1
)


# ============================================================
# 14. Run the DESeq2 analysis
# ============================================================

print("\nRunning PyDESeq2...")
dds.deseq2()


# ============================================================
# 15. Perform statistical testing
# ============================================================

print("\nRunning statistical tests...")

stat_res = DeseqStats(
    dds,
    contrast=["condition", "TNBC", "Paracancerous"],
    n_cpus=1
)

stat_res.summary()


# ============================================================
# 16. Get the results
# ============================================================

results = stat_res.results_df


# ============================================================
# 17. Sort by adjusted p-value
# ============================================================

results = results.sort_values(
    by="padj",
    na_position="last"
)


# ============================================================
# 18. Save complete DE results
# ============================================================

output_file = Path(
    "results/differential_expression_results.csv"
)

results.to_csv(output_file)


# ============================================================
# 19. Select significant DEGs
# ============================================================

significant = results[
    (results["padj"] < 0.05) &
    (abs(results["log2FoldChange"]) >= 1)
]


# ============================================================
# 20. Save significant DEGs
# ============================================================

significant_file = Path(
    "results/significant_DEGs.csv"
)

significant.to_csv(significant_file)


# ============================================================
# 21. Display final results
# ============================================================

print("\nDifferential expression completed!")

print("Total genes tested:", len(results))
print("Significant DEGs:", len(significant))

print(
    "\nComplete results saved to:",
    output_file
)

print(
    "Significant DEGs saved to:",
    significant_file
)

print("\n========== STEP 3 COMPLETED ==========")