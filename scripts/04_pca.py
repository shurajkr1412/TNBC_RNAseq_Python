import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA


# ============================================================
# 1. Define input files
# ============================================================

count_file = Path("results/clean_count_matrix.csv")
metadata_file = Path("results/metadata_for_deseq2.csv")


# ============================================================
# 2. Check that files exist
# ============================================================

if not count_file.exists():
    raise FileNotFoundError(
        f"Count matrix not found: {count_file}"
    )

if not metadata_file.exists():
    raise FileNotFoundError(
        f"Metadata file not found: {metadata_file}"
    )


# ============================================================
# 3. Read count matrix and metadata
# ============================================================

counts = pd.read_csv(
    count_file,
    index_col=0
)

metadata = pd.read_csv(
    metadata_file,
    index_col=0
)


# ============================================================
# 4. Transpose counts
#    Rows = samples
#    Columns = genes
# ============================================================

counts = counts.T


# ============================================================
# 5. Match metadata to count matrix
# ============================================================

metadata = metadata.loc[counts.index]


# ============================================================
# 6. Convert counts to numeric
# ============================================================

counts = counts.astype(float)


# ============================================================
# 7. Normalize counts by library size
# ============================================================

library_sizes = counts.sum(axis=1)

normalized_counts = counts.div(
    library_sizes,
    axis=0
) * 1_000_000


# ============================================================
# 8. Log-transform normalized counts
# ============================================================

log_counts = np.log2(
    normalized_counts + 1
)


# ============================================================
# 9. Select the most variable genes
# ============================================================

gene_variances = log_counts.var(axis=0)

top_genes = gene_variances.nlargest(500).index

pca_data = log_counts[top_genes]


# ============================================================
# 10. Perform PCA
# ============================================================

pca = PCA(
    n_components=2
)

principal_components = pca.fit_transform(pca_data)


# ============================================================
# 11. Create PCA DataFrame
# ============================================================

pca_results = pd.DataFrame(
    principal_components,
    columns=["PC1", "PC2"],
    index=pca_data.index
)

pca_results["condition"] = metadata["condition"]


# ============================================================
# 12. Calculate explained variance
# ============================================================

pc1_variance = (
    pca.explained_variance_ratio_[0] * 100
)

pc2_variance = (
    pca.explained_variance_ratio_[1] * 100
)


# ============================================================
# 13. Create PCA plot
# ============================================================

plt.figure(figsize=(8, 6))

for condition in pca_results["condition"].unique():

    subset = pca_results[
        pca_results["condition"] == condition
    ]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=condition,
        s=100
    )

    for sample in subset.index:

        plt.annotate(
            sample,
            (
                subset.loc[sample, "PC1"],
                subset.loc[sample, "PC2"]
            ),
            xytext=(5, 5),
            textcoords="offset points"
        )


plt.xlabel(
    f"PC1 ({pc1_variance:.2f}% variance)"
)

plt.ylabel(
    f"PC2 ({pc2_variance:.2f}% variance)"
)

plt.title(
    "PCA Plot - TNBC vs Paracancerous"
)

plt.legend()

plt.tight_layout()


# ============================================================
# 14. Save PCA plot
# ============================================================

output_plot = Path(
    "results/PCA_plot.png"
)

plt.savefig(
    output_plot,
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# 15. Save PCA coordinates
# ============================================================

output_table = Path(
    "results/PCA_coordinates.csv"
)

pca_results.to_csv(
    output_table
)


# ============================================================
# 16. Final message
# ============================================================

print("\nPCA analysis completed!")

print(
    f"PC1 explains {pc1_variance:.2f}% of the variance."
)

print(
    f"PC2 explains {pc2_variance:.2f}% of the variance."
)

print(
    f"\nPCA plot saved to: {output_plot}"
)

print(
    f"PCA coordinates saved to: {output_table}"
)

print("\n========== STEP 4 COMPLETED ==========")