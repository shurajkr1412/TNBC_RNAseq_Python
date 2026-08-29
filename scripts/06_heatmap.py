import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path


# ============================================================
# 1. Define input files
# ============================================================

count_file = Path("results/clean_count_matrix.csv")
metadata_file = Path("results/metadata_for_deseq2.csv")
deg_file = Path("results/differential_expression_results.csv")

output_file = Path("results/Heatmap_Top50_DEGs.png")


# ============================================================
# 2. Check that required files exist
# ============================================================

for file in [count_file, metadata_file, deg_file]:
    if not file.exists():
        raise FileNotFoundError(
            f"Required file not found: {file}"
        )


# ============================================================
# 3. Read the data
# ============================================================

counts = pd.read_csv(
    count_file,
    index_col=0
)

metadata = pd.read_csv(
    metadata_file,
    index_col=0
)

deg_results = pd.read_csv(
    deg_file,
    index_col=0
)


# ============================================================
# 4. Select the 50 most significant DEGs
# ============================================================

top_degs = deg_results[
    deg_results["padj"].notna()
].sort_values(
    "padj"
).head(50)


# ============================================================
# 5. Extract expression values for these genes
# ============================================================

heatmap_data = counts.loc[
    top_degs.index
]


# ============================================================
# 6. Normalize each gene for visualization
# ============================================================

log_counts = np.log2(
    heatmap_data + 1
)

gene_means = log_counts.mean(axis=1)

gene_stds = log_counts.std(axis=1)

z_scores = log_counts.sub(
    gene_means,
    axis=0
).div(
    gene_stds,
    axis=0
)


# ============================================================
# 7. Order samples by condition
# ============================================================

sample_order = metadata.sort_values(
    "condition"
).index

z_scores = z_scores[sample_order]


# ============================================================
# 8. Create the heatmap
# ============================================================

plt.figure(
    figsize=(10, 14)
)

sns.heatmap(
    z_scores,
    cmap="vlag",
    center=0,
    xticklabels=True,
    yticklabels=True
)

plt.title(
    "Top 50 Differentially Expressed Genes"
)

plt.xlabel("Samples")

plt.ylabel("Genes")

plt.xticks(
    rotation=45,
    ha="right"
)

plt.tight_layout()


# ============================================================
# 9. Save heatmap
# ============================================================

plt.savefig(
    output_file,
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# 10. Save top DEG list
# ============================================================

top_deg_file = Path(
    "results/Top50_DEGs.csv"
)

top_degs.to_csv(
    top_deg_file
)


# ============================================================
# 11. Final message
# ============================================================

print("\nHeatmap completed!")

print(
    "Number of genes plotted:",
    len(top_degs)
)

print(
    f"Heatmap saved to: {output_file}"
)

print(
    f"Top DEG list saved to: {top_deg_file}"
)

print("\n========== STEP 6 COMPLETED ==========")