import pandas as pd
from pathlib import Path


# ============================================================
# 1. Define input files
# ============================================================

deg_file = Path(
    "results/differential_expression_results.csv"
)

mapping_file = Path(
    "results/GO_enrichment/DEG_gene_symbols.csv"
)


# ============================================================
# 2. Define output file
# ============================================================

output_file = Path(
    "results/final_annotated_DEGs.csv"
)


# ============================================================
# 3. Check that files exist
# ============================================================

if not deg_file.exists():
    raise FileNotFoundError(
        f"DEG results not found: {deg_file}"
    )

if not mapping_file.exists():
    raise FileNotFoundError(
        f"Gene mapping file not found: {mapping_file}"
    )


# ============================================================
# 4. Read DEG results
# ============================================================

deg = pd.read_csv(
    deg_file,
    index_col=0
)


# ============================================================
# 5. Read Ensembl-to-symbol mapping
# ============================================================

mapping = pd.read_csv(
    mapping_file
)


# ============================================================
# 6. Remove Ensembl version numbers
# ============================================================

deg["ensembl_id"] = (
    deg.index
    .astype(str)
    .str.split(".")
    .str[0]
)


# ============================================================
# 7. Prepare mapping table
# ============================================================

mapping["ensembl_id"] = (
    mapping["ensembl_id"]
    .astype(str)
    .str.split(".")
    .str[0]
)

mapping = mapping[
    ["ensembl_id", "symbol"]
].drop_duplicates(
    "ensembl_id"
)


# ============================================================
# 8. Add gene symbols
# ============================================================

deg = deg.merge(
    mapping,
    on="ensembl_id",
    how="left"
)


# ============================================================
# 9. Define DEG category
# ============================================================

deg["category"] = "Not significant"

deg.loc[
    (deg["padj"] < 0.05) &
    (deg["log2FoldChange"] >= 1),
    "category"
] = "Upregulated"

deg.loc[
    (deg["padj"] < 0.05) &
    (deg["log2FoldChange"] <= -1),
    "category"
] = "Downregulated"


# ============================================================
# 10. Move important columns to the front
# ============================================================

columns = [
    "ensembl_id",
    "symbol",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
    "category"
]

deg = deg[
    [column for column in columns if column in deg.columns]
]


# ============================================================
# 11. Sort by adjusted p-value
# ============================================================

deg = deg.sort_values(
    "padj",
    na_position="last"
)


# ============================================================
# 12. Save final annotated table
# ============================================================

deg.to_csv(
    output_file,
    index=False
)


# ============================================================
# 13. Create significant DEG table
# ============================================================

significant = deg[
    (deg["padj"] < 0.05) &
    (abs(deg["log2FoldChange"]) >= 1)
]


significant_file = Path(
    "results/final_significant_DEGs_annotated.csv"
)

significant.to_csv(
    significant_file,
    index=False
)


# ============================================================
# 14. Count upregulated and downregulated genes
# ============================================================

upregulated = (
    significant["category"] == "Upregulated"
).sum()

downregulated = (
    significant["category"] == "Downregulated"
).sum()


# ============================================================
# 15. Display summary
# ============================================================

print("\nFinal DEG annotation completed!")

print(
    "Total genes analyzed:",
    len(deg)
)

print(
    "Significant DEGs:",
    len(significant)
)

print(
    "Upregulated genes:",
    upregulated
)

print(
    "Downregulated genes:",
    downregulated
)

print(
    "\nAnnotated DEG table saved to:",
    output_file
)

print(
    "Annotated significant DEG table saved to:",
    significant_file
)

print(
    "\n========== STEP 10 COMPLETED =========="
)