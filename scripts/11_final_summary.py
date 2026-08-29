import pandas as pd
from pathlib import Path


# ==============================
# PROJECT PATHS
# ==============================

project_dir = Path(__file__).resolve().parent.parent
results_dir = project_dir / "results"

deg_file = results_dir / "final_annotated_DEGs.csv"
sig_deg_file = results_dir / "final_significant_DEGs_annotated.csv"


# ==============================
# LOAD RESULTS
# ==============================

deg = pd.read_csv(deg_file)
sig_deg = pd.read_csv(sig_deg_file)


# ==============================
# BASIC STATISTICS
# ==============================

total_genes = len(deg)
significant_degs = len(sig_deg)


# Detect fold-change column
if "log2FoldChange" in sig_deg.columns:
    fc_column = "log2FoldChange"
elif "log2FC" in sig_deg.columns:
    fc_column = "log2FC"
else:
    raise ValueError("log2FoldChange column not found.")


upregulated = sig_deg[sig_deg[fc_column] > 0]
downregulated = sig_deg[sig_deg[fc_column] < 0]


# ==============================
# TOP UPREGULATED GENES
# ==============================

top_up = (
    upregulated
    .sort_values(by=fc_column, ascending=False)
    .head(10)
)


# ==============================
# TOP DOWNREGULATED GENES
# ==============================

top_down = (
    downregulated
    .sort_values(by=fc_column, ascending=True)
    .head(10)
)


# ==============================
# SAVE SUMMARY
# ==============================

summary_file = results_dir / "final_project_summary.txt"


with open(summary_file, "w", encoding="utf-8") as f:

    f.write("============================================\n")
    f.write("TNBC RNA-seq DIFFERENTIAL EXPRESSION SUMMARY\n")
    f.write("============================================\n\n")

    f.write("PROJECT: TNBC vs Paracancerous Breast Tissue\n\n")

    f.write("1. OVERALL RESULTS\n")
    f.write("--------------------------------------------\n")
    f.write(f"Total genes analyzed: {total_genes}\n")
    f.write(f"Significant DEGs: {significant_degs}\n")
    f.write(f"Upregulated genes: {len(upregulated)}\n")
    f.write(f"Downregulated genes: {len(downregulated)}\n\n")

    f.write("2. TOP 10 UPREGULATED GENES\n")
    f.write("--------------------------------------------\n")

    for _, row in top_up.iterrows():

        gene_name = row.get(
            "Gene_Symbol",
            row.get("gene_symbol", row.get("symbol", row.get("Geneid", "Unknown")))
        )

        f.write(
            f"{gene_name}\t"
            f"log2FC={row[fc_column]:.4f}\n"
        )

    f.write("\n")

    f.write("3. TOP 10 DOWNREGULATED GENES\n")
    f.write("--------------------------------------------\n")

    for _, row in top_down.iterrows():

        gene_name = row.get(
            "Gene_Symbol",
            row.get("gene_symbol", row.get("symbol", row.get("Geneid", "Unknown")))
        )

        f.write(
            f"{gene_name}\t"
            f"log2FC={row[fc_column]:.4f}\n"
        )

    f.write("\n")

    f.write("4. ANALYSIS OUTPUTS\n")
    f.write("--------------------------------------------\n")
    f.write("Differential expression results: differential_expression_results.csv\n")
    f.write("Significant DEGs: significant_DEGs.csv\n")
    f.write("PCA coordinates: PCA_coordinates.csv\n")
    f.write("Annotated DEGs: final_annotated_DEGs.csv\n")
    f.write("Annotated significant DEGs: final_significant_DEGs_annotated.csv\n")
    f.write("PCA plot: PCA_plot.png\n")
    f.write("Volcano plot: Volcano_plot.png\n")
    f.write("MA plot: MA_plot.png\n")
    f.write("Heatmap: Heatmap_Top50_DEGs.png\n")
    f.write("GO enrichment results: GO_enrichment/\n")
    f.write("KEGG enrichment results: KEGG_enrichment/\n")


# ==============================
# TERMINAL OUTPUT
# ==============================

print("\n============================================")
print("STEP 11 COMPLETED")
print("============================================")

print(f"Total genes analyzed: {total_genes}")
print(f"Significant DEGs: {significant_degs}")
print(f"Upregulated genes: {len(upregulated)}")
print(f"Downregulated genes: {len(downregulated)}")

print(f"\nFinal summary saved to:")
print(f"{summary_file}")