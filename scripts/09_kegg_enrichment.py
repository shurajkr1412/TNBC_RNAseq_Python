import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np


# ============================================================
# 1. Define input and output locations
# ============================================================

mapped_file = Path(
    "results/GO_enrichment/DEG_gene_symbols.csv"
)

output_dir = Path(
    "results/KEGG_enrichment"
)

output_dir.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. Check input file
# ============================================================

if not mapped_file.exists():
    raise FileNotFoundError(
        f"Mapped gene file not found: {mapped_file}\n"
        "Run 07_go_enrichment.py first."
    )


# ============================================================
# 3. Read mapped gene symbols
# ============================================================

mapping = pd.read_csv(
    mapped_file
)

gene_symbols = (
    mapping["symbol"]
    .dropna()
    .astype(str)
    .str.upper()
    .drop_duplicates()
    .tolist()
)

print(
    "Genes available for KEGG enrichment:",
    len(gene_symbols)
)


# ============================================================
# 4. Run KEGG enrichment
# ============================================================

print("\nRunning KEGG pathway enrichment...")

kegg = gp.enrichr(
    gene_list=gene_symbols,
    gene_sets="KEGG_2021_Human",
    organism="human",
    outdir=str(output_dir),
    cutoff=0.05
)


# ============================================================
# 5. Get KEGG results
# ============================================================

kegg_results = kegg.results


# ============================================================
# 6. Save complete KEGG results
# ============================================================

kegg_file = (
    output_dir / "KEGG_results.csv"
)

kegg_results.to_csv(
    kegg_file,
    index=False
)


# ============================================================
# 7. Select significant pathways
# ============================================================

significant_kegg = kegg_results[
    kegg_results["Adjusted P-value"] < 0.05
].sort_values(
    "Adjusted P-value"
).head(10)


# ============================================================
# 8. Create top 10 KEGG plot
# ============================================================

if len(significant_kegg) > 0:

    plt.figure(
        figsize=(10, 7)
    )

    plt.barh(
        significant_kegg["Term"].iloc[::-1],
        -np.log10(
            significant_kegg[
                "Adjusted P-value"
            ].iloc[::-1]
        )
    )

    plt.xlabel(
        "-log10 Adjusted P-value"
    )

    plt.ylabel(
        "KEGG Pathway"
    )

    plt.title(
        "Top 10 KEGG Enriched Pathways"
    )

    plt.tight_layout()

    plot_file = (
        output_dir /
        "KEGG_top10.png"
    )

    plt.savefig(
        plot_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

    print(
        f"\nKEGG plot saved to: {plot_file}"
    )

else:

    print(
        "\nNo KEGG pathways passed "
        "adjusted p-value < 0.05."
    )


# ============================================================
# 9. Final message
# ============================================================

print(
    f"\nComplete KEGG results saved to: "
    f"{kegg_file}"
)

print(
    "\n========== STEP 9 COMPLETED =========="
)