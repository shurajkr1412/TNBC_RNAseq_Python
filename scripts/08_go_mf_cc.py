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
    "results/GO_enrichment"
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
    "Gene symbols available for enrichment:",
    len(gene_symbols)
)


# ============================================================
# 4. Function for GO enrichment
# ============================================================

def run_go_enrichment(
    gene_list,
    gene_set,
    output_folder,
    result_file,
    plot_file,
    title
):

    print(f"\nRunning {title} enrichment...")

    enrichment = gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_set,
        organism="human",
        outdir=str(output_folder),
        cutoff=0.05
    )

    results = enrichment.results

    results.to_csv(
        result_file,
        index=False
    )

    significant = results[
        results["Adjusted P-value"] < 0.05
    ].sort_values(
        "Adjusted P-value"
    ).head(10)

    if len(significant) > 0:

        plt.figure(
            figsize=(10, 7)
        )

        plt.barh(
            significant["Term"].iloc[::-1],
            -np.log10(
                significant[
                    "Adjusted P-value"
                ].iloc[::-1]
            )
        )

        plt.xlabel(
            "-log10 Adjusted P-value"
        )

        plt.ylabel(
            title
        )

        plt.title(
            f"Top {title} Enriched Terms"
        )

        plt.tight_layout()

        plt.savefig(
            plot_file,
            dpi=300,
            bbox_inches="tight"
        )

        plt.show()

        print(
            f"Plot saved to: {plot_file}"
        )

    else:

        print(
            "No terms passed adjusted "
            "p-value < 0.05."
        )

    print(
        f"Results saved to: {result_file}"
    )

    return results


# ============================================================
# 5. GO Molecular Function
# ============================================================

mf_results = run_go_enrichment(
    gene_symbols,
    "GO_Molecular_Function_2023",
    output_dir / "GO_MF",
    output_dir / "GO_MF_results.csv",
    output_dir / "GO_MF_top10.png",
    "GO Molecular Function"
)


# ============================================================
# 6. GO Cellular Component
# ============================================================

cc_results = run_go_enrichment(
    gene_symbols,
    "GO_Cellular_Component_2023",
    output_dir / "GO_CC",
    output_dir / "GO_CC_results.csv",
    output_dir / "GO_CC_top10.png",
    "GO Cellular Component"
)


# ============================================================
# 7. Final message
# ============================================================

print(
    "\nGO Molecular Function and "
    "Cellular Component analysis completed!"
)

print(
    "\n========== STEP 8 COMPLETED =========="
)