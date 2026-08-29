import pandas as pd
from pathlib import Path
import mygene
import gseapy as gp
import matplotlib.pyplot as plt


# ============================================================
# 1. Define files
# ============================================================

deg_file = Path(
    "results/significant_DEGs.csv"
)

output_dir = Path(
    "results/GO_enrichment"
)

output_dir.mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# 2. Check input file
# ============================================================

if not deg_file.exists():
    raise FileNotFoundError(
        f"Significant DEG file not found: {deg_file}"
    )


# ============================================================
# 3. Read significant DEGs
# ============================================================

deg = pd.read_csv(
    deg_file,
    index_col=0
)

print("Significant DEG table loaded.")
print("Number of significant DEGs:", len(deg))


# ============================================================
# 4. Extract Ensembl gene IDs
# ============================================================

ensembl_ids = (
    deg.index
    .astype(str)
    .str.split(".")
    .str[0]
    .tolist()
)


# ============================================================
# 5. Convert Ensembl IDs to gene symbols
# ============================================================

print("\nConverting Ensembl IDs to gene symbols...")

mg = mygene.MyGeneInfo()

mapping = mg.querymany(
    ensembl_ids,
    scopes="ensembl.gene",
    fields="symbol",
    species="human",
    as_dataframe=True,
    returnall=False
)


# ============================================================
# 6. Clean the mapping table
# ============================================================

mapping = mapping.reset_index()

mapping = mapping.rename(
    columns={
        "query": "ensembl_id"
    }
)

mapping = mapping[
    ["ensembl_id", "symbol"]
]

mapping = mapping.dropna(
    subset=["symbol"]
)

mapping = mapping.drop_duplicates(
    subset=["ensembl_id"]
)


# ============================================================
# 7. Create gene-symbol list
# ============================================================

gene_symbols = (
    mapping["symbol"]
    .astype(str)
    .str.upper()
    .tolist()
)

print(
    "Gene symbols successfully mapped:",
    len(gene_symbols)
)


# ============================================================
# 8. Save mapped gene list
# ============================================================

mapped_file = output_dir / "DEG_gene_symbols.csv"

mapping.to_csv(
    mapped_file,
    index=False
)


# ============================================================
# 9. Run GO Biological Process enrichment
# ============================================================

print("\nRunning GO Biological Process enrichment...")

go_bp = gp.enrichr(
    gene_list=gene_symbols,
    gene_sets="GO_Biological_Process_2023",
    organism="human",
    outdir=str(output_dir / "GO_BP"),
    cutoff=0.05
)


# ============================================================
# 10. Get GO results
# ============================================================

go_results = go_bp.results


# ============================================================
# 11. Save complete GO results
# ============================================================

go_results_file = (
    output_dir / "GO_BP_results.csv"
)

go_results.to_csv(
    go_results_file,
    index=False
)


# ============================================================
# 12. Select top 10 enriched biological processes
# ============================================================

top_go = go_results[
    go_results["Adjusted P-value"] < 0.05
].sort_values(
    "Adjusted P-value"
).head(10)


# ============================================================
# 13. Create bar plot
# ============================================================

if len(top_go) > 0:

    plt.figure(
        figsize=(10, 7)
    )

    plt.barh(
        top_go["Term"].iloc[::-1],
        -__import__("numpy").log10(
            top_go["Adjusted P-value"].iloc[::-1]
        )
    )

    plt.xlabel(
        "-log10 Adjusted P-value"
    )

    plt.ylabel(
        "GO Biological Process"
    )

    plt.title(
        "Top GO Biological Processes"
    )

    plt.tight_layout()

    plot_file = (
        output_dir /
        "GO_BP_top10.png"
    )

    plt.savefig(
        plot_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.show()

    print(
        f"\nGO plot saved to: {plot_file}"
    )

else:

    print(
        "\nNo GO terms passed adjusted "
        "p-value < 0.05."
    )


# ============================================================
# 14. Final message
# ============================================================

print(
    f"\nComplete GO results saved to: "
    f"{go_results_file}"
)

print(
    f"Mapped gene list saved to: "
    f"{mapped_file}"
)

print("\n========== STEP 7 COMPLETED ==========")