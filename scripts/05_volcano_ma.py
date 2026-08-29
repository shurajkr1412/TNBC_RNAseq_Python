import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# 1. Define input and output files
# ============================================================

results_file = Path(
    "results/differential_expression_results.csv"
)

volcano_file = Path(
    "results/Volcano_plot.png"
)

ma_file = Path(
    "results/MA_plot.png"
)


# ============================================================
# 2. Check that the results file exists
# ============================================================

if not results_file.exists():
    raise FileNotFoundError(
        f"DE results not found: {results_file}\n"
        "Run 03_differential_expression.py first."
    )


# ============================================================
# 3. Read differential-expression results
# ============================================================

results = pd.read_csv(
    results_file,
    index_col=0
)


# ============================================================
# 4. Define significance criteria
# ============================================================

padj_cutoff = 0.05
log2fc_cutoff = 1


# ============================================================
# 5. Create significance categories
# ============================================================

results["significance"] = "Not significant"

results.loc[
    (results["padj"] < padj_cutoff) &
    (results["log2FoldChange"] >= log2fc_cutoff),
    "significance"
] = "Upregulated"

results.loc[
    (results["padj"] < padj_cutoff) &
    (results["log2FoldChange"] <= -log2fc_cutoff),
    "significance"
] = "Downregulated"


# ============================================================
# 6. Calculate -log10 adjusted p-value
# ============================================================

results["minus_log10_padj"] = -np.log10(
    results["padj"].clip(lower=1e-300)
)


# ============================================================
# 7. Count categories
# ============================================================

upregulated = (
    results["significance"] == "Upregulated"
).sum()

downregulated = (
    results["significance"] == "Downregulated"
).sum()

not_significant = (
    results["significance"] == "Not significant"
).sum()


# ============================================================
# 8. Create Volcano Plot
# ============================================================

plt.figure(figsize=(10, 7))

for category in [
    "Not significant",
    "Upregulated",
    "Downregulated"
]:

    subset = results[
        results["significance"] == category
    ]

    plt.scatter(
        subset["log2FoldChange"],
        subset["minus_log10_padj"],
        s=12,
        alpha=0.6,
        label=category
    )


plt.axvline(
    log2fc_cutoff,
    linestyle="--"
)

plt.axvline(
    -log2fc_cutoff,
    linestyle="--"
)

plt.axhline(
    -np.log10(padj_cutoff),
    linestyle="--"
)

plt.xlabel("log2 Fold Change")

plt.ylabel("-log10 Adjusted P-value")

plt.title(
    "Volcano Plot - TNBC vs Paracancerous"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    volcano_file,
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# 9. Create MA Plot
# ============================================================

plt.figure(figsize=(10, 7))

for category in [
    "Not significant",
    "Upregulated",
    "Downregulated"
]:

    subset = results[
        results["significance"] == category
    ]

    plt.scatter(
        subset["baseMean"],
        subset["log2FoldChange"],
        s=12,
        alpha=0.6,
        label=category
    )


plt.axhline(
    0,
    linestyle="--"
)

plt.axhline(
    log2fc_cutoff,
    linestyle="--"
)

plt.axhline(
    -log2fc_cutoff,
    linestyle="--"
)

plt.xscale("log")

plt.xlabel("Mean Expression (baseMean)")

plt.ylabel("log2 Fold Change")

plt.title(
    "MA Plot - TNBC vs Paracancerous"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    ma_file,
    dpi=300,
    bbox_inches="tight"
)

plt.show()


# ============================================================
# 10. Save categorized results
# ============================================================

categorized_file = Path(
    "results/differential_expression_categorized.csv"
)

results.to_csv(
    categorized_file
)


# ============================================================
# 11. Final summary
# ============================================================

print("\nVolcano and MA plots completed!")

print("Upregulated genes:", upregulated)
print("Downregulated genes:", downregulated)
print("Not significant:", not_significant)

print(
    f"\nVolcano plot saved to: {volcano_file}"
)

print(
    f"MA plot saved to: {ma_file}"
)

print(
    f"Categorized results saved to: "
    f"{categorized_file}"
)

print("\n========== STEP 5 COMPLETED ==========")