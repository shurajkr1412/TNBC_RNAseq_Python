# TNBC RNA-seq Differential Expression Analysis

## Project Overview

This project presents a computational RNA-seq analysis of Triple-Negative Breast Cancer (TNBC) compared with paracancerous breast tissue.

The analysis was performed using a Python-based RNA-seq workflow, including PyDESeq2 for differential expression analysis, PCA for sample-level variation, and GO and KEGG enrichment analyses for biological interpretation.

---

## Research Objective

The main objective of this project is to identify genes that are differentially expressed between TNBC and paracancerous breast tissue samples and to characterize their biological functions and associated pathways.

### Specific Objectives

- Process and organize RNA-seq gene count data.
- Prepare sample metadata for differential expression analysis.
- Identify differentially expressed genes.
- Determine significantly upregulated and downregulated genes.
- Perform principal component analysis (PCA).
- Generate MA, volcano, and heatmap visualizations.
- Perform Gene Ontology (GO) enrichment analysis.
- Perform KEGG pathway enrichment analysis.
- Annotate significant DEGs with gene symbols.
- Generate tables and visualizations for downstream interpretation.

---

## Dataset

The dataset used for this project is:

**GSE267442**

The analysis contains:

- 4 TNBC samples
- 4 Paracancerous breast tissue samples
- Human RNA-seq data

Gene-level read counts were used as the input for downstream analysis.

---
## Analysis Workflow

The overall workflow followed in this project was:

```text
RNA-seq Gene Count Data
          |
          v
Count Matrix Preparation
          |
          v
Sample Metadata Preparation
          |
          +----------------------+
          |                      |
          v                      v
Differential Expression     PCA Analysis
Analysis                         |
          |                      v
          v                 PCA Plot
Significant DEGs
          |
          +----------------------+----------------------+
          |                      |                      |
          v                      v                      v
     Volcano Plot             MA Plot              Heatmap
          |
          v
Gene Annotation
          |
          v
Annotated DEGs
          |
          +----------------------+
          |                      |
          v                      v
   GO Enrichment          KEGG Enrichment
          |                      |
          v                      v
Functional Interpretation   Pathway Interpretation




## Tools and Technologies

### Programming Language
- Python

### Python Libraries
- pandas
- NumPy
- SciPy
- matplotlib
- seaborn
- scikit-learn
- Biopython
- gseapy
- PyDESeq2

### Bioinformatics Analysis
- RNA-seq gene count analysis
- Differential expression analysis
- Principal Component Analysis (PCA)
- Gene Ontology (GO) enrichment
- KEGG pathway enrichment
- Gene annotation

### Visualization
- PCA plot
- Volcano plot
- MA plot
- Heatmap
- GO enrichment plots
- KEGG pathway enrichment plot

---

## Project Structure

```text
TNBC_RNAseq_Python/
│
├── data/
│   ├── bam/
│   ├── counts/
│   ├── fastq/
│   └── metadata/
│
├── results/
│   ├── GO_enrichment/
│   ├── KEGG_enrichment/
│   ├── clean_count_matrix.csv
│   ├── differential_expression_results.csv
│   ├── significant_DEGs.csv
|   ├── differential_expression_categorized.csv
│   ├── Top50_DEGs.csv
│   ├── PCA_coordinates.csv
│   ├── PCA_plot.png
│   ├── MA_plot.png
│   ├── Volcano_plot.png
│   ├── Heatmap_Top50_DEGs.png
│   ├── final_annotated_DEGs.csv
│   ├── final_significant_DEGs_annotated.csv
│   └── final_project_summary.txt
│
├── scripts/
│   ├── 01_load_counts.py
│   ├── 02_prepare_metadata.py
│   ├── 03_differential_expression.py
│   ├── 04_pca.py
│   ├── 05_volcano_ma.py
│   ├── 06_heatmap.py
│   ├── 07_go_enrichment.py
│   ├── 08_go_mf_cc.py
│   ├── 09_kegg_enrichment.py
│   ├── 10_annotate_degs.py
│   └── 11_final_summary.py
│
└── README.md










---

## Results

### Differential Expression Analysis

A total of **47,726 genes** were analyzed for differential expression between TNBC and paracancerous breast tissue samples.

Using an adjusted p-value (padj) < 0.05 and an absolute log2 fold change (|log2FC|) ≥ 1, **3,559 genes** were identified as significant differentially expressed genes (DEGs).

Among the significant DEGs:

- **1,878 genes were upregulated**
- **1,681 genes were downregulated**

The complete differential expression results are provided in:

`results/differential_expression_results.csv`

The significant DEGs are provided in:

`results/significant_DEGs.csv`

---

### Principal Component Analysis

Principal Component Analysis (PCA) was performed to evaluate sample-level variation and visualize the separation between TNBC and paracancerous samples.

The first principal component (PC1) explained **44.10%** of the total variance, while the second principal component (PC2) explained **26.07%**.

The PCA plot demonstrated separation between the TNBC and paracancerous sample groups.

**Output:**

`results/PCA_plot.png`

---

### Differential Expression Visualization

Several visualization methods were generated to evaluate the differential expression results:
The volcano plot and MA plot were generated from the differential expression results, while the heatmap shows expression patterns for the top 50 DEGs.

- **Volcano plot** – visualizes statistical significance and magnitude of gene expression changes.
- **MA plot** – displays the relationship between average gene expression and log2 fold change.
- **Heatmap** – visualizes expression patterns of the top differentially expressed genes.

The corresponding files are available in the `results/` directory.

---

### Gene Ontology Enrichment Analysis

Gene Ontology (GO) enrichment analysis was performed to investigate the biological functions associated with the identified significant DEGs.

Enrichment analysis was performed for:

- Biological Process (BP)
- Molecular Function (MF)
- Cellular Component (CC)

The GO enrichment results and plots are available under:

`results/GO_enrichment/`

---

### KEGG Pathway Enrichment Analysis

KEGG pathway enrichment analysis was performed to identify biological pathways associated with the significant DEGs.

The pathway enrichment results and visualization are available under:

`results/KEGG_enrichment/`

---

### Gene Annotation

Significant DEGs were annotated with gene symbols to facilitate biological interpretation.

The final annotated datasets are:

- `results/final_annotated_DEGs.csv`
- `results/final_significant_DEGs_annotated.csv`

---

## Conclusion

This project established a Python-based RNA-seq differential expression workflow for comparing Triple-Negative Breast Cancer (TNBC) with paracancerous breast tissue.

The analysis identified **3,559 significant differentially expressed genes (DEGs)**, including **1,878 upregulated genes** and **1,681 downregulated genes**.

Principal Component Analysis (PCA) demonstrated separation between the TNBC and paracancerous sample groups, indicating differences in their overall gene expression profiles.

Volcano plots, MA plots, and heatmap visualization were generated to examine the expression patterns of the identified DEGs. Gene Ontology (GO) enrichment analysis was performed for Biological Process, Molecular Function, and Cellular Component, while KEGG pathway enrichment was used to investigate pathways associated with the significant DEGs.

The final annotated DEG tables provide gene-level information for further biological interpretation.

Overall, the project demonstrates a reproducible Python-based workflow for RNA-seq differential expression analysis, visualization, functional enrichment, pathway analysis, and DEG annotation.

---

## Future Scope

The analysis can be further extended by:

- Performing additional pathway and network analysis.
- Investigating key hub genes among the significant DEGs.
- Comparing identified genes with known TNBC-associated biomarkers.
- Performing protein-protein interaction (PPI) network analysis.
- Integrating additional transcriptomic or multi-omics datasets.
- Validating candidate genes using independent datasets.



---
## How to Run the Project

```bash
git clone <your-github-repository-url>
cd TNBC_RNAseq_Python

python -m venv .venv

.venv\Scripts\activate

pip install pandas numpy scipy matplotlib seaborn scikit-learn biopython gseapy pydeseq2 anndata statsmodels

python scripts/01_load_counts.py
python scripts/02_prepare_metadata.py
python scripts/03_differential_expression.py
python scripts/04_pca.py
python scripts/05_volcano_ma.py
python scripts/06_heatmap.py
python scripts/07_go_enrichment.py
python scripts/08_go_mf_cc.py
python scripts/09_kegg_enrichment.py
python scripts/10_annotate_degs.py
python scripts/11_final_summary.py

The generated tables and visualizations will be stored in the `results/` directory.
---

## Reproducibility

All analysis scripts used to generate the results are included in the `scripts/` directory.

The project structure separates input data, metadata, analysis scripts, and generated results to make the workflow easier to understand and reproduce.

---

## Requirements

The project was developed using Python and the following libraries:

- Python 3.13
- pandas
- NumPy
- SciPy
- Matplotlib
- Seaborn
- scikit-learn
- gseapy
- Anndata
- Statsmodels
- PyDESeq2

The analysis was performed in a Python virtual environment.

---

## Project Status
The complete Python-based RNA-seq differential expression workflow has been implemented and executed successfully on the GSE267442 dataset.

The final analysis includes:

- Differential expression analysis
- PCA analysis
- Volcano plot
- MA plot
- Top-DEG heatmap
- GO Biological Process enrichment
- GO Molecular Function enrichment
- GO Cellular Component enrichment
- KEGG pathway enrichment
- DEG annotation
- Final project summary


## License
This project is intended for academic and educational purposes.