# TNBC RNA-seq Differential Expression Analysis

## Project Overview

This project performs a Python-based RNA-seq differential expression analysis of **Triple-Negative Breast Cancer (TNBC)** compared with **paracancerous breast tissue**.

The analysis starts from a gene-level count matrix and metadata, followed by data preprocessing, differential expression analysis, visualization, functional enrichment, and annotation of significant genes.

The project demonstrates a reproducible computational biology workflow using Python and commonly used bioinformatics analysis approaches.

---

## Research Objective

The primary objective of this project is to identify genes that are differentially expressed between:

- **TNBC tissue**
- **Paracancerous breast tissue**

The identified differentially expressed genes (DEGs) are further investigated to understand their associated biological functions and pathways.

---

## Analysis Workflow

```text
Gene Count Matrix

        |

        v

Data Loading & Quality Checks

        |

        v

Metadata Preparation

        |

        v

Differential Expression Analysis

        |

        v

DEG Categorization

        |

        +------------------+
        |                  |
        v                  v

       PCA          Volcano / MA Plot

        |                  |
        +--------+---------+
                 |
                 v

          Top DEGs / Heatmap

                 |

                 v

        Functional Enrichment

           /           \

          v             v

         GO            KEGG

     BP / MF / CC

          |

          v

     DEG Annotation

          |

          v

   Final Project Summary
```

## Project Structure

```text
TNBC_RNAseq_Python/
│
├── data/
│   ├── counts/
│   │   └── gene_counts.txt
│   │
│   └── metadata/
│       └── metadata.csv
│
├── results/
│   ├── GO_enrichment/
│   │   ├── GO_BP/
│   │   ├── GO_CC/
│   │   ├── GO_MF/
│   │   ├── GO_BP_results.csv
│   │   ├── GO_CC_results.csv
│   │   ├── GO_MF_results.csv
│   │   └── DEG_gene_symbols.csv
│   │
│   ├── KEGG_enrichment/
│   │   ├── KEGG_results.csv
│   │   └── KEGG_top10.png
│   │
│   ├── PCA_plot.png
│   ├── MA_plot.png
│   ├── Volcano_plot.png
│   ├── Heatmap_Top50_DEGs.png
│   ├── PCA_coordinates.csv
│   ├── Top50_DEGs.csv
│   ├── significant_DEGs.csv
│   ├── differential_expression_results.csv
│   ├── differential_expression_categorized.csv
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
├── .gitignore
├── requirements.txt
└── README.md
```

The large raw count file is excluded from version control using `.gitignore`.

---

## Analysis Steps

### 1. Count Matrix Preparation

The analysis begins with a gene-level RNA-seq count matrix.

The count data are loaded using Python and Pandas and prepared for downstream differential expression analysis.

Large raw sequencing/count files are excluded from Git version control using `.gitignore`.

---

### 2. Metadata Preparation

Sample metadata are prepared to define the experimental conditions:

- **TNBC**
- **Paracancerous breast tissue**

The metadata are used to associate each sample with its corresponding biological condition.

---

### 3. Differential Expression Analysis

Differential expression analysis is performed to identify genes whose expression differs between TNBC and paracancerous breast tissue.

The analysis produces:

- Complete differential expression results
- Significant DEGs
- Categorized DEGs
- Annotated DEGs

---

### 4. Principal Component Analysis

Principal Component Analysis (PCA) is used to visualize the overall structure of the samples and assess separation between TNBC and paracancerous breast tissue samples.

**Output:**

`results/PCA_plot.png`

---

### 5. Volcano Plot and MA Plot

Volcano and MA plots are generated to visualize differential expression results.

**Outputs:**

- `results/Volcano_plot.png`
- `results/MA_plot.png`

---

### 6. Heatmap

A heatmap of the top 50 differentially expressed genes is generated to visualize expression patterns across samples.

**Output:**

`results/Heatmap_Top50_DEGs.png`

---

### 7. Gene Ontology Enrichment

Significant DEGs are investigated using Gene Ontology (GO) enrichment analysis.

Three GO categories are considered:

- **Biological Process (BP)**
- **Molecular Function (MF)**
- **Cellular Component (CC)**

The enrichment results and visualization files are available under:

`results/GO_enrichment/`

---

### 8. KEGG Pathway Enrichment

KEGG pathway enrichment is performed to identify biological pathways associated with the identified DEGs.

**Output:**

`results/KEGG_enrichment/`

---

### 9. DEG Annotation

Significant genes are further annotated to obtain gene symbols and generate final annotated DEG tables.

---

## Results

The final differential expression analysis identified:

- **47,726 genes analyzed**
- **3,559 significant DEGs**
- **1,878 upregulated genes**
- **1,681 downregulated genes**

### Top 10 Upregulated Genes

| Gene | log2FC |
|---|---:|
| IVL | 9.3301 |
| IBSP | 8.1468 |
| CST1 | 7.9489 |
| TMEM151A | 7.6020 |
| MMP1 | 7.4607 |
| CLEC6A | 7.4554 |
| MMP13 | 7.4031 |
| COPDA1 | 7.2562 |
| IL21-AS1 | 7.2189 |
| LOC102723635 | 7.1883 |

### Top 10 Downregulated Genes

| Gene | log2FC |
|---|---:|
| ENSG00000286208 | -8.9361 |
| LINC01087 | -8.7804 |
| ARHGAP36 | -8.3847 |
| CT62 | -8.2681 |
| POTEKP | -8.2277 |
| SERPINA6 | -7.7044 |
| ELOVL2-AS1 | -7.5429 |
| ENSG00000301521 | -7.4885 |
| LOC107986528 | -7.2653 |
| MARCHF11 | -7.1891 |

---

## Results & Visualizations

The analysis generated multiple visualizations to assess sample structure, differential expression, and functional enrichment.

### Principal Component Analysis

PCA was used to visualize sample-level variation and assess separation between TNBC and paracancerous breast tissue samples.

![PCA Plot](results/PCA_plot.png)

### Differential Expression

The Volcano Plot and MA Plot summarize the differential expression results.

#### Volcano Plot

![Volcano Plot](results/Volcano_plot.png)

#### MA Plot

![MA Plot](results/MA_plot.png)

### Top Differentially Expressed Genes

A heatmap of the top 50 differentially expressed genes was generated to visualize expression patterns across samples.

![Top 50 DEGs Heatmap](results/Heatmap_Top50_DEGs.png)

### Gene Ontology Enrichment

GO enrichment analysis was performed across:

- **Biological Process (BP)**
- **Molecular Function (MF)**
- **Cellular Component (CC)**

#### Biological Process

![GO Biological Process](results/GO_enrichment/GO_BP_top10.png)

#### Molecular Function

![GO Molecular Function](results/GO_enrichment/GO_MF_top10.png)

#### Cellular Component

![GO Cellular Component](results/GO_enrichment/GO_CC_top10.png)

### KEGG Pathway Enrichment

KEGG enrichment analysis was performed to identify pathways associated with the identified DEGs.

![KEGG Enrichment](results/KEGG_enrichment/KEGG_top10.png)

---

## Analysis Outputs

The following files are generated by the analysis.

### Differential Expression

- `results/differential_expression_results.csv`
- `results/differential_expression_categorized.csv`
- `results/significant_DEGs.csv`

### Annotation

- `results/final_annotated_DEGs.csv`
- `results/final_significant_DEGs_annotated.csv`

### PCA

- `results/PCA_coordinates.csv`
- `results/PCA_plot.png`

### Differential Expression Visualization

- `results/Volcano_plot.png`
- `results/MA_plot.png`

### Heatmap

- `results/Top50_DEGs.csv`
- `results/Heatmap_Top50_DEGs.png`

### GO Enrichment

- `results/GO_enrichment/GO_BP_results.csv`
- `results/GO_enrichment/GO_CC_results.csv`
- `results/GO_enrichment/GO_MF_results.csv`
- `results/GO_enrichment/GO_BP_top10.png`
- `results/GO_enrichment/GO_CC_top10.png`
- `results/GO_enrichment/GO_MF_top10.png`

### KEGG Enrichment

- `results/KEGG_enrichment/KEGG_results.csv`
- `results/KEGG_enrichment/KEGG_top10.png`

### Final Summary

- `results/final_project_summary.txt`

---

## Technologies & Tools

### Programming

- Python
- Pandas
- NumPy
- Matplotlib
- Seaborn

### Bioinformatics & Statistical Analysis

- RNA-seq differential expression analysis
- Principal Component Analysis (PCA)
- Gene Ontology (GO) enrichment
- KEGG pathway enrichment
- DEG annotation

### Development Environment

- Visual Studio Code
- Python virtual environment
- Git
- GitHub

---

## Reproducibility

The analysis has been organized into sequential Python scripts following a structured workflow.

```text
01 → Load count data
02 → Prepare metadata
03 → Differential expression analysis
04 → PCA
05 → Volcano and MA plots
06 → Heatmap
07 → GO Biological Process enrichment
08 → GO Molecular Function & Cellular Component enrichment
09 → KEGG enrichment
10 → DEG annotation
11 → Final project summary
```

Each step produces outputs that are used by subsequent stages of the analysis.

---

## Project Significance

This project demonstrates an end-to-end RNA-seq bioinformatics workflow, from gene-level count data through differential expression analysis and biological interpretation.

The workflow integrates:

- Data preprocessing
- Experimental metadata handling
- Differential expression analysis
- Statistical visualization
- Gene-level annotation
- Gene Ontology enrichment
- KEGG pathway analysis
- Reproducible project organization

The project provides practical experience with computational analysis of transcriptomic data and demonstrates how differential gene expression can be investigated in the context of TNBC biology.

---

## Author

**Suraj Kumar Chanda**

GitHub: **shurajkr1412**

---

## Disclaimer

This project is intended for educational and research portfolio purposes and demonstrates an RNA-seq data analysis workflow using the available project dataset.
