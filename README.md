
# RNA-SEQ DATA ANALYSIS PIPELINE

This is a reproducible end-to-end bioinformatics pipeline for RNA sequencing data analysis. This pipeline involves the following steps:

- Downloaded raw counts from GEO (GSE60450)
- Cleaned counts matrix by removing non-count columns
- Filtered low expression genes
- Built sample metadata table
- Ran DESeq2 differential expression analysis
- Extracted 4 pairwise comparisons
- VST normalized counts for visualization
- Constructed volcano plots for all 4 comparisons
- Constructed PCA plot
- Generated heatmap of the top 50 variable genes
- Wrapped in Snakemake for reproducibility


## Background:
This project aims to answer biological questions based on the publicly available RNA-seq data GSE60450 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450) from Fu et al. (2015, Nature Cell Biology). The study highlights the important role of the pro-survival protein Mcl-1 in the development and function of the mammary gland. Mcl-1 is a member of the Bcl-2 family of proteins and is an anti-apoptotic protein and a key regulator of cell survival. Mcl-1 has been shown to be indispensable at several developmental stages and it is essential for morphogenesis during puberty, stem and progenitor cell activity, alveologenesis during pregnancy and alveolar cell survival during lactation.
This analysis focuses on comparison of genes between basal and luminal cells under three conditions: virgin, pregnant and lactating. 


## Biological Questions:
This analysis aimed to address the following biological questions:
1. Which genes are differentially expressed between basal and luminal mammary epithelial cells, independent of developmental stage?
2. How does the transcriptional landscape of mammary epithelial cells change across developmental stages (virgin → pregnant → lactating), independent of cell type?

## Data:
The data used for this pipeline is GSE60450 available publicly on Gene Expression Omnibus (GEO). It is a transcriptome analysis of the luminal and basal cell subpopulations in the lactating versus pregnant mammary glands derived from Mus musculus. It contains 12 samples of luminal and basal cells harvested from the mammary glands of virgin, 18.5 day pregnant and 2 day lactating mice with 2 mice per stage. The data was derived from expression profiling by high throughput sequencing.

## Requirements:
- conda or mamba
- All dependencies are specified in `envs/rnaseq.yaml`

## Setup:
```bash
# Clone the repository
git clone <your-repo-url>
cd rnaseq-pipeline

# Create and activate the environment
conda env create -f envs/rnaseq.yaml
conda activate rnaseq

# Download the data
# Go to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60450
# Download GSE60450_Lactation-GenewiseCounts.txt.gz
# Place in data/raw/ and unzip with gunzip
```

## How to Run:
```bash
# Full pipeline
snakemake --cores 1

# Dry run (check without executing)
snakemake --dry-run

# Force rerun everything
snakemake --cores 1 --forceall
```

## Project Structure:
```
rnaseq-pipeline/
├── config/
│ └── config.yaml # Sample metadata and design formula
├── data/
│ └── raw/ # Raw counts file (not tracked by git)
├── envs/
│ └── rnaseq.yaml # Conda environment specification
├── results/
│ ├── figures/ # PCA, volcano plots, heatmap
│ └── *.csv # DESeq2 results tables
├── scripts/
│ ├── deseq2.R # DESeq2 analysis (R)
│ └── plotting.py # Visualization (Python)
├── Snakefile # Workflow definition
└── README.md
```

## Tools Used
| Tool | Language | Purpose |
|---|---|---|
| DESeq2 | R | Differential expression analysis |
| yaml | R | Config file parsing |
| ggplot2 | R | Plotting (R side) |
| pheatmap | R | Heatmap (R side) |
| Snakemake | Python | Workflow management |
| pandas | Python | Data manipulation |
| seaborn | Python | Statistical visualization |
| scikit-learn | Python | PCA |
| matplotlib | Python | Plotting |
| numpy | Python | Numerical computing |


## References
1. Fu NY, Rios AC, Pal B, Soetanto R, Lun AT, Liu K, Beck T, Best SA,
   Vaillant F, Bouillet P, Strasser A, Preiss T, Smyth GK, Lindeman GJ,
   Visvader JE. (2015). EGF-mediated induction of Mcl-1 at the switch to
   lactation is essential for alveolar cell survival. Nature Cell Biology,
   17, 365–375. https://doi.org/10.1038/ncb3117
   PMID: 25730472

2. Visvader JE. (2007). Origins of breast cancer subtypes and therapeutic
   implications. Nature Clinical Practice Oncology, 4(9), 516–525.
   https://doi.org/10.1038/ncponc0908
   PMID: 17728710

## Author
Jayasri Ramesh

