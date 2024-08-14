# DESeq2 Analysis
## Overview
Description of the analysis of RNA-seq data comparing FEN-1 sensitive and FEN-1 resistance cell lines. The analysis includes steps for data import, differential expression analysis, variant analysis result visualization. GSEA was conducted using the GUI based tool found at: https://www.gsea-msigdb.org/gsea/index.jsp. GSEA results were visualized in GUI cytoscape with EnrichMap package.

The results directory contains all final relevant results obtained. Furhter intermediate files can be obtained upon request.

All scripts were included in the "scripts" directory

**Run list**
1) trim_QC.sh
2) STAR.sh
    If not indexed, run STAR_indexing.sh prior running STAR.sh. This indexes the reference genome. 
3) Salmon.sh
    If not indexed, run Salmon-index.sh prior running Salmon.sh. This indexes the reference genome
----> 1-2-3 allow to then run the DeSeq2 analysisi in the R notebook
4) VarScan.sh
5) merge_vcf.sh (merge SNPs and INDEL from VarScan)
6) VCF_annotator tool (see [repo](https://github.com/marcellobeltrami/VCF_annotator) for instructions)
7) MAF_tools.rmd (change inside file paths based on location of step 5 results. This file is generated running the previous step. The exapnded for furhter analysis is also included in repo as MAFTools_og.rmd).

----> this allows to obtain variant calling files which were summarized as plots in the manuscript. These plots were generated using VCF Annotator and MAFTools_og.rmd. 

**Data Import and Preparation and Database Construction**

  - Create or load a txdb database for importing .sf files using tximport. 
  - Define functions to import gene counts and create metadata tables

**Differential Expression Analysis done with DESeq2 Analysis** 

  - Create DESeq2 dataset and run differential expression analysis and Generates volcano plots and filter results.

  - Genes of interest for GBM are then further investigated.

    

**Author**

Marcello Beltrami, Msc Bioinformatics 23/24
