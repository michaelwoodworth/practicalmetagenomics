# Practical Metagenomics
The primary motivation for this repository to document reading outlines, resources, and questions as an introduction to practical analysis of metagenomes.

At each step in analysis, we post the tools, their citation (when available), and annotated code used with (in most cases) the expected input and output for metagenomic data. We use the genomic and metagenomic analyses from the paper **Randomized trial: Fecal microbiota transplant promotes reduction of antimicrobial resistance by strain replacement** as a practical example. 

Data for the analyses in the PREMIX paper that will be used for this repository are available through the NCBI Bioproject (accession PRJNA728680).

Subsets of these data and other files for this course are on [OneDrive](https://emory-my.sharepoint.com/:f:/g/personal/mwoodwo_emory_edu/EiKUkQ__b2lLjGX8h-VmZZAB5Hx1y0kPhSvmcwR59xg97g?e=eZk9ad), please email me for access if you have trouble.

---

## Overview
We will meet weekly from January 27 until April 28 to cover the following topics:

| Week | Dates | Topic | Paper(s)/Link |
| --- | --- | --- | --- |
| 1	| Jan 22 - Jan 28 | [Understand Premix dataset and log in to server](pages/23.01.27.md) | Woodworth et al (under review)
| 2	| Jan 29 - Feb 4 | [Metagenome quality control](pages/23.02.03.md) | [PMID: 24695404 (Trimmomatic)](https://pubmed.ncbi.nlm.nih.gov/24695404/)
| 3	| Feb 5 - Feb 11 | [Short read taxonomic classification](pages/23.02.10.md) | [Kraken2](https://pubmed.ncbi.nlm.nih.gov/31779668/) and [bracken](https://peerj.com/articles/cs-104/)
| 4	| Feb 12 - Feb 18 | [Intro to diversity metrics](pages/23.02.17.md) | ([Vegan](https://github.com/vegandevs/vegan) and [Phyloseq](https://joey711.github.io/phyloseq/index.html) packages in R)
| 5	| Feb 19 - Feb 25 | [Dimensionality reduction (PCA, PCoA, TSNE)](pages/23.02.24.md) |	https://joey711.github.io/phyloseq/index.html
| 6 | Feb 26 - Mar 4 | [Contig Assembly and binning](pages/23.03.03.md)	| [PMID: 28298430](https://doi.org/10.1101/gr.213959.116)
| 7	| Mar 5 - March 11 | [Contig gene prediction and annotation](pages/23.03.10.md)	| PMID: 20211023 PMID: 34597405
| 8	| Mar 12 - Mar 18 | Metagenomic plotting and data vis	| https://ggplot2.tidyverse.org
| 9	| Mar 19 - Mar 25 | Analysis with QIIME 2	| https://qiime2.org
| 10 | Mar 26 - Apr 1 | Analysis with QIIME 2	| https://qiime2.org
| 11 | Apr 2 - Apr 8 | How to work with 16S data and perform 16S sequencing| 	https://qiime2.org
| 12 | Apr 9 - Apr 15 | Final Project and intro to Python | 	https://swcarpentry.github.io/python-novice-inflammation/
| 13 | Apr 16 - Apr 22 | Final Project
| 14 | Apr 23 - April 29 | Final Project
