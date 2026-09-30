# AP-1 Primes Transcriptional Noise-Led Cancer Treatment Adaptation Through KDM5 Chromatin Rewiring.


## Abstract

How drug-tolerant persister (DTP) cells escape quiescence to drive tumour relapse remains a central unresolved question in cancer evolution. Here, we use neuroblastoma and hepatoblastoma, paediatric cancers with low genetic instability, as tractable model systems to study the principles and molecular mechanisms that govern the earliest stages of adaptation to treatment. We identify transcriptional noise, stochastic variability in gene expression, as a latent cellular property amplified post treatment, during awakening from persistence, and drives adaptive regrowth after therapy withdrawal.

Through functional assays, lineage tracing and single-cell multiomics, we show that treatment enriches for progenitor-like tolerant states with limited detectable clonal selection under the resolution of our barcode system. Upon drug withdrawal, awakening proceeds via a stochastic process driven by noise-enabled plasticity in cell identity programmes.

Mechanistically, AP-1 activity primes KDM5A recruitment to cell-state identity genes during awakening, driving H3K4me3 peak narrowing at their promoters and the transcriptional-noise surge that enables MES-to-ADRN plasticity. Pharmacological or genetic suppression of AP-1 or KDM5 suppresses transcriptional noise and phenotypic plasticity, blocks DTP exit and prevents tumour recovery in both neuroblastoma and hepatoblastoma models.

Together, these results identify KDM5-mediated transcriptional noise, operating within AP-1-engaged chromatin landscapes, as a key driver of early adaptation to treatment. Targeting this axis offers a strategy to exploit DTP awakening as an evolutionary bottleneck to prevent treatment adaptation and future cancer relapse.

## Code

Analysis code for the manuscript. Each folder has its own README.

The following analyses use the pipelines published in
[PCEBrunaLab/Roux-et-al.-2024](https://github.com/PCEBrunaLab/Roux-et-al.-2024):

- scRNA-seq pre-processing: [`Single_Cell_RNA_Sequencing`](https://github.com/PCEBrunaLab/Roux-et-al.-2024/tree/main/Single_Cell_RNA_Sequencing)
- MuTrans analysis: [`MuTrans_Analysis`](https://github.com/PCEBrunaLab/Roux-et-al.-2024/tree/main/MuTrans_Analysis) (untreated and cisplatin conditions)
- Cellecta barcode analysis (bulk DNA): [`DNA_Cellecta_Barcodes`](https://github.com/PCEBrunaLab/Roux-et-al.-2024/tree/main/DNA_Cellecta_Barcodes)

Pre-processing of the PDX bulk RNA-seq data used the Institut Curie RNA-seq pipeline
([bioinfo-pf-curie/RNA-seq v4.1.0](https://github.com/bioinfo-pf-curie/RNA-seq/tree/v4.1.0);
[doi:10.5281/zenodo.13744441](https://doi.org/10.5281/zenodo.13744441)), with the `--trimming` and `--pdx` options.

| Folder | Content |
|---|---|
| `BASiCS/` | Transcriptional noise quantification from scRNA-seq (BASiCS) analysis |
| `Multiome/` | Single-nucleus multiome (snRNA-seq and snATAC-seq) analysis |
| `snRNA_clinical/` | snRNA-seq of clinical neuroblastoma samples |
| `snRNAseq_cell_lines/` | snRNA-seq of SK-N-SH and HuH6 treated with DMSO, KDM5-C70 or T-5224 |
| `CUTTag/` | CUT&Tag processing and analysis |
| `CUTRUN/` | CUT&RUN processing and analysis |
