# AP-1 Primes Transcriptional Noise-Led Cancer Treatment Adaptation Through KDM5 Chromatin Rewiring.


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
