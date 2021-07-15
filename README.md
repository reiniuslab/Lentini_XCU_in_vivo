# Allele-resolved single-cell multi-omics uncovers the dynamics and transcriptional kinetics of X-chromosome upregulation
*Antonio Lentini<sup>1</sup>, Huaitao Cheng<sup>2</sup>, Joyce Carol Noble<sup>1</sup>, Natali Papanicolaou<sup>1</sup>, Christos Coucoravas<sup>1</sup>, Nathanael Andrews<sup>2</sup>, Qiaolin Deng<sup>3</sup>, Martin Enge<sup>2</sup> and Björn Reinius<sup>1</sup>*

<sup>1</sup>Department of Medical Biochemistry and Biophysics, Karolinska Institutet, Stockholm, Sweden.
<sup>2</sup>Department of Oncology and Pathology, Karolinska Institutet, Stockholm, Sweden.
<sup>2</sup>Department of Physiology and Pharmacology, Karolinska Institutet, Stockholm, Sweden.
##

This repository contains scripts and additional data needed to reproduce findings described in [PAPER](), [PRE-PRINT](https://www.biorxiv.org/content/10.1101/2021.07.14.452323v1).

Raw and processed data is available through ArrayExpress: [[Smart-seq3]](), [[Joint Smart-seq3+scATAC]]() and [[Allelic dilution series]]().

``x_escapees.tsv`` is compiled from ([1](https://dx.doi.org/10.1038/ng.3678), [2](https://dx.doi.org/10.1101/gr.103200.109), [3](https://dx.doi.org/10.1038/nsmb.3365), [4](https://dx.doi.org/10.1186/1471-2164-11-614))

``xy_homologs.tsv`` is modified from Table1 in ([5](https://dx.doi.org/doi:10.1016/j.cell.2014.09.052)) 

## Abstract
X-chromosome inactivation (XCI) and upregulation (XCU) are the major opposing chromosome-wide modes of gene regulation that collectively achieve dosage compensation in mammals, but the regulatory link between the two remains elusive. Here, we use allele-resolved single-cell RNA-seq combined with chromatin accessibility profiling to finely dissect the separate effects of XCI and XCU on RNA levels during mouse development. We uncover that balanced X dosage is flexibly attained through expression tuning by XCU in a sex- and lineage-specific manner along varying degrees of XCI and across developmental and cellular states. Male blastomeres achieve XCU upon zygotic genome activation while females experience two distinct waves of XCU, upon imprinted- and random XCI, and ablation of Xist impedes female XCU. Contrary to widely established models of mammalian dosage compensation, naïve female embryonic cells carrying two active X chromosomes do not exhibit upregulation but express both alleles at basal level, yet collectively exceeding the RNA output of a single hyperactive allele. We show, in vivo and in vitro, that XCU is kinetically driven by X-specific modulation of transcriptional burst frequency, coinciding with increased compartmentalization of the hyperactive allele. Altogether, our data provide unprecedented insights into the dynamics of mammalian XCU, prompting a revised model of the chain in events of allelic regulation by XCU and XCI in unitedly achieving stable cellular levels of X-chromosome transcripts.
##

### Session info
*R version 3.6.1 (2019-07-05)*

| Package | Version | Source |
| :------ | :------ | :----- |
| `ArchR` | 1.0.1 | Github (GreenleafLab/ArchR@968e442) |
| `Biobase` | 2.44.0 | Bioconductor |
| `BiocGenerics` | 0.30.0 | Bioconductor |
| `BiocParallel` | 1.18.1 | Bioconductor |
| `biomaRt` | 2.40.3 | Bioconductor |
| `boot` | 1.3-23 | CRAN (R 3.6.1) |
| `circlize` | 0.4.6 | CRAN (R 3.6.1) |
| `ComplexHeatmap` | 2.0.0 | Bioconductor |
| `cowplot` | 1.0.0 | CRAN (R 3.6.1) |
| `data.table` | 1.12.2 | CRAN (R 3.6.1) |
| `DelayedArray` | 0.10.0 | Bioconductor |
| `DESeq2` | 1.24.0 | Bioconductor |
| `GenomeInfoDb` | 1.20.0 | Bioconductor |
| `GenomicRanges` | 1.36.0 | Bioconductor |
| `ggbeeswarm` | 0.6.0 | CRAN (R 3.6.1) |
| `ggplot2` | 3.2.1 | CRAN (R 3.6.1) |
| `hdf5r` | 1.3.2.9000 | Github (hhoeflin/hdf5r@d38b053) |
| `IRanges` | 2.18.1 | Bioconductor |
| `loomR` | 0.2.1.9000 | Github (mojaveazure/loomR@1eca16a) |
| `magrittr` | 1.5 | CRAN (R 3.6.1) |
| `MAST` | 1.10.0 | Bioconductor |
| `matrixStats` | 0.54.0 | CRAN (R 3.6.1) |
| `mclust` | 5.4.5 | CRAN (R 3.6.1) |
| `princurve` | 2.1.4 | CRAN (R 3.6.1) |
| `R6` | 2.4.0 | CRAN (R 3.6.1) |
| `RColorBrewer` | 1.1-2 | CRAN (R 3.6.0) |
| `rtracklayer` | 1.44.2 | Bioconductor |
| `S4Vectors` | 0.22.0 | Bioconductor |
| `scater` | 1.12.2 | Bioconductor |
| `scran` | 1.12.1 | Bioconductor |
| `SingleCellExperiment` | 1.6.0 | Bioconductor |
| `slingshot` | 1.2.0 | Bioconductor |
| `SummarizedExperiment` | 1.14.1 |  Bioconductor |
| `tximport` | 1.12.3 | Bioconductor |
