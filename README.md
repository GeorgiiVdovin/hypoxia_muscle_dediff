# Hypoxia inhibits newt skeletal muscle dedifferentiation

Code behind the figures of the manuscript. Open `A1 Analysis.Rproj` (renv restores the packages from `renv.lock`).
Every notebook sets the project root itself and runs top to bottom in a fresh R session. Plots are saved next to the
code that makes them.

Raw sequencing data, quantification files and the Del-Tsonis transcriptome annotation are not in this repository
(see Data availability in the paper). The Notophthalmus A1 single-cell and bulk data can be browsed at
[georgii-vdovin.shinyapps.io/nviridescens-a1-browser](https://georgii-vdovin.shinyapps.io/nviridescens-a1-browser/)
([code](https://github.com/GeorgiiVdovin/nviridescens-a1-browser)).

## Figures and the code behind them

| Panels | Code | Plot files |
| --- | --- | --- |
| Fig. 1A-C, Supp. Fig. 1A, C, D, 2A | `A1_single_cell_data/notebooks/a1_dge.Rmd` | in `A1_single_cell_data/plots`: 1A `Fig1A UMAP Clusters.png`, 1B `Fig1B Markers Dotplot.png`, 1C `Fig1C Marker Violins.png`, S1A `SuppFig1A Cluster Markers Heatmap_v2.png`, S1C `SuppFig1C Cluster GO BP_v2.png`, S1D `SuppFig1D Mapping Rate TE_v2.png`, S2A `SuppFig2A CC Phases per Cluster.png` |
| Supp. Fig. 1B, Supplementary Table 1 | `A1_single_cell_data/notebooks/a1_monocle3.Rmd` | `A1_single_cell_data/plots/SuppFig1B Monocle3 Trajectory_v2.png`, `A1_single_cell_data/results/A1 Trajectory Genes.csv` |
| Fig. 1D, E, Supp. Fig. 2B-E | `A1_cross_species_comparisons/top_tpm_*` | in each folder: 1D `Log2 Median TPM A1 Samples_v2.png`, 1E `Log2 Median TPM 3 C2C12 Samples_v2.png`, S2B `Log2 Median TPM A1 Samples_Supp_v2.png`, S2C `Log2 Median TPM 3 MuSC Samples_v2.png`, S2D `Log2 Median TPM 3 MEF Samples_v2.png`, S2E `Log2 Median TPM 4 Mouse VSMC Samples_v2.png` |
| Fig. 2C, D, F, 3C, E-G, 4C-F, 5A-D, Supp. Fig. 3, Supp. Fig. 4B-F | `A1_diff_dediff.Rmd` (Fig. 2F Myf5, Myog and Myf6 from `sandberg_reference_realigned_bulk/sandberg_muscle_genes.Rmd`) | in the project root: 2C `MT PCA.png`, 2D `Cell Division GO Gene Set.png`, 2F `Muscle Differentiation_v2.png`, 3C `Tricarboxylic Acid Cycle GO Gene Set.png`, 3E `PPARGC1A_v2.png`, 3F `LDHA_v2.png`, 3G `PKM_v2.png`, 4C `BMP-SMAD's_v2.png`, 4D `TGFbeta-SMADs and Downstream Targets_v2.png`, 4E `Overlapping SMAD1 SMAD3_v2.png`, 4F `TGFbeta Signalling Inhibitors_v2.png`, 5A `HIF3A Volcano 8h Dediff 0.001_v2.png`, 5B `Hypoxic Signaling Genes_v2.png`, 5C `HIF3A_v2.png`, 5D `HIF2A_5D.png`, S3A `Top 20 Annotated Unique DEG 8h Dediff.png`, S3B `Cnetplot 0.001 Myotubes 8h Dediff GO BP.png`, S3C `Myotubes Single Clusters 0.001_v2.png`, S3D `0.001 Myotubes 8h Dediff GO BP_v2.png`, S4B `Volcano 8h Dediff 0.001_v2.png`, S4C `MMP GO Gene Set.png`, S4D `Electron Transport Chain GO Gene Set_v2.png`, S4E `HIF1B_v2.png`, S4F `EGLN2_v2.png` |
| Fig. 3A, B, 5E, Supp. Fig. 4A | `IPA_data/IPA_LFC_0.5_p_0.001_plots/ipa_z-score_plots.Rmd` | in the same folder: 3A `TGF-β, BMP and p53 bar.png`, 3B `Metabolic Changes bar.png`, 5E `Cellular Response to Hypoxia bar.png`, S4A `Top 20 IPA Upstream Regulators.png` |
| Statistics of Fig. 4G-I, 5F-H (Supplementary Table 2) | `logistic_regression/logit_link_myotubes.Rmd` | printed only |
| Fig. 6A, C, D | `figure6/figure6.Rmd` | in `figure6/plots`: 6A `A_culture_differentiation.png`, 6C `C_regeneration_heatmap.png`, 6D `D_hif_timecourse.png` |
| Fig. 6B | `figure6/figure6.Rmd`, data from `in_vivo_mouse_muscle_development/dossantos_multiome.Rmd` | `figure6/plots/B_mouse_development.png` |
| Fig. 6E | `figure6/figure6.Rmd`, data from `in_vivo_newt_regeneration/notophthalmus_sc.Rmd` | `figure6/plots/E_blastema_dotplot.png` |

Not code: Prism plots (Fig. 2A, B, E, 4G-I, 5F-H), IPA networks (Fig. 4A, B), metabolic assays (Fig. 3D, H),
microscopy (Fig. 5I, J) and the model (Fig. 6F).

## Folders and notebooks

`A1_diff_dediff.Rmd` - bulk RNA-seq time course of myotube differentiation (0-120h) and dedifferentiation (128-168h):
tximport of the salmon files, DESeq2 with apeglm (128h vs 120h, 168h vs 120h, 120h vs 0h), k-means clusters, GO
enrichment, heatmaps, gene and gene set plots. `salmon_MT_Katia_reference/samples.txt` is the sample sheet.

`sandberg_reference_realigned_bulk/sandberg_muscle_genes.Rmd` - the same time course aligned to the Sandberg transcriptome
(Abdullayev et al., 2013), for Myf5, Myog and Myf6, which are missing from the Del-Tsonis transcriptome. It writes their vst
values to `muscle_diff.csv` for the Fig. 2F heatmap.

`A1_single_cell_data` - single-cell RNA-seq of A1 cells. `a1_preprocessing.Rmd`: QC, removal of contaminating cells,
clustering. `a1_dge.Rmd`: markers, GO, cluster plots. `a1_monocle3.Rmd`: trajectory. `raw_data` holds the Cell Ranger
matrix and the per-cell mapping statistics from the Cell Ranger BAM.

`A1_cross_species_comparisons` - top_tpm directories compare bulk data of A1 cells and mouse cell types (C2C12, MuSC,
MEF, VSMC) to determine the A1 cell type. The mouse notebooks read `tx2gene_mouse.csv` (from gencode.vM35).

`IPA_data` - preparation of the DGE tables for Ingenuity Pathway Analysis (IPA), and the bar plots of IPA pathway z-scores
and the upstream regulator heatmap.

`logistic_regression` - deblinded nuclei counts and binomial GLMs with Tukey contrasts behind the stars of Fig. 4G-I and
5F-H. It prints the results of Supplementary Table 2.

`figure6` - Figure 6 (differentiation and regeneration). The Data chunks pull the plotted values from the
upstream analyses into `figure6/data`; the Panel chunks read only `figure6/data`, so the panels can be rebuilt from this
repository. Sources: mouse C2C12 differentiation (GSE148667, Denes et al., 2019), human myoblast lines AB1167 and AB678
(GSE287954, Donega et al., 2025), KD3 myoblast and myotube nuclei (GSE79920, Zeng et al., 2016), mouse hindlimb muscle
development (GSE211543, Dos Santos et al., 2023), N. viridescens 18 dpa limb blastema (PRJNA904299, Subramanian et al., 2023)
and P. waltl limb regeneration (Elewa et al., 2017; DESeq2 results from Yassouf et al., 2025, figshare
10.6084/m9.figshare.22748393). The Donega transcripts are mapped to genes with the GENCODE v50 annotation and the authors'
SQANTI table (their Supplementary Data 6). The C2C12 sample sheet is
`revision/A1_sc_cross_species/mouse_differentiation_regeneration/GSE148667_C2C12_3/GSE148667_metadata.csv`.

`in_vivo_mouse_muscle_development/dossantos_multiome.Rmd` - mouse hindlimb muscle 10x multiome (E14.5, E18.5, P5, adult;
Dos Santos et al., 2023, GSE211543, deposited Cell Ranger ARC matrices), processed following the authors' script: QC,
SCTransform and RPCA integration of the RNA, LSI of the ATAC peaks, weighted nearest neighbours, cell types from the paper's
markers and chromVAR motif activity (JASPAR 2024). As the fragment files are not deposited, the peak counts of each library
are summed over the peaks called in all four libraries. It writes the myogenic nuclei and their motif activity for Fig. 6B.

`in_vivo_newt_regeneration` - download, quantification (kallisto) and processing of the N. viridescens blastema
single-cell data used in Fig. 6E.

`Nvir_Annotations` - annotation of the Del-Tsonis transcriptome (Tsissios et al., 2023): mitochondrial genes for
single-cell QC, blastp annotation against Swiss-Prot with GO terms, and the org.Nvi.eg.db package for clusterProfiler.

`tx2gene_katia` - transcript to gene mapping for tximport, also used for org.Nvi.eg.db.
