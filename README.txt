A1 Single Cell Data - processing and ShinyApp of A1 SC data from Fig 1.

A1_cross_species_comparisons - top_tpm directories are analysis of bulk data from different mouse cell lines, 
to determine A1 cell type (Fig 1). 

Combined Annotations Nvir - combined transcriptome annotation tables from 3 assemblies (Sandberg lab (Abdullayev et 
al., 2013), Simon lab (Subramanian et al., 2024), Del-Rio Tsonis lab (Tsissios et al., 2023)) for exploration.

IPA_data - scripts to prepare DGE tables for input into Ingenuity Pathway Analysis (IPA), and scripts to generate 
heatmaps using IPA-aquired data from Fig. 3,5 & Supp. Fig. 3.

logistic_regression

Katia_Reference_OG - original, unmodified Del-Tsonis transcriptome and annotation (Tsissios et al., 2023).

Nvir_Annotations - code to create the annotation for the Del-Tsonis transcriptome, introducing mitochondrial genes
for SC RNAseq quality control, blastp annotations against Swissprot with GO terms, as well as generation of the org.Nvi.eg.db 
annotation package for use with ClusterProfiler and other enrichment packages. "Nvir Latest..." folder contains the same 
annotation in a separate directory.

salmon_MT_Katia_reference - salmon quantification files from aligning the myotube differentiation timecourse reads against
the Del-Tsonis transcriptome using Galaxy.

tx2gene_katia - transcript to gene mapping file for importing pseudocounts using tximport. Also used in generating the 
org.Nvi.eg.db annotation package.