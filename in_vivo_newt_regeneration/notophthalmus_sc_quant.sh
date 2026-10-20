#!/usr/bin/env bash
# Notophthalmus adult limb blastema at 18 dpa, Smart-seq2 (Subramanian et al. 2023, SRP409251), quantified against the
# Katia transcriptome and tx2gene used for the A1 bulk data. Needs notophthalmus_sc_download.sh to have run.
# kb_python 0.30.2 (with its kallisto 0.52.0) in a micromamba Python 3.11 env, ~/nvir_blastema/env.
# Run from Windows, with RStudio closed (indexing peaks at ~10 GB RAM, ~60 GB disk in WSL):
#   wsl -e bash "/mnt/c/.../A1 Analysis/in_vivo_newt_regeneration/notophthalmus_sc_quant.sh"
set -euo pipefail

W=~/nvir_blastema
P="$(cd "$(dirname "$0")/.." && pwd)" # project root (A1 Analysis)
OUT="$P/in_vivo_newt_regeneration/data/notophthalmus_sc"
KALLISTO=$W/env/lib/python3.11/site-packages/kb_python/bins/linux/kallisto/kallisto
cd $W

# 1. Use every FASTQ that passes its ENA md5 and report the rest (2 runs are broken on ENA's side)
(md5sum -c md5.txt 2>/dev/null || true) | grep ": OK$" | sed 's|^fastq/||; s|.fastq.gz: OK$||' | sort > ok.txt
cut -f1 runs.tsv | sort | comm -13 ok.txt - > skipped.txt
echo "Using $(wc -l < ok.txt) of $(wc -l < runs.tsv) cells, skipped: $(tr '\n' ' ' < skipped.txt)"

# 2. Index the Katia transcriptome once, t2g from the Katia tx2gene
[ -f katia.fa ] || cp "$P/Katia_Reference_OG/NVir_transcriptome_katia.fasta" katia.fa
tail -n +2 "$P/tx2gene_katia/tx2gene_nvir_katia.csv" | tr -d '"' | awk -F, -v OFS='\t' '{print $1, $2, $2}' > t2g.txt
[ -f katia.idx ] || $KALLISTO index -i katia.idx -t 14 katia.fa

# 3. Pseudoalign all cells in one batch, then EM per cell (single-end, salmon's default fragment length).
#    kb names batch sample i by its base-4 index over 16 nt (AAAAAAAAAAAAAAAA, ...C, ...G, ...T, ...CA); on-listing all of
#    them keeps every cell, otherwise the bustools whitelist knee silently drops the low-read ones
awk -v W=$W '{print $1 "\t" W "/fastq/" $1 ".fastq.gz"}' ok.txt > batch.txt
awk -v n=$(wc -l < ok.txt) 'BEGIN{split("A C G T", b, " "); for (i = 0; i < n; i++) {s = ""; k = i; for (p = 0; p < 16; p++) {s = b[k % 4 + 1] s; k = int(k / 4)}; print s}}' > onlist.txt
$W/env/bin/kb count -i katia.idx -g t2g.txt -x SMARTSEQ2 --parity single --tcc --fragment-l 250 --fragment-s 25 \
  -w onlist.txt -t 14 -m 8G -o kb_out batch.txt

# 4. Gene-level estimated counts and TPM per cell back to the project (rows = cells in matrix.cells, columns = genes.txt).
#    Matrix rows are the barcodes with reads (cells without any are dropped), so map them back to run IDs in matrix order
paste kb_out/matrix.sample.barcodes kb_out/matrix.cells > bc2cell.tsv
awk 'NR==FNR{n[$1]=$2; next} {print n[$1]}' bc2cell.tsv kb_out/counts_unfiltered/cells_x_tcc.barcodes.txt > cells_in_matrix.txt
mkdir -p "$OUT"
cp kb_out/quant_unfiltered/matrix.abundance.gene.mtx kb_out/quant_unfiltered/matrix.abundance.gene.tpm.mtx \
   kb_out/quant_unfiltered/genes.txt kb_out/run_info.json skipped.txt "$OUT/"
cp cells_in_matrix.txt "$OUT/matrix.cells"
echo "Done: $(wc -l < cells_in_matrix.txt) of $(wc -l < ok.txt) cells had reads, $(grep p_pseudoaligned kb_out/run_info.json)"

# 5. ERCC spike-ins per cell, for the ERCC QC in notophthalmus_sc.Rmd. ERCC92 sequences from Thermo Fisher
[ -f ERCC92.fa ] || { curl -sL -o ERCC92.zip https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip; python3 -m zipfile -e ERCC92.zip .; }
grep '^>' ERCC92.fa | sed 's/^>//' | awk -v OFS='\t' '{print $1, $1, $1}' > ercc_t2g.txt
[ -f ercc.idx ] || $KALLISTO index -i ercc.idx ERCC92.fa
$W/env/bin/kb count -i ercc.idx -g ercc_t2g.txt -x SMARTSEQ2 --parity single -w onlist.txt -t 14 -o kb_ercc batch.txt
awk 'NR==FNR{n[$1]=$2; next} {print n[$1]}' bc2cell.tsv kb_ercc/counts_unfiltered/cells_x_genes.barcodes.txt > "$OUT/ercc.cells"
cp kb_ercc/counts_unfiltered/cells_x_genes.mtx "$OUT/ercc.mtx"
