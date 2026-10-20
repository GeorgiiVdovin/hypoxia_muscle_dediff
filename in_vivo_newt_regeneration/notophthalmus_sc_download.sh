#!/usr/bin/env bash
# Notophthalmus adult limb blastema at 18 dpa, Smart-seq2 (Subramanian et al. 2023, SRP409251 / PRJNA904299).
# Downloads the 768 single-cell FASTQs (single-end, 24 GB) from ENA into WSL, outside Nextcloud. Re-running resumes.
# 2 runs are broken on ENA's side and fail their md5 in every round; notophthalmus_sc_quant.sh skips them
# Run from Windows:  wsl -e bash "/mnt/c/.../A1 Analysis/in_vivo_newt_regeneration/notophthalmus_sc_download.sh"
set -euo pipefail

W=~/nvir_blastema # working folder in WSL, also used by notophthalmus_sc_quant.sh
mkdir -p $W/fastq
cd $W

# Run list with ENA FASTQ links, md5s and sizes, and the md5sum check file ("<md5>  fastq/<run>.fastq.gz")
curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA904299&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv&limit=0" |
  tail -n +2 > runs.tsv
awk -F'\t' '{print $3"  fastq/"$1".fastq.gz"}' runs.tsv > md5.txt

# Fetch whatever is missing or fails its md5, up to 5 rounds, 4 in parallel (ENA drops connections under more)
for round in 1 2 3 4 5; do
  md5sum -c md5.txt 2>/dev/null | grep -v ": OK$" | sed 's|^fastq/||; s|.fastq.gz:.*||' > todo.txt || true # run IDs to fetch
  [ -s todo.txt ] || break
  echo "Round $round: $(wc -l < todo.txt) files to fetch"

  # Full-size files that fail md5 are corrupt and start again, shorter ones resume (-C -). -L: some runs moved on ENA (301)
  grep -Ff todo.txt runs.tsv | while IFS=$'\t' read -r id url md5 bytes; do
    f=fastq/$id.fastq.gz
    [ -f $f ] && [ $(stat -c %s $f) -ge $bytes ] && rm $f
    echo "https://$url"
  done | xargs -P 4 -n 1 curl -sSL --retry 10 --retry-all-errors --retry-delay 10 -C - --output-dir fastq -O || true
done

md5sum -c --quiet md5.txt && echo "All $(wc -l < runs.tsv) FASTQs downloaded, md5 OK"
