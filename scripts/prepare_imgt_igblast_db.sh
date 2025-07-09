#!/bin/bash

# Exit on error
set -e

# Set working directory
mkdir -p imgt_refs && cd imgt_refs
BLASTPATH="/home/ckibet@iavi.org/repos/ncbi-igblast-1.22.0/bin/"
# Define IMGT FASTA URLs
# declare -A urls=(
#   [IGHV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta"
#   [IGHD]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta"
#   [IGHJ]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHJ.fasta"
#   [IGKV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKV.fasta"
#   [IGLV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLV.fasta"
#   [IGLJ]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLJ.fasta"
# )

# echo "[Step 1] Downloading IMGT FASTA files..."
# for gene in "${!urls[@]}"; do
#     echo "Downloading ${gene}..."
#     curl -s -O "${urls[$gene]}"
# done

# Check for required tools
if ! command -v $BLASTPATH/edit_imgt_file.pl &> /dev/null; then
    echo "[ERROR] edit_imgt_file.pl not found. Please add IgBLAST's bin directory to your PATH."
    exit 1
fi

if ! command -v $BLASTPATH/makeblastdb &> /dev/null; then
    echo "[ERROR] makeblastdb not found. Please install BLAST+ tools and add to PATH."
    exit 1
fi

echo "[Step 2] Cleaning headers with edit_imgt_file.pl..."
for gene in "${!urls[@]}"; do
    in_file="${gene}.fasta"
    out_file="${gene}_igblast.fasta"
    echo "Processing $in_file -> $out_file"
    perl $BLASTPATH/edit_imgt_file.pl "$in_file" > "$out_file"
done

echo "[Step 3] Combining files into V, D, J sets..."
cat IGHV_igblast.fasta IGKV_igblast.fasta IGLV_igblast.fasta > all_V.fasta
cp IGHD_igblast.fasta all_D.fasta
cat IGHJ_igblast.fasta IGLJ_igblast.fasta > all_J.fasta

echo "[Step 4] Creating BLAST databases..."
$BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in all_V.fasta -out igblast_db/ig_v
$BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in all_D.fasta -out igblast_db/ig_d
$BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in all_J.fasta -out igblast_db/ig_j

echo "[✅ DONE] IgBLAST-ready database created in 'igblast_db/' directory."
ls -lh igblast_db
