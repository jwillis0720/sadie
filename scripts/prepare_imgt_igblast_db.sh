#!/bin/bash

# Exit on error
set -e

# Set working directory
mkdir -p imgt_refs && cd imgt_refs

BLASTPATH="../../src/sadie/reference/bin/linux"
#Define IMGT FASTA URLs
declare -A urls=(
  [IGHV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHV.fasta"
  [IGHD]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHD.fasta"
  [IGHJ]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGHJ.fasta"
  [IGKV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKV.fasta"
  [IGLV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLV.fasta"
  [IGLJ]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGLJ.fasta"
  [IGKJ]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/Homo_sapiens/IG/IGKJ.fasta"
)


echo "[Step 1] Downloading IMGT FASTA files (if not present)..."
for gene in "${!urls[@]}"; do
    filename="${gene}.fasta"
    if [ -f "$filename" ]; then
        echo "✔️  $filename already exists, skipping download."
    else
        echo "⬇️  Downloading $filename..."
        curl -s -o "$filename" "${urls[$gene]}"
        echo "✅  Downloaded $filename"
    fi
done


# Check for required tools
if ! command -v ../edit_imgt_file.pl &> /dev/null; then
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
    perl ../edit_imgt_file.pl "$in_file" > "$out_file"
done

echo "[Step 3] Combining files into V, D, J sets..."
cat IGHV_igblast.fasta IGKV_igblast.fasta IGLV_igblast.fasta > human_V.fasta
cp IGHD_igblast.fasta human_D.fasta
cat IGHJ_igblast.fasta IGLJ_igblast.fasta IGKJ_igblast.fasta > human_J.fasta

python ../convert_fasta.py human_V.fasta
python ../convert_fasta.py human_D.fasta
python ../convert_fasta.py human_J.fasta

echo "[Step 4] Creating BLAST databases..."
$BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in human_V.fasta -out igblast_db/human_V
$BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in human_D.fasta -out igblast_db/human_D
$BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in human_J.fasta -out igblast_db/human_J

cp human_V.fasta igblast_db/human_V.fasta
cp human_D.fasta igblast_db/human_D.fasta
cp human_J.fasta igblast_db/human_J.fasta

echo "[✅ DONE] IgBLAST-ready database created in 'igblast_db/' directory."
ls -lh igblast_db
