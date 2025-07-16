#!/bin/bash

# Exit on error
set -e

# Set working directory
gcwd=$(dirname "$0")
echo "Current working directory: $gcwd"

BLASTPATH="${gcwd}/../src/sadie/reference/bin/linux"
ABS_BLASTPATH="$(cd "$(dirname "$BLASTPATH")" && pwd)/$(basename "$BLASTPATH")"
echo "$ABS_BLASTPATH"

cd ${gcwd}

#Define IMGT FASTA URLs
organism="Macaca mulatta"  # Change this to the desired organism,
# dictionary to change the scientific name to mouse
# e.g., "Mus musculus" for mouse

declare -A species_dict=(
  ["Canis lupus familiaris"]="dog"
  ["Homo sapiens"]="human"
  ["Macaca mulatta"]="macaque"
  ["Mus musculus"]="mouse"
  ["Oryctolagus cuniculus"]="rabbit"
  ["Rattus norvegicus"]="rat"
)

# Current working directory
gcwd=$(pwd)

organism_link=$(echo $organism | sed 's/ /_/g')

mkdir -p ${species_dict[$organism]} && cd ${species_dict[$organism]}

declare -A urls=(
  [IGHV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/${organism_link}/IG/IGHV.fasta"
  [IGHD]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/${organism_link}/IG/IGHD.fasta"
  [IGHJ]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/${organism_link}/IG/IGHJ.fasta"
  [IGKV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/${organism_link}/IG/IGKV.fasta"
  [IGLV]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/${organism_link}/IG/IGLV.fasta"
  [IGLJ]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/${organism_link}/IG/IGLJ.fasta"
  [IGKJ]="https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/${organism_link}/IG/IGKJ.fasta"
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
if ! command -v $gcwd/edit_imgt_file.pl &> /dev/null; then
    echo "[ERROR] edit_imgt_file.pl not found. Please add IgBLAST's bin directory to your PATH."
    exit 1
fi

if ! command -v $ABS_BLASTPATH/makeblastdb &> /dev/null; then
    echo "[ERROR] makeblastdb not found. Please install BLAST+ tools and add to PATH."
    exit 1
fi

echo "[Step 2] Cleaning headers with edit_imgt_file.pl..."
for gene in "${!urls[@]}"; do
    in_file="${gene}.fasta"
    out_file="${gene}_igblast.fasta"
    echo "Processing $in_file -> $out_file"
    perl $gcwd/edit_imgt_file.pl "$in_file" > "$out_file"
done

echo "[Step 3] Combining files into V, D, J sets..."
cat IGHV_igblast.fasta IGKV_igblast.fasta IGLV_igblast.fasta > ${species_dict[$organism]}_V.fasta
cp IGHD_igblast.fasta ${species_dict[$organism]}_D.fasta
cat IGHJ_igblast.fasta IGLJ_igblast.fasta IGKJ_igblast.fasta > ${species_dict[$organism]}_J.fasta

python $gcwd/convert_fasta.py ${species_dict[$organism]}_V.fasta
python $gcwd/convert_fasta.py ${species_dict[$organism]}_D.fasta
python $gcwd/convert_fasta.py ${species_dict[$organism]}_J.fasta

echo "[Step 4] Creating BLAST databases..."
$ABS_BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in ${species_dict[$organism]}_V.fasta -out igblast_db/${species_dict[$organism]}_V
$ABS_BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in ${species_dict[$organism]}_D.fasta -out igblast_db/${species_dict[$organism]}_D
$ABS_BLASTPATH/makeblastdb -dbtype nucl -parse_seqids -in ${species_dict[$organism]}_J.fasta -out igblast_db/${species_dict[$organism]}_J

cp ${species_dict[$organism]}_V.fasta igblast_db/${species_dict[$organism]}_V.fasta
cp ${species_dict[$organism]}_D.fasta igblast_db/${species_dict[$organism]}_D.fasta
cp ${species_dict[$organism]}_J.fasta igblast_db/${species_dict[$organism]}_J.fasta

echo "[✅ DONE] IgBLAST-ready database created in 'igblast_db/' directory."
ls -lh igblast_db

cp igblast_db/* ${gcwd}/../src/sadie/airr/data/germlines/Ig/blastdb/${species_dict[$organism]}
cd "$gcwd"
echo "Current working directory: $(pwd)"
echo "Script completed successfully."