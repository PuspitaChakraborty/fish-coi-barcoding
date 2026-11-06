#!/usr/bin/env bash

# Project: Fish COI DNA Barcoding
# Purpose: Download annotated GenBank records from NCBI

set -Eeuo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
INPUT="$PROJECT_ROOT/data/accessions.csv"
OUTPUT="$PROJECT_ROOT/data/raw"
LOG="$PROJECT_ROOT/logs/genbank_download.log"

mkdir -p "$OUTPUT" "$PROJECT_ROOT/logs"

echo "Starting NCBI GenBank record retrieval..."

while IFS=, read -r target_species genbank_organism accession expected_length notes
do
    # Skip CSV header
    [[ "$target_species" == "target_species" ]] && continue

    # Validate accession format
    if [[ ! "$accession" =~ ^[A-Z]{1,4}[0-9]{5,9}(\.[0-9]+)?$ ]]; then
        echo "ERROR: Invalid accession: $accession"
        exit 1
    fi

    FILE="$OUTPUT/${accession}.gb"

    if [[ -s "$FILE" ]]; then
        echo "Already exists, skipping: $accession"
        continue
    fi

    echo "Downloading: $target_species ($accession)"

    if curl --fail --silent --show-error --location \
        --retry 3 --retry-delay 2 \
        --get "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
        --data-urlencode "db=nuccore" \
        --data-urlencode "id=$accession" \
        --data-urlencode "rettype=gb" \
        --data-urlencode "retmode=text" \
        -o "$FILE" &&
        grep -q "^LOCUS" "$FILE" &&
        grep -q "^ACCESSION" "$FILE"; then

        echo "SUCCESS: $accession" | tee -a "$LOG"
    else
        rm -f "$FILE"
        echo "FAILED: $accession" | tee -a "$LOG"
        exit 1
    fi

    sleep 0.4

done < "$INPUT"

echo "NCBI retrieval completed."
echo "Records stored in: $OUTPUT"
