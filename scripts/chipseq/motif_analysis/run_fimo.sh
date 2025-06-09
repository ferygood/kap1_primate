#!/bin/bash

# Usage:
# ./run_fimo.sh --input-folder <meme_folder> --oc <fimo_output_root> --fasta <te_fasta_file>

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input-folder) MEME_FOLDER="$2"; shift ;;
        --oc) OUTPUT_ROOT="$2"; shift ;;
        --fasta) FASTA="$2"; shift ;;
        *) echo "❌ Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

# Check for required arguments
if [[ -z "$MEME_FOLDER" || -z "$OUTPUT_ROOT" || -z "$FASTA" ]]; then
    echo "❗ Usage: $0 --input-folder <meme_folder> --oc <fimo_output_root> --fasta <te_fasta_file>"
    exit 1
fi

# Create output root directory if it doesn't exist
mkdir -p "$OUTPUT_ROOT"

# Loop through all .meme files and run FIMO
for MEME_FILE in "$MEME_FOLDER"/*.meme; do
    BASENAME=$(basename "$MEME_FILE" .meme)
    OUTDIR="$OUTPUT_ROOT/$BASENAME"
    mkdir -p "$OUTDIR"
    echo "🔍 Running FIMO on $MEME_FILE using $FASTA → $OUTDIR"
    fimo --oc "$OUTDIR" "$MEME_FILE" "$FASTA"
done

echo "✅ All FIMO jobs completed."
