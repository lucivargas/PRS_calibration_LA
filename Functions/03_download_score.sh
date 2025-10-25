#!/bin/bash

# Script to download PGS scoring file from the PGS Catalog
# Author: Luciana Vargas
# Date: 06/05/2025

# Call:
# sh ./Functions/03_download_score.sh --pgs-id $PGS_ID --genome-build $GENOME_BUILD --out-file $PGS_PATH


# Parse named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --pgs-id)
      PGS_ID="$2"
      shift 2
      ;;
    --genome-build)
      GENOME_BUILD="$2"
      shift 2
      ;;
    --out-file)
      OUTPUT_FILE="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 --pgs-id PGS000123 --genome-build GRCh38 --out-file /path/to/output"
      exit 1
      ;;
  esac
done


# Define the PGSCatalog API URL
API_URL="https://www.pgscatalog.org/rest/score/${PGS_ID}/"

# Query the API and extract the FTP link
FTP_URL=$(curl -s "$API_URL" | grep -oP "(?<=\"$GENOME_BUILD\":{\"positions\":\")[^\"]+")

# Check if FTP_URL is valid
if [ -z "$FTP_URL" ]; then
  echo " Error: Could not retrieve scoring file URL for $PGS_ID and genome build $GENOME_BUILD"
  exit 1
fi

# Download the scoring file
echo " Downloading scoring file from: $FTP_URL"
wget --no-verbose "$FTP_URL" -O "$OUTPUT_FILE"

echo " ✅ Download complete: ${PGS_ID}.txt.gz saved to $OUTPUT_FILE"

## End
