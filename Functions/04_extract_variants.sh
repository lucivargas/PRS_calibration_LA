#!/bin/bash

# Script to find extract variant information from VCF
# Author: Luciana Vargas
# Date: 06/06/2025

# Call:
# sh ./Functions/04_extract_variants.sh --gen-data $TARGET_DATA --out-file $OUTPUT_FILE


#####################
## Define variables
#####################

# Parse named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --gen-data)
      GEN_DATA="$2"
      shift 2
      ;;
    --chr)
      CHR="$2"
      shift 2
      ;;
    --out-file)
      OUTPUT_FILE="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 --gen-data target_hg38.vcf.gz --chr 22 --out-file /path/to/output"
      exit 1
      ;;
  esac
done


#########################################
## Extract variants from VCF 
#########################################

echo " Extracting variants from $GEN_DATA..."

# Filter variants based on missingness (>90% of samples) and imputation quality (>90%) if it exists
# Then, extract CHR, POS, REF and ALT 
echo -e "CHROM\tPOS\tREF\tALT" > "$OUTPUT_FILE"
if bcftools view -h "$GEN_DATA" | grep -q '^##INFO=<ID=R2'; then
    echo "  R2 tag exists — applying imputation quality (>90%) filter..."
    bcftools view -i 'F_MISSING<0.9  && INFO/R2>0.9' "$GEN_DATA" -Ou | \
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' | awk -v chr="$CHR" '$1 == chr || $1 == "chr"chr { sub(/^chr/, "", $1); print $1 "\t" $2 "\t" $3 "\t" $4 }' >> "$OUTPUT_FILE"
else
    echo "  R2 tag does NOT exist — skipping R2 filter"
    bcftools view -i 'F_MISSING<0.9' "$GEN_DATA" -Ou | \
      bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' | awk -v chr="$CHR" '$1 == chr || $1 == "chr"chr { sub(/^chr/, "", $1); print $1 "\t" $2 "\t" $3 "\t" $4 }' >> "$OUTPUT_FILE"
fi


# Print message
GEN_COUNT=$(($(wc -l < "$OUTPUT_FILE") - 1))
echo " $GEN_COUNT variants saved to $OUTPUT_FILE"

## End
