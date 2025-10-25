#!/bin/bash

# Script to filter target and reference data to keep PRS variants
# Author: Luciana Vargas
# Date: 06/24/2025

# Call:
# sh ./Functions/06_filter_data.sh --target-data $TARGET_DATA --ref-data $REF_DATA --chr $CHR \
#                                  --variant-list $OVERLAP_LIST --out $OUTPUT_FOLDER --threads $THREADS



#####################
## Define variables
#####################

# Parse named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --target-data)
      TARGET_DATA="$2"
      shift 2
      ;;
    --ref-data)
      REF_DATA="$2"
      shift 2
      ;;
    --chr)
      CHR="$2"
      shift 2
      ;;
    --pgs-id)
      PGS_ID="$2"
      shift 2
      ;;
    --variant-list)
      OVERLAP_LIST="$2"
      shift 2
      ;;
    --sample-list)
      SAMPLE_OVERLAP_LIST="$2"
      shift 2
      ;; 
    --out)
      OUTPUT_FOLDER="$2"
      shift 2
      ;;
    --threads)
      THREADS="$2"
      shift 2
      ;; 
    *)
      echo "Unknown option: $1"
      echo "Usage: $0 --target-data target_hg38.vcf.gz --ref-data ref_hg38.vcf.gz --chr 22 --pgs-id PGS000123 --variant-list overlap_PRS_target_ref_chr22.txt --sample-list overlap_sample_list_chr22.txt --out /path/to/output --threads 16"
      exit 1
      ;;
  esac
done


#############################
## Create intermediate files
#############################

# Create list of variants with chr prefix for compatibility
OVERLAP_LIST_CHR="${OVERLAP_LIST%.txt}_chr.txt"
if [[ ! -f $OVERLAP_LIST_CHR ]]; then
  awk -F'\t' 'BEGIN{OFS="\t"} {print "chr"$1, $2}' "$OVERLAP_LIST" > "$OVERLAP_LIST_CHR"
fi


####################################
## Extract PRS variants from target
####################################

# Print message
echo "Filtering target data..."

# Define output file
TARGET_TEMP="$OUTPUT_FOLDER/gen/temp_target_chr${CHR}_${PGS_ID}.vcf.gz"
TARGET_FILTERED="$OUTPUT_FOLDER/gen/target_chr${CHR}_${PGS_ID}.vcf.gz"

# Check if the target data uses chr prefix
if [ $(bcftools query -f '%CHROM\n' $TARGET_DATA | head -n 10 | grep -cv 'chr') -eq 0 ]; then
  VAR_LIST="$OVERLAP_LIST_CHR"
else
  VAR_LIST="$OVERLAP_LIST"
fi


# If the file is in dosage format, use bcftools + plink for extraction
if [[ $(bcftools view -h $TARGET_DATA | grep -c 'ID=DS') -gt 0 ]]; then

  # Filter target data, and remove any multiallelic variants (dosage extraction won't work w/ multiallelics)
  bcftools view -R "$VAR_LIST" -S "$SAMPLE_OVERLAP_LIST" "$TARGET_DATA" --threads $THREADS -Ou |  \
  bcftools norm -m -any -Oz -o "$TARGET_TEMP" --threads $THREADS 

  # Index filtered data
  bcftools index --tbi $TARGET_TEMP --threads $THREADS

  echo -e "Detected dosage information in target data. Extracting Minimac3-style dosages from filtered data..."
  
  # Extract dosage only
  plink2 --vcf $TARGET_TEMP  dosage=DS \
    --export vcf 'bgz' 'vcf-dosage=DS-only' \
    --out ${TARGET_FILTERED%.vcf.gz} --threads $THREADS
  bcftools index --tbi $TARGET_FILTERED --threads $THREADS
  
  # Remove temporary file
  rm $TARGET_TEMP
  rm $TARGET_TEMP.tbi
  
else

  # Filter target data
  bcftools view -R "$VAR_LIST" -S "$SAMPLE_OVERLAP_LIST" "$TARGET_DATA" --threads $THREADS -Oz -o "$TARGET_TEMP"
  
  # Index filtered data
  bcftools index --tbi $TARGET_TEMP --threads $THREADS

  # If not dosage, just save temporarily filtered file to final filtered file
  mv $TARGET_TEMP $TARGET_FILTERED
  mv $TARGET_TEMP.tbi $TARGET_FILTERED.tbi
  
  # Remove temporary file
  if [[ -f $TARGET_FILTERED && -f $TARGET_FILTERED.tbi ]]; then
    rm $TARGET_TEMP
    rm $TARGET_TEMP.tbi
  fi
  
fi

# Check how many variants are left in filtered file
TARGET_COUNT=$(bcftools index --nrecords $TARGET_FILTERED --threads $THREADS)
echo "After filtering based on PRS variants, $TARGET_COUNT variants remain in target data for chromosome $CHR"


#######################################
## Extract PRS variants from reference
#######################################

# Filter reference
echo "Filtering reference data..."
REF_FILTERED="$OUTPUT_FOLDER/gen/ref_chr${CHR}_${PGS_ID}.vcf.gz"

# Check if the target data uses chr prefix
if [ $(bcftools query -f '%CHROM\n' $REF_DATA | head -n 10 | grep -cv 'chr') -eq 0 ]; then
  VAR_LIST="$OVERLAP_LIST_CHR"
else
  VAR_LIST="$OVERLAP_LIST"
fi

# Filter
bcftools view -R $VAR_LIST -Oz -o $REF_FILTERED $REF_DATA --threads $THREADS
bcftools index --tbi $REF_FILTERED --threads $THREADS

# Check how many variants are left in filtered file
REF_COUNT=$(bcftools index --nrecords $REF_FILTERED --threads $THREADS)
echo "After filtering based on PRS variants, $REF_COUNT variants remain in reference data for chromosome $CHR"



###################
## Harmonize data
###################



#############################
## Clean intermediate files
#############################

rm $OVERLAP_LIST_CHR


## End
