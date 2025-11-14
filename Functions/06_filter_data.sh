#!/bin/bash

# Script to filter target and reference data to keep PRS variants
# Author: Luciana Vargas
# Date: 06/24/2025

# Call:
# sh ./Functions/06_filter_data.sh --base-dir $BASE_DIR \
#                                  --target-data $TARGET_DATA --ref-data $REF_DATA --chr $CHR --pgs-id $PGS_ID \
#                                  --variant-list $OVERLAP_LIST --sample-list $SAMPLE_OVERLAP_LIST \
#                                  --out $OUTPUT_FOLDER --threads $THREADS


#####################
## Define variables
#####################

# Parse named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --base-dir)
      $BASE_DIR="$2"
      shift 2
      ;;
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
      echo "Usage: $0 --base-dir . --target-data target_hg38.vcf.gz --ref-data ref_hg38.vcf.gz --chr 22 --pgs-id PGS000123 --variant-list overlap_PRS_target_ref_chr22.txt --sample-list overlap_sample_list_chr22.txt --out /path/to/output --threads 16"
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


################################################################
## Create temporary file w/ PRS variants extracted from TARGET
################################################################

# Print message
echo "Filtering target data..."

# Define output files
TARGET_TEMP="$OUTPUT_FOLDER/gen/temp_target_chr${CHR}_${PGS_ID}.vcf.gz"
TARGET_FILTERED="$OUTPUT_FOLDER/gen/target_chr${CHR}_${PGS_ID}.vcf.gz"

# Check if the target data uses chr prefix
if [ $(bcftools query -f '%CHROM\n' $TARGET_DATA | head -n 10 | grep -cv 'chr') -eq 0 ]; then
  VAR_LIST="$OVERLAP_LIST_CHR"
else
  VAR_LIST="$OVERLAP_LIST"
fi


# If the file is in dosage format, use extra step w/ bcftools 
if [[ $(bcftools view -h $TARGET_DATA | grep -c 'ID=DS') -gt 0 ]]; then

  # Filter target data, and remove any multiallelic variants (dosage extraction won't work w/ multiallelics)
  bcftools view -R "$VAR_LIST" -S "$SAMPLE_OVERLAP_LIST" "$TARGET_DATA" --threads $THREADS -Ou |  \
    bcftools norm -m -any -Oz -o "$TARGET_TEMP" --threads $THREADS 
  bcftools index -f --tbi "$TARGET_TEMP" --threads $THREADS
  
else

  # Filter target data
  bcftools view -R "$VAR_LIST" -S "$SAMPLE_OVERLAP_LIST" "$TARGET_DATA" \
                --threads $THREADS -Oz -o "$TARGET_TEMP"
  bcftools index -f --tbi "$TARGET_TEMP" --threads $THREADS
  
fi


#######################################
## Extract PRS variants from reference
#######################################

# Print message
echo "Filtering reference data..."

# Define output files
REF_TEMP="$OUTPUT_FOLDER/gen/temp_ref_chr${CHR}_${PGS_ID}.vcf.gz"
REF_FILTERED="$OUTPUT_FOLDER/gen/ref_chr${CHR}_${PGS_ID}.vcf.gz"

# Check if the target data uses chr prefix & define appropriate variant list
if [ $(bcftools query -f '%CHROM\n' $REF_DATA | head -n 10 | grep -cv 'chr') -eq 0 ]; then
  VAR_LIST="$OVERLAP_LIST_CHR"
else
  VAR_LIST="$OVERLAP_LIST"
fi

# Filter variants
bcftools view -R $VAR_LIST -Oz -o $REF_TEMP $REF_DATA --threads $THREADS
bcftools index -f --tbi $REF_TEMP --threads $THREADS


###################################################
## Find batch effect between TARGET and REFERENCE 
###################################################

# Print message
echo " Finding any batch effects between target and reference data after adjusting for PC1 (p < 5e-08)"

# Merge TARGET and REFERENCE data
MERGED_VCF="$OUTPUT_FOLDER/gen/gwas/merged_chr${CHR}_${PGS_ID}.vcf.gz"
bcftools merge $TARGET_TEMP $REF_TEMP -Oz -o $MERGED_VCF --threads $THREADS
bcftools index -f --tbi $MERGED_VCF --threads $THREADS

# Convert merged data to plink binaries & Calculate PCA
MERGED_BFILE="$OUTPUT_FOLDER/gen/gwas/merged_chr${CHR}_${PGS_ID}"
plink2 --vcf $MERGED_VCF \
       --pca 10 \
       --make-pgen --out $MERGED_BFILE --threads $THREADS

# Create a phenotype file for the batch effect (captures REF vs. TARGET)
PHENO="$OUTPUT_FOLDER/gen/gwas/pheno_chr${CHR}_${PGS_ID}.txt"
echo "IID batch" > $PHENO
bcftools query -l $TARGET_TEMP | awk '{print $1, 1}' >> $PHENO
bcftools query -l $REF_TEMP | awk '{print $1, 2}' >> $PHENO

# Create a covariate file for PC1 to use as a covariate (captures pop. structure)
COVAR="$OUTPUT_FOLDER/gen/gwas/covar_chr${CHR}_${PGS_ID}.txt"
awk '{print $1, $2}' $MERGED_BFILE.eigenvec > $COVAR

# Run a GWAS on the merged dataset to find variants with batch effect, adjusting for pop stratification (PC1)
HARM_BFILE="$OUTPUT_FOLDER/gen/gwas/postgwas_chr${CHR}_${PGS_ID}"
plink2 --pfile $MERGED_BFILE \
       --pheno $PHENO --pheno-name batch \
       --covar $COVAR --covar-name PC1 \
       --glm --out $HARM_BFILE --threads $THREADS

# Calculate allele frequencies
plink2 --pfile $MERGED_BFILE \
       --pheno $PHENO --pheno-name batch \
       --keep-if 'batch==1' --freq --out "${HARM_BFILE}_target"
plink2 --pfile $MERGED_BFILE \
       --pheno $PHENO --pheno-name batch \
       --keep-if 'batch==2' --freq --out  "${HARM_BFILE}_reference"

# Extract batch-associated variants from GWAS result & update snp lists
OVERLAP_LIST="$OUTPUT_FOLDER/snp_lists/overlap_${PGS_ID}_chr${CHR}_variant_list.txt"
EXCLUDE_LIST="$OUTPUT_FOLDER/snp_lists/excluded_${PGS_ID}_chr${CHR}_variant_list.txt"
Rscript "$BASE_DIR/Functions/06_01_process_gwas.R" \
        --gwas-out $HARM_BFILE --pheno-file $PHENO --pca-file "$MERGED_BFILE.eigenvec" \
        --overlap-list $OVERLAP_LIST --exclude-list $EXCLUDE_LIST 

# Path to list of variants w/ significant batch effect from GWAS
BATCH_VARS="$OUTPUT_FOLDER/gen/gwas/postgwas_chr${CHR}_${PGS_ID}_batch_variants.txt"

# Create list of batch variants with chr prefix for compatibility
BATCH_VARS_CHR="${BATCH_VARS%.txt}_chr.txt"
awk -F'\t' 'BEGIN{OFS="\t"} {print "chr"$1, $2}' "$BATCH_VARS" > "$BATCH_VARS_CHR"


######################################
## Remove batch effect from REFERENCE 
######################################

# Check if the reference data uses chr prefix & define appropriate variant list
if [ $(bcftools query -f '%CHROM\n' $REF_TEMP | head -n 10 | grep -cv 'chr') -eq 0 ]; then
  VAR_LIST="$BATCH_VARS_CHR"
else
  VAR_LIST="$BATCH_VARS"
fi

# Filter data from REFERENCE (skip if no batch effect found)
if [[ $(wc -l < "$VAR_LIST") -gt 0 ]];then
  bcftools view -T ^"$VAR_LIST" "$REF_TEMP" -Oz -o "$REF_FILTERED" --threads $THREADS 
  bcftools index -f --tbi "$REF_FILTERED" --threads $THREADS
else
  cp "$REF_TEMP" "$REF_FILTERED"
fi

# Check how many variants are left in REFERENCE filtered file
REF_COUNT=$(bcftools index --nrecords $REF_FILTERED --threads $THREADS)
echo "After filtering based on PRS variants, $REF_COUNT variants remain in reference data for chromosome $CHR"


######################################
## Remove batch effect from TARGET 
######################################

# Check if the target data uses chr prefix & define appropriate variant list
if [ $(bcftools query -f '%CHROM\n' $TARGET_TEMP | head -n 10 | grep -cv 'chr') -eq 0 ]; then
  VAR_LIST="$BATCH_VARS_CHR"
else
  VAR_LIST="$BATCH_VARS"
fi

# If in dosage format, format to keep dosage only & remove SNPs w/ batch effect
if [[ $(bcftools view -h $TARGET_DATA | grep -c 'ID=DS') -gt 0 ]]; then
  
  echo -e "Detected dosage information in target data. Extracting Minimac3-style dosages from filtered data..."
  
  # Create BED-like file with variant list
  BATCH_VARS_PLINK="${VAR_LIST%.txt}.bed"
  awk '{print $1, $2, $2}' $VAR_LIST > $BATCH_VARS_PLINK
  
  # Filter data from TARGET & Extract dosage only
  plink2 --vcf $TARGET_TEMP  dosage=DS \
    --exclude range $BATCH_VARS_PLINK \
    --export vcf 'bgz' 'vcf-dosage=DS-only' \
    --out ${TARGET_FILTERED%.vcf.gz} --threads $THREADS
  bcftools index -f --tbi $TARGET_FILTERED --threads $THREADS

else
  
  # Filter data from TARGET (skip if no batch effect found)
  if [[ $(wc -l < "$VAR_LIST") -gt 0 ]];then
    bcftools view -T ^"$VAR_LIST" "$TARGET_TEMP" -Oz -o "$TARGET_FILTERED" --threads $THREADS 
    bcftools index -f --tbi "$TARGET_FILTERED" --threads $THREADS
  else
    cp "$TARGET_TEMP" "$TARGET_FILTERED"
  fi
  
fi

# Check how many variants are left in TARGET filtered file
TARGET_COUNT=$(bcftools index --nrecords $TARGET_FILTERED --threads $THREADS)
echo "After filtering based on PRS variants, $TARGET_COUNT variants remain in target data for chromosome $CHR"


########################################
## Clean intermediate & temporary files
########################################

# # Remove intermediate files
# rm -f $BATCH_VARS_PLINK
# rm $OVERLAP_LIST_CHR
# rm $BATCH_VARS_CHR
# 
# # If pipeline was successful, remove temporary files
# if [[ -f $TARGET_FILTERED && -f $TARGET_FILTERED.tbi && -f $REF_FILTERED && -f $REF_FILTERED.tbi ]]; then
#   rm $TARGET_TEMP*
#   rm $REF_TEMP*
#   rm $MERGED_VCF
#   rm $MERGED_BFILE.bed
#   rm $MERGED_BFILE.bim
#   rm $MERGED_BFILE.fam
#   rm $PHENO
#   rm $COVAR
# fi


## End
