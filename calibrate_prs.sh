#!/bin/bash

# Example code on how to run the pipeline
# sh calibrate_prs.sh \
#	    --program-path . \
# 		--target-data /path/to/target \
#  		--target-lai /path/to/lai \
# 		--ref-data /path/to/ref \
# 		--ref-map /path/to/refmap \
#     --chr 22 \
# 		--pgs-id PGS000123 \
#     --genome-build GRCh38 \
# 		--genome-dist genetic_map_hg38_withX.txt.gz \
# 		--out /path/to/save/results \
#     --threads 16


#######################
## Prepare environment
#######################


# Parse named arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --program-path)
      BASE_DIR="$2"
      shift 2
      ;;
    --target-data)
      TARGET_DATA="$2"
      shift 2
      ;;
    --target-lai)
      TARGET_LAI="$2"
      shift 2
      ;;
    --ref-data)
      REF_DATA="$2"
      shift 2
      ;;
    --ref-map)
      REF_MAP="$2"
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
    --genome-build)
      GENOME_BUILD="$2"
      shift 2
      ;;
    --genome-dist)
      GENOME_DIST="$2"
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
    --resume)
      RESUME="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: $0  --program-path . 
      --target-data /path/to/target 
      --target-lai /path/to/lai
      --ref-data /path/to/ref
      --ref-map /path/to/map
      --chr 22
      --pgs-id PGS000123
      --genome-build GRCh38
      --genome-dist /path/to/gendist
      --out /path/to/output
      --threads 16
      --resume true"
      exit 1
      ;;
  esac
done


# If not specified, get base dir 
if [[ -z "$BASE_DIR+x" ]]; then
  BASE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
fi


# If not specified, set number of cores to use
if [[ -z "${THREADS+x}" ]]; then
  THREADS=$(yq '.threads' $BASE_DIR/config.yml)
fi


# Load values from config.yml
NSIMULATIONS=$(yq '.nsimulations' $BASE_DIR/config.yml)
NBOOTSTRAPS=$(yq '.nbootstraps' $BASE_DIR/config.yml)
MIN_OVERLAP=$(yq '.min_overlap' $BASE_DIR/config.yml)
RANDSEED=$(yq '.random_seed' $BASE_DIR/config.yml)


# Check that simulation and bootstrap parameters meet minimum requiment
SIMSUM=$(($NSIMULATIONS + $NBOOTSTRAPS))
if [[ $NSIMULATIONS -lt 50 ]]; then
  echo "The number of simulations does not meet the minimum requiment of 50, please change config file." | tee -a $MAIN_LOG
  exit 1
elif [[ $SIMSUM -lt 100 ]]; then
  echo "The number of simulations and bootstraps does sum to 100, the minimum requiment for percentile calculation, please change config file." |  tee -a $MAIN_LOG
  exit 1
fi  


## Create directories to save output
mkdir -p $OUTPUT_FOLDER/snp_lists
mkdir -p $OUTPUT_FOLDER/gen
mkdir -p $OUTPUT_FOLDER/gen/gwas
mkdir -p $OUTPUT_FOLDER/gen/lai
mkdir -p $OUTPUT_FOLDER/pgs
mkdir -p $OUTPUT_FOLDER/pgs/simulations
mkdir -p $OUTPUT_FOLDER/logs


#######################
## Configure logs
#######################

# Start log files
TIMESTAMP_BEGIN=$(date "+%Y%m%d_%H%M")
MAIN_LOG="$OUTPUT_FOLDER/logs/calibrate_chr${CHR}_$TIMESTAMP_BEGIN.log"
echo "[$(date '+%Y-%m-%d %H:%M:%S')]" > $MAIN_LOG

# Define function for logging
run_logged() {
  
  local LABEL=$1
  shift
  
  # Define named pipe objects in environment
  local OUT_PIPE
  OUT_PIPE=$(mktemp -u)
  mkfifo "$OUT_PIPE" 
  
  # Redirect pipes to logs 
  tee -a "$MAIN_LOG" < "$OUT_PIPE" &
  
  # Initiate exit code capturing
  local status
  
  {
    # Print start message 
    echo -e "\n[$(date '+%Y-%m-%d %H:%M:%S')] [$LABEL] Starting..."
    
    # Print messages from script, or exit if failed
    "$@" || {
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] [$LABEL] FAILED: exit code $?"
        exit 1   # kill the whole script
    }
    
    # Store the exit code of the inner command
    status=$?
    
    # Print end message 
    if [[ $status -eq 0 ]]; then
      echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] [$LABEL] Finished.\n"
    else 
      echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] [$LABEL] FAILED: exit code $status \n"
    fi
    
    # Write exit status to a temporary file
    echo "$status" > "$OUT_PIPE.status"
    
  # Save messages to log files
  }  > "$OUT_PIPE" 2>&1
  
  wait
  
  # Read status from file, clean up file, and return status to console
  if [[ -f "$OUT_PIPE.status" ]]; then
    status=$(<"$OUT_PIPE.status")
    rm -f "$OUT_PIPE.status"
  else
    status=1
  fi
  
  # Clean environment
  rm "$OUT_PIPE"
  
  # Return status
  return $status
  
}


###################
## Initialize logs
###################


# Make sure RESUME flag is correctly specified
RESUME=$(echo "$RESUME" | tr '[:upper:]' '[:lower:]')
if [[ "$RESUME" =~ ^(t|true|yes|1)$ ]]; then
  RESUME="true"
elif [[ "$RESUME" =~ ^(f|false|no|0)$ ]]; then
  RESUME="false"
else
  echo "Resume flag not recognized, please specify if \"true\" or \"false\" and rerun the script. Setting it to \"false\"." | tee -a $MAIN_LOG
  RESUME="false"
fi


# Log the arguments passed to function 
echo "Options in effect: 
      --program-path $BASE_DIR
      --target-data $TARGET_DATA
      --target-lai $TARGET_LAI 
      --ref-data $REF_DATA
      --ref-map $REF_MAP
      --chr $CHR
      --pgs-id $PGS_ID
      --genome-build $GENOME_BUILD
      --genome-dist $GENOME_DIST
      --out $OUTPUT_FOLDER
      --threads $THREADS
      --resume $RESUME" | tee -a $MAIN_LOG


## Check that number of threads set match those available in system
AVAILABLE_THREADS=$(nproc)
if [[ $THREADS -gt $AVAILABLE_THREADS ]]; then
 echo "The number of threads is set to $THREADS, but there are only $AVAILABLE_THREADS available. Changing threads argument to $AVAILABLE_THREADS." |  tee -a $MAIN_LOG
 THREADS=$AVAILABLE_THREADS
fi


# Log number of available cores
echo -e "Using $THREADS available cores for computation\n" | tee -a $MAIN_LOG


## Check that input files exist
echo -e "Checking input files..." | tee -a $MAIN_LOG 
for file in "$TARGET_DATA" "$TARGET_LAI" "$REF_DATA" "$REF_MAP" "$GENOME_DIST"; do
  if [[ ! -f "$file" ]]; then
      echo " ERROR: File not found: $file" | tee -a $MAIN_LOG 
      exit 1
  fi
done
echo -e "All input files exist.\n" | tee -a $MAIN_LOG 
echo -e "\nStarting main script...\n"  | tee -a $MAIN_LOG

#####################################################
## STEP 1: Check if required packages are installed
#####################################################

run_logged "1_CHECK_PACKAGES" sh "$BASE_DIR/Functions/01_check_requirements.sh" 
if [ $? -ne 0 ]; then
    echo " ERROR: exiting main script due to missing dependencies" | tee -a $MAIN_LOG 
    exit 1
fi


#############################
## STEP 2: Check input data
#############################


# Check if reference VCF is phased 
IS_PHASED=$(bcftools view -H "$REF_DATA" | head -n 100 | awk '{if ($0 ~ /\|/) count++} END {print count}')
if [ "$IS_PHASED" -eq 0 ]; then
    run_logged "2_CHECK_VCF" echo -e " ERROR: Reference VCF is not phased. Please make sure data is phased before proceeding with the calibration. \n ❌ Exiting main script due to reference data not being phased" 
    exit 1
else
    run_logged "2_CHECK_VCF" echo -e  " Reference VCF is phased."
fi


# Check if samples match between LAI and target VCF
run_logged "2_CHECK_LAI" Rscript "$BASE_DIR/Functions/02_check_target.R" \
                         --chr $CHR --target-data $TARGET_DATA --target-lai $TARGET_LAI --out $OUTPUT_FOLDER || {
    echo " ERROR: 2_CHECK_LAI failed at line $LINENO" | tee -a "$MAIN_LOG"
    exit 1
    }

# This script has saved the LAI codes for the data, define location of file:
LAI_CODES="$OUTPUT_FOLDER/gen/$(basename $TARGET_LAI).codes"

# This script has saved a list of overlapping samples between VCF and LAI, define location of file:
SAMPLE_OVERLAP_LIST="$OUTPUT_FOLDER/snp_lists/overlap_sample_list_chr${CHR}.tsv"


###########################################
## STEP 3: Download score from PGSCatalog
###########################################

PGS_PATH="$OUTPUT_FOLDER/pgs/${PGS_ID}_hmPOS_${GENOME_BUILD}.txt.gz"
if [[ ! -f "$PGS_PATH" ]]; then
  run_logged "3_DOWNLOAD_PGS" sh "$BASE_DIR/Functions/03_download_score.sh" \
                                --pgs-id $PGS_ID --genome-build $GENOME_BUILD --out-file $PGS_PATH || {
    echo " ERROR: 3_DOWNLOAD_PGS failed at line $LINENO" | tee -a "$MAIN_LOG"
    exit 1
    }
else
  run_logged "3_DOWNLOAD_PGS" echo " PGS scoring file already exists, skipping download"
fi

# Make sure this step produced the file needed for next steps
if [[ ! -f "$PGS_PATH" ]]; then
  exit 1
fi
  

#########################################################
## STEP 4: Extract list of variants present in VCF files
#########################################################

# Extract variant list from TARGET data  
TARGET_LIST="$OUTPUT_FOLDER/snp_lists/target_chr${CHR}_variant_list.txt"
(
  if [[ ! -f "$TARGET_LIST" || "$RESUME" == "false" ]]; then
      run_logged "4_GET_VARS_TARGET" sh "$BASE_DIR/Functions/04_extract_variants.sh" \
                                      --gen-data $TARGET_DATA --chr $CHR --out-file $TARGET_LIST || {
      echo " ERROR: 4_GET_VARS_TARGET failed at line $LINENO" | tee -a "$MAIN_LOG"
      exit 1
  }
  else
      run_logged "4_GET_VARS_TARGET" echo " List of variants in target data for chr$CHR already exists, skipping extraction"
  fi
    ) &
    
# Extract variant list from REFERENCE data  
REF_LIST="$OUTPUT_FOLDER/snp_lists/ref_chr${CHR}_variant_list.txt"
(
  if [[ ! -f "$REF_LIST" || "$RESUME" == "false" ]]; then
      run_logged "4_GET_VARS_REF" sh "$BASE_DIR/Functions/04_extract_variants.sh" \
                                      --gen-data $REF_DATA --chr $CHR --out-file $REF_LIST || {
    echo " ERROR: 4_GET_VARS_REF failed at line $LINENO" | tee -a "$MAIN_LOG"
    exit 1
  }
  else
      run_logged "4_GET_VARS_REF" echo " List of variants in reference data for chr$CHR already exists, skipping extraction"
  fi
    ) &

# Wait for both jobs to finish
wait

# Make sure this step produced the files needed for next steps
if [[ ! -f "$TARGET_LIST" || ! -f "$REF_LIST" ]]; then
  exit 1
fi

############################################################
## STEP 5: Identify variants present in PGS, TARGET and REF
############################################################

OVERLAP_LIST="$OUTPUT_FOLDER/snp_lists/overlap_${PGS_ID}_chr${CHR}_variant_list.txt"

if [[ ! -f "$OVERLAP_LIST"  || "$RESUME" == "false" ]]; then
    run_logged "5_GET_OVERLAP_VARS" Rscript "$BASE_DIR/Functions/05_find_overlapping_variants.R" \
                                    --target-list $TARGET_LIST --ref-list $REF_LIST \
                                    --pgs-path $PGS_PATH --chr $CHR --out $OUTPUT_FOLDER \
                                    --program-path $BASE_DIR || {
    echo " ERROR: 5_GET_OVERLAP_VARS failed at line $LINENO " | tee -a "$MAIN_LOG"
    exit 1
    }
else
    run_logged "5_GET_OVERLAP_VARS" echo " List of overlapping variants between PGS, target and reference data already exists, skipping. "
fi

# Make sure this step produced the file needed for next steps
if [[ ! -f "$OVERLAP_LIST" ]]; then
  exit 1
fi

######################################################
## STEP 6: Filter TARGET and REF to keep PGS variants
######################################################

# Define output files
TARGET_FILTERED="$OUTPUT_FOLDER/gen/target_chr${CHR}_${PGS_ID}.vcf.gz"
REF_FILTERED="$OUTPUT_FOLDER/gen/ref_chr${CHR}_${PGS_ID}.vcf.gz"

# If any of the output files do NOT exist, or if set to NOT resume pipeline, then run code
if [[ ! -f "$TARGET_FILTERED"  || ! -f "$REF_FILTERED" ||  "$RESUME" == "false" ]]; then
    run_logged "6_FILTER_VARS" sh "$BASE_DIR/Functions/06_filter_data.sh" \
                               --base-dir $BASE_DIR \
                               --target-data $TARGET_DATA --ref-data $REF_DATA \
                               --chr $CHR --pgs-id $PGS_ID \
                               --variant-list $OVERLAP_LIST --sample-list $SAMPLE_OVERLAP_LIST \
                               --out $OUTPUT_FOLDER --threads $THREADS || {
    echo " ERROR: 6_FILTER_VARS failed at line $LINENO " | tee -a "$MAIN_LOG"
    exit 1
    }
else
    run_logged "6_FILTER_VARS" echo " Filtered target and reference genetic files already exist, skipping filtering. "
fi

# Make sure this step produced the files needed for next steps
if [[ ! -f "$TARGET_FILTERED" || ! -f "$REF_FILTERED" ]]; then
  exit 1
fi


########################################
## STEP 7: Apply PGS to reference data 
########################################

# Define output file
PGS_REF="$OUTPUT_FOLDER/pgs/${PGS_ID}_ref_chr$CHR.tsv" 

# Get number of samples in reference data
N_REF=$(bcftools query -l $REF_FILTERED | wc -l)

# If there is a previous scored file, check if file is complete, or remove it and rerun step
if [[ -f "$PGS_REF" ]]; then
  N_REF_SCORE=$(wc -l < "$PGS_REF" | awk '{print int($1/2)}')
  if [[ $N_REF_SCORE -lt $N_REF ]]; then
    echo -e " Incomplete reference scoring file found, removing $PGS_REF \n" | tee -a "$MAIN_LOG"
    rm $PGS_REF
  fi
fi

# If file does NOT exist, or if set to NOT resume pipeline, then run code
if [[ ! -f "$PGS_REF"  || "$RESUME" == "false" ]]; then
    run_logged "7_SCORE_REF" Rscript "$BASE_DIR/Functions/07_apply_pgs_ref.R" \
                            --ref-data $REF_FILTERED --ref-map $REF_MAP --chr $CHR \
                            --pgs-id $PGS_ID --pgs-path $PGS_PATH --variant-list $OVERLAP_LIST \
                            --out-file $PGS_REF --threads $THREADS --program-path $BASE_DIR || {
    echo " ERROR: 7_SCORE_REF failed at line $LINENO" | tee -a $MAIN_LOG
    exit 1
    }
else
    run_logged "7_SCORE_REF" echo " Scored reference genetic data already exist, skipping scoring"
fi

# Make sure this step produced the file needed for next steps
if [[ ! -f $PGS_REF ]]; then
  exit 1
fi

####################################
## STEP 8: Apply PGS to target data 
####################################

# Define output file
PGS_TARGET="$OUTPUT_FOLDER/pgs/${PGS_ID}_target_chr$CHR.tsv"

# Get number of samples in target data
N_TARGET=$(bcftools query -l $TARGET_FILTERED | wc -l)

# If there is a previous scored file, check if file is complete, or remove it and rerun step
if [[ -f "$PGS_TARGET" ]]; then
  N_TARGET_SCORE=$(wc -l < "$PGS_TARGET")
  if [[ $N_TARGET_SCORE -lt $N_TARGET ]]; then
    echo -e " Incomplete target scoring file found, removing $PGS_TARGET \n" | tee -a "$MAIN_LOG"
    rm $PGS_TARGET
  fi
fi

# If file does NOT exist, or if set to NOT resume pipeline, then run code
if [[ ! -f "$PGS_TARGET"  || "$RESUME" == "false" ]]; then
  run_logged "8_SCORE_TARGET" Rscript "$BASE_DIR/Functions/08_apply_pgs_target.R" \
                              --gen-data $TARGET_FILTERED --chr $CHR \
                              --pgs-id $PGS_ID --pgs-path $PGS_PATH --variant-list $OVERLAP_LIST \
                              --out-file $PGS_TARGET --out $OUTPUT_FOLDER --threads $THREADS \
                              --program-path $BASE_DIR || {
    echo " ERROR: 8_SCORE_TARGET failed at line $LINENO" | tee -a $MAIN_LOG
    exit 1
    }
else
    run_logged "8_SCORE_TARGET" echo " Scored target genetic data already exist, skipping scoring"
fi

# Make sure this step produced the file needed for next steps
if [[ ! -f $PGS_TARGET ]]; then
  exit 1
fi

########################
# Start PRS calibration
########################

# We will run steps 9 and 10 for each sample in our target data. 
# To make the process faster, we will use the GNU Parallel tool

# Extract sample names
SAMPLES=($(bcftools query -l "$TARGET_FILTERED"))

                
##################################
## STEP 9: Match LAI and PRS SNPs 
##################################

# Description: LAI can be inferred up to array data, but imputed data is often needed for PRS calculation. 
# To have both LAI and genotype data for up to as many variants in the PRS as possible, 
# we will add a step to expand LAI to all variants present in the PRS and genetic data.
  
# Define function to match LAI
match_lai() {
    
    # Define sample
    local SAMPLE=$1
    
    # Define output file
    LAI_SAMPLE="$OUTPUT_FOLDER/gen/lai/lai_${PGS_ID}_${SAMPLE}_chr${CHR}.tsv"
    
    # If file does NOT exist, or if set to NOT resume pipeline, then run code
    if [[ ! -f $LAI_SAMPLE  || "$RESUME" == "false" ]]; then
        Rscript "$BASE_DIR/Functions/09_match_lai_to_prs.R" \
                --sample $SAMPLE --target-lai $TARGET_LAI --chr $CHR \
                --genomeDist $GENOME_DIST \
                --pgs-id $PGS_ID --variant-list $OVERLAP_LIST \
                --out $OUTPUT_FOLDER
    else
        echo " Expanded LAI for sample $SAMPLE already exists, skipping" 
    fi
    
    # Make sure this step produced the file needed for next steps
    if [[ ! -f $LAI_SAMPLE ]]; then
        echo "ERROR: Missing $LAI_SAMPLE" | tee -a "$MAIN_LOG"
        return 1  
    fi
    
} 

# Export relevant variables and functions to parallel
export BASE_DIR SAMPLES TARGET_LAI LAI_CODES CHR \
       GENOME_DIST PGS_ID OVERLAP_LIST \
       OUTPUT_FOLDER RESUME MAIN_LOG
export -f match_lai 

# Run STEP 9 for all target samples in parallel
run_logged "9_MATCH_LAI"  parallel -j $THREADS match_lai ::: "${SAMPLES[@]:9000:9200}" || {
    echo "ERROR: 9_MATCH_LAI parallel step failed" 
    exit 1
}


##################################################
## STEP 10: Create simulated reference population
##################################################

run_simulation(){
    
    # Define sample
    local SAMPLE=$1
  
    # Define input file
    LAI_SAMPLE="$OUTPUT_FOLDER/gen/lai/lai_${PGS_ID}_${SAMPLE}_chr${CHR}.tsv"
    
    # Define output file
    SIM_SAMPLE="$OUTPUT_FOLDER/pgs/simulations/simulation_${PGS_ID}_${SAMPLE}_chr${CHR}.tsv"
  
    # If file does NOT exist, or if set to NOT resume pipeline, then run code
    if [[ ! -f "$SIM_SAMPLE"  || "$RESUME" == "false" ]]; then
        Rscript $BASE_DIR/Functions/10_simulations.R \
                --chr $CHR --sample $SAMPLE --lai-sample $LAI_SAMPLE --lai-codes $LAI_CODES \
                --pgs-id $PGS_ID --pgs-target $PGS_TARGET --pgs-ref $PGS_REF --ref-map $REF_MAP \
                --nsimulations $NSIMULATIONS --nbootstraps $NBOOTSTRAPS --seed $RANDSEED \
                --out $OUTPUT_FOLDER
    else
        echo " Simulations for sample $SAMPLE already exist, skipping"
    fi


    # Make sure this step produced the file needed for next steps
    if [[ ! -f $SIM_SAMPLE ]]; then
        echo "ERROR: Missing $SIM_SAMPLE" | tee -a "$MAIN_LOG"
        return 1  
    fi
}
    
# Export relevant variables and functions to parallel
export BASE_DIR SAMPLES LAI_CODES CHR \
       PGS_ID PGS_REF PGS_TARGET REF_MAP \
       NSIMULATIONS NBOOTSTRAPS RANDSEED \
       OUTPUT_FOLDER RESUME MAIN_LOG
export -f run_simulation 

# Run STEP 10 for first 200 target samples in parallel
run_logged "10_SIMULATION"  parallel --ungroup  -j $THREADS run_simulation ::: "${SAMPLES[@]:9000:9200}" || {
    echo "ERROR: 10_SIMULATION parallel step failed"
    exit 1
}


##################################################
## STEP 11: Calculate percentiles to final report
##################################################

# Define output file
PERC_OUT="$OUTPUT_FOLDER/${PGS_ID}_calibrated_chr$CHR.tsv"


# Run STEP 11 for all samples with simulation result
if [[ ! -f "$PERC_OUT"  || "$RESUME" == "false" ]]; then
    run_logged "11_CALC_PERCENTILES" Rscript $BASE_DIR/Functions/11_calculate_percentiles.R \
                                     --chr $CHR --pgs-id $PGS_ID --lai-codes $LAI_CODES --out $OUTPUT_FOLDER
else
    run_logged "11_CALC_PERCENTILES" echo " Calibrated percentiles for $PGS_ID on chromosome $CHR already exist, skipping."
fi


# Make sure this step produced the file needed for next steps
if [[ ! -f $PERC_OUT ]]; then
    echo "ERROR: Missing $PERC_OUT" | tee -a "$MAIN_LOG"
    return 1   
fi


########################
## Print final message
########################

echo -e "Main script ended. Results can be found at $OUTPUT_FOLDER\n"  | tee -a $MAIN_LOG

            
## End
