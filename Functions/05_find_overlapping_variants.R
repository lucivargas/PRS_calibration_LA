# Call:
# Rscript ./Functions/05_find_overlapping_variants.R --target-list $TARGET_LIST --ref-list $REF_LIST \
#                                                    --pgs-path $PGS_PATH --chr $CHR --out $OUTPUT_FOLDER

###################
## Set environment
###################

# Load libraries
suppressWarnings({
  suppressPackageStartupMessages({
    library(tidyverse)
    library(argparse)
  })
})

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--target-list", required = TRUE, help = "Target genetic data")
parser$add_argument("--ref-list", required = TRUE, help = "Reference genetic data")
parser$add_argument("--pgs-path", required = TRUE, help = "Path to PGS scoring file")
parser$add_argument("--chr", required = TRUE, help = "Chromosome")
parser$add_argument("--out", required = TRUE, help = "Output folder")
parser$add_argument("--program-path", required = TRUE, help = "Path to program")
args <- parser$parse_args()

# Define variables
TARGET_LIST <- args$target_list
REF_LIST <- args$ref_list
CHR <- args$chr
PGS_PATH <- args$pgs_path
OUTPUT_FOLDER <- args$out
BASE_DIR <- args$program_path

# Import pre-built functions
source(paste0(BASE_DIR, "/Functions/utils.R"))

########################
## Import variant lists
########################

## Import PGS list
pgs <- read_tsv(PGS_PATH, comment="#", col_types="cccccccccccc")
if(!"other_allele" %in% names(pgs)){pgs$other_allele <- pgs$hm_inferOtherAllele}
pgs <- pgs %>% 
  select(hm_chr, hm_pos, effect_allele, other_allele) %>%
  rename(chr="hm_chr", pos="hm_pos") %>%
  mutate(chr = gsub("^chr", "", chr)) %>%
  filter(chr == CHR) %>%
  mutate(across(c(chr, pos), as.numeric))

# Print PGS summary message
PGS_ID <- str_extract(basename(PGS_PATH),"PGS[0-9]*")
PGS_COUNT <- nrow(pgs)
cat(paste0(" The scoring file for ",PGS_ID," has ",PGS_COUNT," variants on chromosome ", CHR), sep="\n")

# Import variants from target VCF
target <- read_tsv(TARGET_LIST, show_col_types = FALSE)

# Print target message
TARGET_COUNT=nrow(target)
cat(paste0(" The target data has ", TARGET_COUNT," variants on chromosome ", CHR), sep="\n")

# Import variants from reference VCF
ref <- read_tsv(REF_LIST, show_col_types = FALSE)

# Print ref message
REF_COUNT=nrow(ref)
cat(paste0(" The reference data has ", REF_COUNT," variants on chromosome ", CHR), sep="\n")


####################################
## Identify AT/CGs in PGS, if any
####################################

pgs$AMBIGUOUS_SNP <- NA

for(i in 1:nrow(pgs)){
  alleles <- c(pgs$effect_allele[i], pgs$other_allele[i])
  is_ambiguous <- all(c("A", "T") %in% alleles) | all(c("C", "G") %in% alleles) 
  pgs$AMBIGUOUS_SNP[i] <- ifelse(is_ambiguous, TRUE, FALSE)
}

cat(paste0(" There are ", sum(pgs$AMBIGUOUS_SNP, na.rm=F),
           " ambiguous (A/T or C/G) variants in the score"), sep="\n")


##############################################
## Identify missing PGS vars in target or ref
##############################################

# Print message
cat(paste0(" Finding overlapping variants between ",PGS_ID,", reference, and target data on chromosome ",
           CHR,"..."), sep="\n")

# Create columns to save lookup results
pgs$MISSING_POS_FROM_REF <- NA
pgs$INCONSISTENT_ALLELE_FROM_REF <- NA
pgs$MISSING_POS_FROM_TARGET <- NA
pgs$INCONSISTENT_ALLELE_FROM_TARGET <- NA

# Create progress bar
pb <- txtProgressBar(min = 0, max = nrow(pgs), style = 3)

# Iterate variants in PGS
for(i in 1:nrow(pgs)){
  
  # Extract PGS variant position and alleles
  pos <- pgs$pos[i]
  pgs_alleles <-  c(pgs$effect_allele[i], pgs$other_allele[i])
  
  # Check if PGS variant is missing from the reference
  pgs$MISSING_POS_FROM_REF[i] <- ifelse(pos %in% ref$POS, FALSE, TRUE)
  
  # If variant is found, check if both PGS alleles are present (verbatim or flipped)
  if(!pgs$MISSING_POS_FROM_REF[i]){
    pgs$INCONSISTENT_ALLELE_FROM_REF[i] <- alelle_mismatch(pos, pgs_alleles, ref)
  } else{
    pgs$INCONSISTENT_ALLELE_FROM_REF[i] <- NA
  }
  
  # Check if PGS variant is missing from the target
  pgs$MISSING_POS_FROM_TARGET[i] <- ifelse(pos %in% target$POS, FALSE, TRUE)
  
  # If variant is found, check if both PGS alleles are present (verbatim or flipped)
  if(!pgs$MISSING_POS_FROM_TARGET[i]){
    pgs$INCONSISTENT_ALLELE_FROM_TARGET[i] <- alelle_mismatch(pos, pgs_alleles, target)
  } else{
    pgs$INCONSISTENT_ALLELE_FROM_TARGET[i] <- NA
  }
  
  # Print progress
  if (interactive()){
    setTxtProgressBar(pb, i)
  }
}

# Close progress bar
setTxtProgressBar(pb, i)
close(pb)


########################
## Create output lists
########################

# Where to save lists
OVERLAP_LIST <- paste0(OUTPUT_FOLDER, "/snp_lists/overlap_",PGS_ID,"_chr",CHR,"_variant_list.txt")
EXCLUDE_LIST <- paste0(OUTPUT_FOLDER, "/snp_lists/excluded_",PGS_ID,"_chr",CHR,"_variant_list.txt")

# Create list with variants to keep
pgs_overlap <- pgs %>%
  filter(! AMBIGUOUS_SNP ) %>%
  filter(! (MISSING_POS_FROM_REF | INCONSISTENT_ALLELE_FROM_REF )) %>%
  filter(! (MISSING_POS_FROM_TARGET | INCONSISTENT_ALLELE_FROM_TARGET )) %>%
  select(chr, pos)

# Print message
perc_overlap <- round((nrow(pgs_overlap) / PGS_COUNT)*100 , digits = 2)
cat(paste0(" Found ", nrow(pgs_overlap), " variants (", perc_overlap,
           "%) that overlap between PGS, target, and reference data on chromosome " ,CHR), sep="\n")

# Save list of variants to keep
write_tsv(pgs_overlap, file=OVERLAP_LIST, col_names = F)

# Create list with variants to exclude
pgs_exclude <- pgs %>%
  filter( AMBIGUOUS_SNP | MISSING_POS_FROM_REF | INCONSISTENT_ALLELE_FROM_REF | MISSING_POS_FROM_TARGET | INCONSISTENT_ALLELE_FROM_TARGET )
cat(paste0(" ", nrow(pgs_exclude), " variants were excluded."), sep="\n")

# Save list of variants to exclude
write_tsv(pgs_exclude, file=EXCLUDE_LIST)


## End
