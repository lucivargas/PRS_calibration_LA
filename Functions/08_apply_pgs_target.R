# Script to create a weighted PRS table for target samples
# Author: Luciana Vargas
# Date: 06/25/2025

# Example usage:
# Rscript $BASE_DIR/Functions/08_apply_pgs_target.R --gen-data $TARGET_FILTERED --chr $CHR \
#                       --pgs-id $PGS_ID --pgs-path $PGS_PATH --variant-list $OVERLAP_LIST \
#                       --out-file $PGS_TARGET --out $OUTPUT_FOLDER --threads $THREADS --program-path $BASE_DIR


#######################
## Prepare environment
#######################

# Load packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--gen-data", required = TRUE, help = "Genetic data to apply PGS on")
parser$add_argument("--chr", required = TRUE, help = "Chromosome")
parser$add_argument("--pgs-id", required = TRUE, help = "PGS ID")
parser$add_argument("--pgs-path", required = TRUE, help = "Path to PGS scoring file")
parser$add_argument("--variant-list", required = TRUE, help = "List of overlapping variants between PGS, REF, and TARGET data")
parser$add_argument("--out", required = TRUE, help = "Output folder")
parser$add_argument("--out-file", required = TRUE, help = "Output file")
parser$add_argument("--threads", required = TRUE, help = "Number of threads")
parser$add_argument("--program-path", required = TRUE, help = "Path to program")
args <- parser$parse_args()

# Define variables
GEN_DATA <- args$gen_data
CHR <- args$chr
PGS_ID <- args$pgs_id
PGS_PATH <- args$pgs_path
OVERLAP_LIST <- args$variant_list
OUTPUT_FILE <- args$out_file
OUTPUT_FOLDER <- args$out
THREADS <- args$threads
BASE_DIR <- args$program_path

# Import pre-built functions
source(paste0(BASE_DIR, "/Functions/utils.R"))


###############
## Import data
###############

## Import list of variants present in target genetic data and PGS 
variants <- read.table(OVERLAP_LIST, sep="\t", header = FALSE, col.names = c("chr", "pos")) %>%
  mutate(across(c(chr, pos), as.numeric)) 

## Import PGS weights
cat(" Importing PGS weights...", sep="\n")
pgs <- read_tsv(PGS_PATH, comment="#", col_types="cccccccccccc")

## Import and parse VCF
cat(" Importing genetic data...", sep="\n")
gen <- fread(cmd=paste0("zcat ",GEN_DATA," | grep -v '^##'"), sep = "\t", header = TRUE)
gen_header <- names(gen)


###############
## Format data
###############

cat(paste0(" Formatting data..."), sep="\n")

# Format PGS table
if(!"chr_name" %in% names(pgs)){pgs$chr_name <- pgs$hm_chr}
if(!"other_allele" %in% names(pgs)){pgs$other_allele <- pgs$hm_inferOtherAllele}
pgs <- pgs %>% 
  # Remove positions not harmonized to current genome build
  filter(!is.na(hm_pos)) %>%
  # Keep only chromosome of interest
  filter(chr_name == CHR) %>%
  # Format columns to numeric
  mutate(across(c(chr_name, hm_pos, effect_weight), as.numeric))  %>%
  # Filter to keep variants present in target genetic data
  inner_join(., variants, by = c("chr_name" = "chr", "hm_pos" = "pos")) 

# Format VCF, to keep only SNPs in PGS
gen <- gen[POS %in% pgs$hm_pos]

# Remove chr prefix from gen data
gen$ID <- gsub("^chr","", gen$ID)

# Merge VCF with PGS to get weight and risk allele
pgs$POS <- pgs$hm_pos
gen <- merge(gen, pgs, by = "POS")

# Extract genotype columns from gen data
geno_cols <- gen_header[10:length(gen_header)]


########################################
## Clean alternative alleles from data
########################################

# Some variants may have two rows for the same SNP in the data, with alternative alleles for each
# We will need to identify which one is relevant for the PGS, and remove the other

# Identify instances with two alternative forms in the same position
duplicated_vars <- which(gen$POS %in% gen$POS[duplicated(gen$POS)])

# Iterate variants and list rows that should be removed (those with allele forms that do not match PGS)
if(length(duplicated_vars)>0){
  remove_var <- c()
  for(each_var in duplicated_vars){
    
    # Get alleles for index of variant 
    data_alleles <- c(gen$REF[each_var], gen$ALT[each_var])
    pgs_alleles <-  c(gen$effect_allele[each_var], gen$other_allele[each_var])
    
    # If no PGS allele is present in the data, flip data
    if(sum(data_alleles %in% pgs_alleles) == 0){
      data_alleles <- reverse_complement(data_alleles)
    }
    
    # Check if both alleles match, or if is missing, and add to removal list
    if(!all(data_alleles %in% pgs_alleles)){
      remove_var <- c(remove_var, each_var)
    }
  }
  
  # Remove those from genetic data
  gen <- as.data.frame(gen)
  gen <- gen[-remove_var,]
}


####################################
## Apply PGS in batches of samples
####################################

# Start new results file with just the header
pgs_header <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(pgs_header) <- c("Sample", paste0(PGS_ID,"_chr",CHR))
write_tsv(pgs_header, file=OUTPUT_FILE, append = F)

# Detect available cores
detected_cores <- max(1, detectCores() - 1)
n_cores <- ifelse(detected_cores < THREADS , detected_cores, THREADS)
cat(paste0(" Using ", n_cores, " cores for parallelization with PSOCK clusters..."), sep="\n")

# Define type of cluster to create
# cl_type <- ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK")
# FORK doesnt always work in linux, so I'll force PSOCK for now

# Register parallel backend
cl <- makeCluster(as.numeric(n_cores), type = "PSOCK")
print(cl)
registerDoParallel(cl)

# Create batches, including 1,000 samples at a time
sample_count <- length(geno_cols) 
batches <- data.frame(start=seq(from=1, to=sample_count, by=100),
                      stringsAsFactors = F)
batches$end <- c(batches$start[-1] - 1,  sample_count)

# Print message before start
cat(paste0(" Applying ", PGS_ID, " to ", sample_count, " samples in target data..."), sep="\n")

# Create progress bar
pb <- txtProgressBar(min = 0, max = nrow(batches), style = 3)

# Apply PGS in each batch
for(j in 1:nrow(batches)){
  
  # Subset data
  index_start <- batches$start[j]
  index_end <- batches$end[j]
  batch_cols <- geno_cols[index_start:index_end]
  gen_batch <- as.data.frame(gen[, batch_cols])
  
  # Export variant columns
  info <- gen %>%
    select(ID, REF, ALT, effect_allele, other_allele, effect_weight)
  
  # Perform SNP weighting using parallel processing
  pgs_matrix <- foreach(i = 1:nrow(gen_batch), 
                        .combine='cbind', # .export=c("info","reverse_complement"),
                        .packages=c("data.table", "tidyverse")) %dopar% {
                          
    # Initialize object to save weights
    mat <- matrix(NA, nrow = length(batch_cols), ncol = 1)                   
                          
    # Get alleles from PGS and data
    pgs_alleles <-  c(info$effect_allele[i], info$other_allele[i])
    data_alleles <-  c(info$REF[i], info$ALT[i])
    
    # If none of the alleles in data match PGS, flip to reverse complement
    if(!all(pgs_alleles %in% data_alleles)){
      data_alleles <- reverse_complement(data_alleles)
    }
    
    # Extract possible alleles at SNP
    ref <- data_alleles[1]
    alt <- data_alleles[2]
    risk <- info$effect_allele[i]
    weight <- as.numeric(info$effect_weight[i])
    
    # Determine which gen allele is the PGS risk allele (0=REF, 1=ALT)
    risk_code <- ifelse(risk == alt, 1,
                        ifelse(risk == ref, 0, NA))
    
    # Get genotype dosages for all individuals
    dosages <- as.numeric(gen_batch[i,])
    
    # Dosage in Minimac3: Posterior probability of the alternate allele at that site
    # Determine if ALT is risk allele, if not, flip dosage to reflect risk allele probability
    if(alt == risk){
      dosage_aligned <- dosages
    } else if(ref == risk){
      dosage_aligned <-  2 - dosages
    }
    
    # Check if all alleles match 
    if(!all(pgs_alleles %in% data_alleles)){
      cat(paste0("Error in calculating PGS on target data, SNP ID ", info$ID[i], ", alleles do not match between PGS and data"), 
          file=paste0(OUTPUT_FOLDER,"/logs/error.log"), append=TRUE)
    } else{
      
      # Multiply dosage of risk allele by weight
      vals <- dosage_aligned * weight
      
      # Add weighted values to batch matrix
      mat[,1] <- vals
      
    }
    
    # Return scores
    return(mat)
  }
  
  # Assign SNP and haplotype IDs to batch matrix
  pgs_matrix <- as.matrix(pgs_matrix)
  rownames(pgs_matrix) <- batch_cols
  colnames(pgs_matrix) <- gen$ID
  
  # Format & Save SNP-based PGS matrix
  pgs_matrix %>%
    # Convert matrix to data frame
    as.data.frame(., stringsAsFactors = FALSE) %>%
    # Convert all to numbers
    mutate(across(everything(), as.numeric)) %>%
    # Sum all variants weights to one value per chromosome
    rowwise() %>%
    mutate(total = sum(c_across(everything()), na.rm=T)) %>%
    ungroup() %>%
    # Keep only the total sum
    select(total) %>%
    # Format sample names
    add_column(Samples=batch_cols, .before=1) %>%
    # Save weighted PGS dataframe
    write_tsv(., file=OUTPUT_FILE, append = T)
  
  # Print progress
  setTxtProgressBar(pb, j)
}

# Close progress bar
close(pb)

# Cleanup
stopCluster(cl)

# Print final message
cat(paste0(PGS_ID, " succesfully applied to data. File saved at ", OUTPUT_FILE), sep="\n")

## End
