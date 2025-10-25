# Script to create simulated reference population based on individual local ancestry pattern
# Author: Luciana Vargas
# Date: 06/20/2025

#######################
## Prepare environment
#######################

# Load packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(boot))

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--chr", required = TRUE, help = "Chromosome")
parser$add_argument("--sample", required = TRUE, help = "Sample ID")
parser$add_argument("--lai-codes", required = TRUE, help = "Codes for panels used in Local Ancestry Inference")
parser$add_argument("--lai-sample", required = TRUE, help = "Matched LAI for PGS variants for Sample ID")
parser$add_argument("--pgs-id", required = TRUE, help = "PGS ID")
parser$add_argument("--pgs-ref", required = TRUE, help = "Weighted PGS matrix based on reference haplotypes")
parser$add_argument("--pgs-target", required = TRUE, help = "Weighted PGS matrix based on target samples")
parser$add_argument("--ref-map", required = TRUE, help = "Map of ancestry codes for reference samples")
parser$add_argument("--nsimulations", required = TRUE, help = "Number of simulations")
parser$add_argument("--nbootstraps", required = TRUE, help = "Number of bootstraps")
parser$add_argument("--seed", required = TRUE, help = "Random seed")
parser$add_argument("--out", required = TRUE, help = "Output folder")
args <- parser$parse_args()

# Define variables
CHR <- args$chr
SAMPLE <- args$sample 
LAI_CODES <- args$lai_codes
LAI_SAMPLE <- args$lai_sample
PGS_ID <- args$pgs_id
PGS_REF <- args$pgs_ref
PGS_TARGET <- args$pgs_target
REF_MAP <- args$ref_map
NSIMULATIONS <- args$nsimulations
NBOOTSTRAPS <- args$nbootstraps
RANDSEED <- args$seed
OUTPUT_FOLDER <- args$out

# Set seed
set.seed(RANDSEED)

# Print message before start
cat(paste0(" Creating reference simulations for target sample ",SAMPLE,"..."), sep="\n")


###############
## Import data
###############

# Import LAI for sample i
lai <- read_tsv(file=LAI_SAMPLE, show_col_types = F)

# Import LAI codes
lai_codes <- read_tsv(LAI_CODES, show_col_types = FALSE)

# Import PGS from target
score_name <- paste0(PGS_ID, "_chr", CHR)
target <- fread(PGS_TARGET,  select = c("Sample", score_name)) %>%
  # Filter sample of interest
  filter(Sample == SAMPLE)

# Import PGS from reference
ref <- read_tsv(PGS_REF, show_col_types = F)

# Import reference map
map <- read_tsv(REF_MAP, show_col_types = F) %>%
  rename(Sample="#Sample")


###############
## Format data
###############

# Join reference data and map
ref <- ref %>%
  mutate(Sample = gsub("\\_[12]{1}$", "", Haplotypes)) %>%
  left_join(., map, by="Sample") %>%
  select(c("Haplotypes","Sample", "Panel", everything()))

# Add LAI ancestry codes to match samples in reference data
ref <- ref %>%
  left_join(lai_codes, ., by="Panel") 

# Extract list of positions in the ref data
# Remove SNP ID prefixes (chr)
ref_pos <- gsub("^[0-9]{1,2}\\:", "",names(ref))
# Remove SNP ID suffixes (alleles)
ref_pos <- gsub("[ATCG:]*$", "",ref_pos)

#####################
## Start simulation
#####################

# Create simulation codes for each haplotype
max_digits <- nchar(NSIMULATIONS)
sim_codes <- paste0("sim",  sprintf(paste0("%0",max_digits,"d"), 1:NSIMULATIONS))
sim_codes_haplo <- expand.grid(c("1", "2"), sim_codes)
sim_codes_haplo <- paste(sim_codes_haplo$Var2, sim_codes_haplo$Var1, sep = "_")

# Create a dataframe to store results
sim <- data.frame(Sample=sim_codes_haplo, stringsAsFactors = F)
sim[[score_name]] <- 0

# Add sample of interest
sim <- rbind(target, sim)

# Iterate each haplotype (haplo1 or haplo2) in sample i
for(each_haplo in c(1, 2)){
  
  # Define haplotype name
  haplo_name <- paste0("haplo", each_haplo)
  
  # Find LAI blocks (continuous segments where ancestry assignments are homogeneous)
  lai1 <- lai %>%
    mutate(block = cumsum(
      lag(haplo1, default = first(haplo1)) != haplo1)) %>%
    select(chr, pos, all_of(haplo_name), block)
  
  # Define simulation codes corresponding to haplotype each_haplo
  haplo_codes <- paste0(sim_codes, "_", each_haplo)
  
  # Iterate each LAI block in each haplotype
  blocks <- unique(lai1$block)
  for(each_block in blocks){
    
    # Get index of LAI rows in the block
    index_lai <- lai1$block == each_block
    
    # Extract all positions within block (block_pos)
    block_pos <- lai1$pos[index_lai]
    
    # Extract LAI assignment for the block
    block_lai <- unique(unlist(lai1[index_lai, haplo_name]))
    
    # Find columns in PGS_REF corresponding to the block
    ref_cols <- which(ref_pos %in% block_pos)
    ref_index <- which(ref$Code == block_lai)
    
    # From PGS_REF, extract block positions & ancestry-matched individuals 
    ref1 <- as.data.frame(ref[ref_index, ref_cols])
    
    # Sample WITH replacement NSIMULATIONS from reference data (for this haplotype and block)
    n_sim <- min(NSIMULATIONS, nrow(ref1))
    sim_index <- sample(1:nrow(ref1), size = n_sim, replace = TRUE)
    block_pgs <- as.data.frame(ref1[sim_index,])
    row.names(block_pgs) <- haplo_codes
     
    # Sum scores for variants in the block
    block_pgs <- rowSums(block_pgs)
    
    # Sum score in this block to overall pgs in results matrix
    results_index <- which(sim$Sample %in% haplo_codes)
    old_pgs <- unlist(sim[results_index, ..score_name])
    new_pgs <- old_pgs + block_pgs
    
    # Add sum to results matrix
    sim[results_index, score_name] <- new_pgs
    
  }
  
}


###################
## Start bootstrap
###################

NBOOTSTRAPS <- as.numeric(NBOOTSTRAPS)

if(NBOOTSTRAPS > 0){
  
  # Create bootstrap codes for each haplotype
  max_digits <- nchar(NBOOTSTRAPS)
  boot_codes <- paste0("boot", sprintf(paste0("%0",max_digits,"d"), 1:NBOOTSTRAPS))
  
  # For each haplotype j
  boot <- c()
  for(each_haplo in c(1, 2)){
    
    # Identify each haplotype in simulated data
    haplo_index <- grep(paste0("\\_", each_haplo), sim$Sample)
    
    # Sample with replacement NBOOTSTRAP values from simulated data
    boot_index <- sample(haplo_index, size = NBOOTSTRAPS, replace = TRUE)
    
    # Extract sample from simulated data
    add <- sim[boot_index,]
    add$Sample <- paste0(boot_codes, "_", each_haplo)
    
    # Join results to main frame
    boot <- rbind(boot, add)
  }
  
  # Order bootstrap samples
  boot <- arrange(boot, Sample)
  
  # Add bootstraps to simulated data
  results <- rbind(sim, boot)
  
} else {
  results <- sim
}


################
## Save results
################

# Save simulation results
SIM_SAMPLE <- paste0(OUTPUT_FOLDER, "/pgs/simulations/simulation_",PGS_ID,"_",SAMPLE,"_chr",CHR,".tsv")
write_tsv(results, file=SIM_SAMPLE)

# Print final message
cat(" Done.", sep="\n")

## End
