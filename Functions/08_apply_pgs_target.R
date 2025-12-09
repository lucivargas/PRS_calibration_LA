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
suppressWarnings({
  suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(tidyverse)
    library(doParallel)
    library(foreach)
  })
})

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
parser$add_argument("--mem-mb", required = TRUE, help = "Memory available (in MB)")
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
MEMORY <- args$mem_mb
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

# Export variant columns into new object
info <- gen %>%
  select(ID, REF, ALT, effect_allele, other_allele, effect_weight)

# Convert all genotype columns to numeric type
gen <- as.data.table(gen)
gen[, (geno_cols) := lapply(.SD, as.numeric), .SDcols = geno_cols]


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


#############################################
## Register parallelization & create batches
#############################################

# Detect number of available cores
detected_cores <- max(1, detectCores() - 1)
THREADS <- as.numeric(THREADS)
n_cores <- min(THREADS, detected_cores)
cat(paste0(" Using ", n_cores, " cores for parallelization with PSOCK clusters..."), sep="\n")

# Define type of cluster to create
# cl_type <- ifelse(.Platform$OS.type == "unix", "FORK", "PSOCK")
# FORK doesnt always work in linux, so I'll force PSOCK for now

# Register parallel backend
cl <- makeCluster(as.numeric(n_cores), type = "PSOCK")
cl_info <- capture.output(print(cl))
cat(paste0(" ", cl_info), sep="\n")
registerDoParallel(cl)

# Find the main R memory usage
# pid <- Sys.getpid()
# rss_bytes <- as.numeric(system(sprintf("ps -o rss= -p %d", pid), intern = TRUE)) * 1024
# 
# # Get total RAM available to our processes
# safety_margin <- 0.7
# MEMORY <- as.numeric(MEMORY)
# MEMORY_bytes <- MEMORY * (1024^2)
# total_ram_available <- ( MEMORY_bytes * safety_margin ) - rss_bytes
# cat(paste0(" Total RAM available in the ", n_cores, " cores: ", 
#            round(total_ram_available / (1024^3), digits=2), " GB"), sep="\n")
# 
# # Estimate batch size as a function of number of cores, data size and RAM needed to process our data
# bytes_in_gen <- as.numeric(object.size(gen))
# total_ram_available <- total_ram_available - bytes_in_gen
# bytes_per_column <- bytes_in_gen / length(geno_cols)
# bytes_per_column_incl_temp_objects <- bytes_per_column * 20  # accounts for memory used in temp objects when scoring
# batch_size <- (total_ram_available * safety_margin ) / bytes_per_column_incl_temp_objects
# batch_size <- floor(batch_size * safety_margin)
batch_size <- 1000

# Create batches 
sample_count <- length(geno_cols) 
batches <- data.frame(start=seq(from=1, to=sample_count, by=batch_size), stringsAsFactors = F)
batches$end <- c(batches$start[-1] - 1,  sample_count)
cat(paste0(" Data will be scored in ",nrow(batches)," batch(es) of ", batch_size, " samples."), sep="\n")


####################################
## Apply PGS in batches of samples
####################################

# Print message before start
cat(paste0(" Applying ", PGS_ID, " to ", sample_count, " samples in target data..."), sep="\n")

# Start new results file with just the header
pgs_header <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(pgs_header) <- c("Sample", paste0(PGS_ID,"_chr",CHR))
write_tsv(pgs_header, file=OUTPUT_FILE, append = F)

# Create progress bar
pb <- txtProgressBar(min = 0, max = nrow(batches), style = 3)

# Apply PGS in each batch
for(j in 1:nrow(batches)){
  
  cat(paste0(" Processing batch ", j, "/", nrow(batches), "..."), sep="\n")
  
  # Subset data
  index_start <- batches$start[j]
  index_end <- batches$end[j]
  batch_cols <- geno_cols[index_start:index_end]
  gen_batch <- gen[, batch_cols, drop = FALSE]
  
  # Perform SNP weighting using parallel processing
  pgs_list <- foreach(i = 1:nrow(gen_batch), 
                      .export = c("gen_batch", "batch_cols", "reverse_complement", "info"),
                      .packages=c("data.table", "tidyverse", "dplyr")) %dopar% {
                          
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
    dosages <- as.numeric(gen_batch[i, ])
    
    # Dosage in Minimac3: Posterior probability of the alternate allele at that site
    # Determine if ALT is risk allele, if not, flip dosage to reflect risk allele probability
    if(alt == risk){
      dosage_aligned <- dosages
    } else if(ref == risk){
      dosage_aligned <-  2 - dosages
    }
    
    # Check if all alleles match 
    if(!all(pgs_alleles %in% data_alleles)){
      
      # If not, return NA
      cat(paste0(" Error in calculating PGS on target data for SNP ID ", info$ID[i], 
                 ": alleles do not match between PGS (", paste(pgs_alleles, collapse = "/"),
                 ") and data (", paste(data_alleles, collapse = "/"),
                 ")"), sep = "\n")
      vals <- NA
      
    } else{
      
      # If so, multiply dosage of risk allele by weight
      vals <- dosage_aligned * weight
      
    }
    
    # Return scores
    return(vals)
  }
  
  # Bind the results into a matrix
  pgs_matrix <- do.call(cbind, pgs_list)
  
  # Assign SNP and haplotype IDs to batch matrix
  rownames(pgs_matrix) <- batch_cols
  colnames(pgs_matrix) <- gen_batch$ID
  
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
