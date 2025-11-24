# Script to calculate PGS percentiles bases on simulated reference population 
# Author: Luciana Vargas
# Date: 08/06/2025


#######################
## Prepare environment
#######################

# Load packages
suppressWarnings({
  suppressPackageStartupMessages({
    library(argparse)
    library(data.table)
    library(tidyverse)
  })
})

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--chr", required = TRUE, help = "Chromosome")
parser$add_argument("--pgs-id", required = TRUE, help = "PGS ID")
parser$add_argument("--lai-codes", required = TRUE, help = "Codes for panels used in Local Ancestry Inference")
parser$add_argument("--out", required = TRUE, help = "Output folder")
args <- parser$parse_args()

# Define variables
CHR <- args$chr
PGS_ID <- args$pgs_id
LAI_CODES <- args$lai_codes
OUTPUT_FOLDER <- args$out

# Print message before start
cat(paste0(" Calculating percentiles of PGS distribution for each target sample..."), sep=" ")


###############
## Import data
###############

# Define location of results
sim_files <- list.files(path = paste0(OUTPUT_FOLDER, "pgs/simulations/"),
                        pattern = paste0("simulation\\_",PGS_ID,".*chr",CHR,"\\.tsv"), 
                        full.names = T)

# Read in ancestry codes for LAI intervals (query_results.msp.codes)
lai_codes <- read_tsv(LAI_CODES, show_col_types = FALSE)


############################
## Get average LAI per chr
############################

# Create object to save results
lai <- c()

for(each_sim in sim_files){
  
  # Extract sample name  for sample i
  sample_name <- gsub(paste0("^simulation\\_PGS[0-9]*\\_|\\_chr.*\\.tsv$"), "", basename(each_sim))
  
  # Define path to their expanded LAI result
  lai_file <- paste0(OUTPUT_FOLDER, "gen/lai//lai_",PGS_ID,"_", sample_name, "_chr",CHR,".tsv")
  
  # Import LAI results (per SNP)
  lai_sample <- read_tsv(lai_file, col_types = "cccc")
  
  # Iterate LAI panels
  for(i in 1:nrow(lai_codes)){
    
    # Extract panel and code 
    panel_name <- lai_codes$Panel[i]
    panel_code <- lai_codes$Code[i]
    
    # Get average LAI across chr
    haplo1 <- sum(lai_sample$haplo1 == panel_code)/nrow(lai_sample)
    haplo2 <- sum(lai_sample$haplo2 == panel_code)/nrow(lai_sample)
    global_lai <- mean(haplo1, haplo2)
    
    # Summarize results 
    add <- data.frame(Sample=sample_name,
                      panel_name=panel_name,
                      panel_perc=global_lai)
    
    # Add to main frame
    lai <- rbind(lai, add)
  }
}

# Format LAI
lai <- lai %>%
  pivot_wider(values_from = panel_perc, names_from=panel_name, names_prefix = "LAI_perc_")

#####################
## Summarize results
#####################

# Extract results from simulation files
sim_df <- data.frame(Sample=character(0), CHR=character(0), PGS_ID=character(0), PGS_RAW=character(0), 
                     PGS_SIM=character(0), PGS_BOOT=character(0), PGS_CAL=character(0))

for(each_sim in sim_files){
  
  # Import simulation results for sample i
  result <- read_tsv(each_sim, show_col_types = FALSE)
  
  # Define column of interest
  pgs_col <- paste0(PGS_ID, "_chr", CHR)
  
  # Collapse the haplotypes into one score per sample
  result <- result %>%
    separate(col = Sample, into = c("Sample_ID", "Haplotype"), remove = F) %>%
    group_by(Sample_ID) %>%
    mutate("{pgs_col}" := sum(.data[[pgs_col]])) %>%
    ungroup 
  result$Sample[-1] <- result$Sample_ID[-1]
  result <- result %>% select(-Sample_ID, -Haplotype)
  
  # Extract sample name  for sample i
  sample_name <- gsub(paste0("^simulation\\_PGS[0-9]*\\_|\\_chr.*\\.tsv$"), "", basename(each_sim))
  
  # Extract raw PGS
  index_sample <- which(result$Sample == sample_name)
  pgs_raw <- unlist(result[index_sample, pgs_col])
  
  # Find percentile among simulated samples
  index_sim <- grep("^sim", result$Sample)
  pgs_sim_dist <- ecdf(unlist(result[index_sim, pgs_col]))
  pgs_sim <- pgs_sim_dist(pgs_raw) 
  
  # Find percentile among bootstraped samples
  index_boot <- grep("^boot", result$Sample)
  pgs_boot_dist <- ecdf(unlist(result[index_boot, pgs_col]))
  pgs_boot <- pgs_boot_dist(pgs_raw) 
  
  # Find percentile among both simulated and bootstraped samples
  pgs_cal_dist <- ecdf(unlist(result[c(index_sim,index_boot), pgs_col]))
  pgs_cal <- pgs_cal_dist(pgs_raw) 
  
  # Add results to main frame
  add <- data.frame(Sample=sample_name, CHR=CHR, PGS_ID=PGS_ID, 
                    PGS_RAW=pgs_raw, 
                    PGS_SIM=pgs_sim, PGS_BOOT=pgs_boot, PGS_CAL=pgs_cal)
  sim_df <- rbind(sim_df, add)
  
}

# Join global LAI to results
sim_df <- left_join(sim_df, lai, by="Sample")


# Also report the PGS percentile when compared to the overall population
pgs_pop_dist <- ecdf(sim_df$PGS_RAW)
pgs_pop <- pgs_pop_dist(sim_df$PGS_RAW) 
sim_df <- sim_df %>%
  add_column(PGS_POP=pgs_pop, .after="PGS_RAW")

# Save results
write_tsv(sim_df, file=paste0(OUTPUT_FOLDER, "/", PGS_ID, "_calibrated_chr", CHR, ".tsv"))

## End
