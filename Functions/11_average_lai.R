# Script to average LAI assignments (across all chromosomes)
# Author: Luciana Vargas
# Date: 08/29/2025

#######################
## Prepare environment
#######################

# Load libraries
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--input", required = TRUE, help = "Input files samplesheet")
parser$add_argument("--out", required = TRUE, help = "Output folder")
args <- parser$parse_args()

# Define variables
INPUT <- args$input
OUTPUT_FOLDER <- args$out


###############
## Import data
###############

# Read in samplesheet
info <- read_csv(INPUT)

# Read in sample names
sample_files <- list.files(path = paste0(OUTPUT_FOLDER,"/snp_lists/"), 
                           pattern = "^overlap_sample_list", full.names = T)
sample_names <- c()
for(each_file in sample_files){
  add <- unlist(read_tsv(each_file, col_names = F))
  sample_names <- unique(c(sample_names, add))
}


###########################
## Average LAI into global
###########################

lai_header <- c("CHR", "LAI_panel", "LAI_code", sample_names)
lai_global <- as.data.frame(matrix(nrow=0, ncol=length(lai_header)))
names(lai_global) <- lai_header

# Iterate LAI results for each chromosome
for(i in 1:nrow(info)){
  
  # Extract LAI file path and chromosome nome
  each_lai <- info$TARGET_LAI[i]
  each_chr <- info$CHR[i]
  
  # Print message
  cat(paste0(" Averaging LAI assignments for chromosome ",each_chr,"..." ))
  
  # Import lai result  
  cols_info <- c("#chm",	"spos",	"epos",	"sgpos",	"egpos",	"n snps")
  cols_sample <- paste0(sample_names, c(".0", ".1"))
  lai <- fread(each_lai, select = c(cols_info, cols_sample))
  
  # Import lai codes
  LAI_CODES <- paste0(OUTPUT_FOLDER,"/gen/", basename(each_lai), ".codes")
  lai_codes <- read_tsv(LAI_CODES)
  
  # Create object to save counts
  lai_header <- c("CHR", "LAI_panel", "LAI_code", sample_names)
  add <- as.data.frame(matrix(nrow=nrow(lai_codes), ncol=length(lai_header)))
  names(add) <- lai_header
  add$CHR <- rep(each_chr, nrow(lai_codes))
  add$LAI_panel <- lai_codes$Panel
  add$LAI_code <- lai_codes$Code
  
  # Iterate samples & count
  for(each_sample in sample_names){
    
    # Find index of sample in data
    sample_haplos <- paste0(each_sample, c(".0", ".1"))
    
    # Iterate haplotype in each sample
    for(each_haplo in sample_haplos){
      
      # Iterate LAI panels
      for(j in 1:nrow(add)){
        
        # Extract panel and code 
        panel_name <- add$LAI_panel[j]
        panel_code <- add$LAI_code[j]
        
        # Find index of intervals where the lai code is present
        index_lai <- which(lai[[each_haplo]] == panel_code)
        
        # Sum the number of SNPs in all intervals with that lai code assignment
        lai_snp_count <- sum(lai[index_lai,"n snps"])
        
        # Add result to main dataframe
        add[j, each_sample] <- sum(add[j, each_sample], lai_snp_count, na.rm = T)
        
      }
    }
  }
  
  # Add chr specific results to main frame
  lai_global <- rbind(lai_global, add)
  
}


##################
## Format results
##################

# Sum counts across chromosomes
lai_global <- lai_global %>%
  select(-LAI_code) %>%
  group_by(LAI_panel) %>%
  summarise( across(-(1), ~ sum(.x, na.rm = TRUE))  ) %>%
  ungroup() 

# Format and get percentages
lai_global_perc <- lai_global %>%
  # Transpose data
  pivot_longer(cols = -(1), names_to = "Sample", values_to = "Count") %>%
  pivot_wider(names_from = "LAI_panel", values_from="Count") %>%
  # Calculate percentages of each panel
  rowwise() %>% 
  mutate(total = sum(c_across(2:ncol(.)), na.rm = TRUE)) %>%
  mutate(across(2:(ncol(.) - 1), ~ .x / total * 100, .names = "LAI_perc_{.col}")) %>%
  ungroup() %>%
  # Clean up final df
  select(Sample, starts_with("LAI_perc"))


################
## Save results
################

# Save results
lai_file <- paste0(OUTPUT_FOLDER, "/gen/target_LAI_global_percentages.tsv")
write_tsv(lai_global_perc, file=lai_file)

## End
