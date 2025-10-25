# Script to match local ancestry inference to genetic data 
# Description: this script allows for local ancestry inference (LAI) that was done on array data to be
# matched up to more complete genetic data (e.g., imputed) by using genetic distance map to fill the gaps
# Author: Luciana Vargas
# Date: 06/18/2025


#######################
## Prepare environment
#######################

# Load packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--sample", required = TRUE, help = "Sample name")
parser$add_argument("--target-lai", required = TRUE, help = "Path to LAI file")
parser$add_argument("--chr", required = TRUE, help = "Chromosome")
parser$add_argument("--pgs-id", required = TRUE, help = "PGS ID")
parser$add_argument("--variant-list", required = TRUE, help = "List of overlapping variants between PRS and target data")
parser$add_argument("--genomeDist", required = TRUE, help = "File with genetic distances")
parser$add_argument("--out", required = TRUE, help = "Output folder")
args <- parser$parse_args()

# Define variables
SAMPLE <- args$sample
TARGET_LAI <- args$target_lai
CHR <- args$chr
PGS_ID <- args$pgs_id
OVERLAP_LIST <- args$variant_list
GENOME_DIST <- args$genomeDist
OUTPUT_FOLDER <- args$out

# Print message before start
cat(paste0(" Matching ", PGS_ID, " to local ancestry for target sample ",SAMPLE,"..."), sep=" ")


###############
## Import data
###############

## Read in local ancestry inference (LAI) intervals assignment (query_results.msp) for sample i
# Define column names to read in
cols_info <- c("#chm",	"spos",	"epos",	"sgpos",	"egpos",	"n snps")
cols_sample <- paste0(SAMPLE, c(".0", ".1"))
# Import lai
lai <- fread(TARGET_LAI, select = c(cols_info, cols_sample))
# Format lai
lai <- as.data.frame(lai)
lai <- lai %>%
  mutate(`#chm` = gsub("^chr", "", `#chm`)) %>%
  mutate(across(everything(), as.numeric)) %>%
  rename(chr = `#chm`, n_snps=`n snps`,
         haplo1 := !!sym(cols_sample[1]),
         haplo2 := !!sym(cols_sample[2]))

## Read in list of PRS variants
vars <- readLines(OVERLAP_LIST)

## Read in genetic distances map
gen_dist <- fread(cmd=paste0("gzip -dc ",GENOME_DIST), header = TRUE)
gen_dist <- as.data.frame(gen_dist)
# Keep only data for chromosome of interest
gen_dist <- gen_dist %>%
  filter(chr==CHR) %>%
  mutate(across(everything(), as.numeric))


########################################################################
## Collapse intervals if LAI assignment doesn't change for individual i 
########################################################################

lai <- lai %>%
  mutate(group = cumsum(
    lag(haplo1, default = first(haplo1)) != haplo1 |
      lag(haplo2, default = first(haplo2)) != haplo2)) %>%
  group_by(chr, haplo1, haplo2, group) %>%
  summarise(spos = min(spos),
            epos = max(epos),
            sgpos = min(sgpos),
            egpos = max(egpos),
            n_snps = sum(n_snps),
            .groups = "drop") %>%
  select(chr, spos, epos, sgpos, egpos, n_snps, haplo1, haplo2) %>%
  arrange(chr, spos)

###############################################
## Create data frame for individual i and PRS p 
###############################################

# Create new data frame to store info 
df <- data.frame(variant=vars, stringsAsFactors = F)
df <- df %>%
  separate(col=variant, into=c("chr", "pos")) %>%
  mutate(across(c(chr, pos), as.numeric)) %>%
  mutate(haplo1=NA, haplo2=NA) 
  

# Iterate each SNP in the PRS score, and find LAI assignment
for(v in 1:nrow(df)){
  
  # Extract SNP position
  query_pos <- df$pos[v]
  
  # Find whether SNP is WITHIN existing lai interval
  match_index <- which(lai$spos <= query_pos & query_pos <= lai$epos)
  
  # Find whether SNP is OUTSIDE all existing lai interval (flanking tails)
  # That is, before first interval or after last interval 
  # If so, just choose the last interval
  flanking_prev <- all(lai$epos > query_pos)
  flanking_last <- all(lai$spos < query_pos)
  if(flanking_prev | flanking_last){
    match_index <- ifelse(flanking_prev, 1,
                          ifelse(flanking_last, nrow(lai)))
  }
  
  # If WITHIN lai interval, extract LAI assignment
  if(length(match_index) > 0) {
    
    # If the position has more then 1 interval match, choose the interval with more SNPs
    if(length(match_index) > 1) {
      lai_sizes <- lai$n_snps[match_index]
      larger_lai <- which(lai_sizes == max(lai_sizes))
      match_index <- match_index[larger_lai]
      
      # If it is still tied, pick the first one
      if(length(match_index) > 1) {
        match_index <- match_index[1]
      }
    }
    
    # Extract LAI assignment
    haplo1 <- lai[match_index, "haplo1"]
    haplo2 <- lai[match_index, "haplo2"]
    
  # If OUTSIDE lai interval, find which interval is closer, and extract LAI assignment
  } else {
    
    ## Extract flanking positions to the gap between intervals
    match_index <- NA
    
    # Find the end position of the last interval
    prev_interval_index <- max(which(lai$epos < query_pos))
    prev_interval_epos <- lai$epos[prev_interval_index]
    
    # Find the start position of the next interval
    next_interval_index <- min(which(lai$spos > query_pos))
    next_interval_spos <- lai$spos[next_interval_index]
    
    ## Check the genetic distance map for distances
    # Find genomic position closest to query SNP 
    closest_query <- gen_dist$position[which.min(abs(gen_dist$position - query_pos))]
    query_gpos <-  gen_dist$`Genetic_Map(cM)`[gen_dist$position == closest_query]
    
    # Find genomic position closest to previous interval
    closest_prev <- gen_dist$position[which.min(abs(gen_dist$position - prev_interval_epos))]
    prev_gpos <-  gen_dist$`Genetic_Map(cM)`[gen_dist$position == closest_prev]
    
    # Find genomic position closest to next interval
    closest_next <- gen_dist$position[which.min(abs(gen_dist$position - next_interval_spos))]
    next_gpos <-  gen_dist$`Genetic_Map(cM)`[gen_dist$position == closest_next]
    
    # Get genetic distances to query SNP
    gdist_to_prev <- abs(query_gpos - prev_gpos)
    gdist_to_next <- abs(query_gpos - next_gpos)
    
    # Get physical distances to query SNP (as a tie criteria)
    dist_to_prev <- abs(query_pos - prev_interval_epos)
    dist_to_next <- abs(query_pos - next_interval_spos)
    
    ## Identify which interval is closest to query SNP
    is_previous <- case_when(gdist_to_prev < gdist_to_next ~ TRUE,
                             gdist_to_prev > gdist_to_next ~ FALSE,
                             dist_to_prev  < dist_to_next  ~ TRUE,
                             dist_to_prev  > dist_to_next  ~ FALSE,
                             TRUE ~ NA)
    if(is.na(is_previous)){break}
    
    # Get index of closest interval
    match_index <- ifelse(is_previous, prev_interval_index, next_interval_index)
    
    # Extract LAI assignment
    haplo1 <- lai[match_index, "haplo1"]
    haplo2 <- lai[match_index, "haplo2"]
  }
  
  # Add LAI assignment to data frame
  df$haplo1[v] <- haplo1
  df$haplo2[v] <- haplo2
}

# Format dataframe
df <- as.data.frame(df)
df$haplo1 <- as.numeric(df$haplo1)
df$haplo2 <- as.numeric(df$haplo2)

# Save dataframe
LAI_SAMPLE=paste0(OUTPUT_FOLDER,"/gen/lai/lai_",PGS_ID,"_",SAMPLE,"_chr",CHR,".tsv")
write_tsv(df, file=LAI_SAMPLE)

# Print final message
cat(" Done.", sep="\n")

## End
