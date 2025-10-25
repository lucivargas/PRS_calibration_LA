# Script to check if local ancestry inference (LAI) samples match to genetic data samples
# Author: Luciana Vargas
# Date: 06/18/2025

# Load packages
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyverse))

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--chr", required = TRUE, help = "Chromosome")
parser$add_argument("--target-data", required = TRUE, help = "Path to VCF file")
parser$add_argument("--target-lai", required = TRUE, help = "Path to LAI file")
parser$add_argument("--out", required = TRUE, help = "Output folder")
args <- parser$parse_args()

# Define variables
CHR <- args$chr
TARGET_DATA <- args$target_data
TARGET_LAI <- args$target_lai
OUTPUT_FOLDER <- args$out

# Print message
cat(" Checking target data...", sep="\n")


################################
## Extract LAI code assignments
################################

## Read in ancestry codes for LAI intervals (query_results.msp)
lai_codes <- readLines(TARGET_LAI, n=1)

# Remove prefix
lai_codes <- gsub("\\#Subpopulation order\\/codes\\: ", "", lai_codes)

# Split populations
lai_codes <- str_split(lai_codes, pattern="\t")[[1]]

# Split each code by '=' and convert to data frame
lai_codes <- do.call(rbind, strsplit(lai_codes, "="))
lai_codes <- as.data.frame(lai_codes, stringsAsFactors = FALSE)
names(lai_codes) <-  c("Panel", "Code")

# Save data frame for future use
write_tsv(lai_codes, file=paste0(OUTPUT_FOLDER, "/gen/",basename(TARGET_LAI),".codes"))


#######################################
## Check that all VCF samples have LAI
#######################################

# Import samples names from VCF
sample_list_vcf <- system2("bcftools", args = c("query", "-l", TARGET_DATA), stdout = TRUE)

# Import samples names from LAI
sample_list_lai <- readLines(TARGET_LAI, n = 2)[2]
sample_list_lai <- str_split(sample_list_lai, pattern="\\t")[[1]]
sample_list_lai <- unique(gsub("\\.[01]{1}$", "", sample_list_lai))

# Get list of samples that overlap between VCF and LAI
index_overlap <- which(sample_list_vcf %in% sample_list_lai)
sample_overlap <- sample_list_vcf[index_overlap]


# Check that all samples in the genetic file have LAI information
if(all(sample_list_vcf %in% sample_list_lai)){
  
  # Print success message
  cat(paste0(" ✅ Passed sample check\n", 
             "     All samples in the VCF file were found in the LAI file\n"), sep= "\n")
  
} else {
  
  # Calculate the number of samples that do overlap, if any
  overlap_n <- sum(sample_list_vcf %in% sample_list_lai)
  
  # Print WARNING message
  cat(paste0("  WARNING: Not all samples in the VCF file were found in the LAI file.\n",
      " VCF file contains ",length(sample_list_vcf), " samples, but only ",
      overlap_n, " samples were found in LAI file."), sep="\n")
  
  
  # If no samples overlap, err out
  if( overlap_n < 1){
    
    cat(paste0(" ERROR: No samples overlap between VCF file and LAI file.\n",
               " Please ensure samples match before running the pipeline."), sep="\n")
    
    # Return non-zero exit 
    quit(status = 1)  
    
  }
  
}


# Save list of overlapping samples between VCF and LAI
sample_overlap <- as.data.frame(sample_overlap)
write_tsv(sample_overlap, file=paste0(OUTPUT_FOLDER, "/snp_lists/overlap_sample_list_chr",CHR,".tsv"), col_names = F)


# Return zero exit 
quit(status = 0)


## End
