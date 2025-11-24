# Script to create a weighted PRS table for reference haplotypes
# Author: Luciana Vargas
# Date: 06/19/2025

# Calle:
# Rscript $BASE_DIR/Functions/07_apply_pgs_ref.R --ref-data $REF_FILTERED --chr $CHR \
#                             --pgs-id $PGS_ID --pgs-path $PGS_PATH --variant-list $OVERLAP_LIST \
#                             --out-file $PGS_REF --threads $THREADS --program-path $BASE_DIR


#######################
## Prepare environment
#######################

# Load packages
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--ref-data", required = TRUE, help = "Genetic data to apply PGS on")
parser$add_argument("--ref-map", required = TRUE, help = "Map of ancestry codes for reference samples")
parser$add_argument("--chr", required = TRUE, help = "Chromosome")
parser$add_argument("--pgs-id", required = TRUE, help = "PGS ID")
parser$add_argument("--pgs-path", required = TRUE, help = "Path to PGS scoring file")
parser$add_argument("--variant-list", required = TRUE, help = "List of overlapping variants between PRS ref and target data")
parser$add_argument("--out-file", required = TRUE, help = "Output file")
parser$add_argument("--threads", required = TRUE, help = "Number of threads")
parser$add_argument("--program-path", required = TRUE, help = "Path to program")
args <- parser$parse_args()

# Define variables
REF_DATA <- args$ref_data
REF_MAP <- args$ref_map
CHR <- args$chr
PGS_ID <- args$pgs_id
PGS_PATH <- args$pgs_path
OVERLAP_LIST <- args$variant_list
OUTPUT_FILE <- args$out_file
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
gen <- fread(cmd=paste0("zcat ", REF_DATA," | grep -v '^##'"), sep = "\t", header = TRUE)
gen_header <- names(gen)

# Import reference map
map <- read_tsv(REF_MAP, show_col_types = F) %>%
  rename(Sample="#Sample")

###############
## Format data
###############

cat(paste0(" Formatting data..."), sep="\n")

# Format PGS table
if(!"chr_name" %in% names(pgs)){pgs$chr_name <- pgs$hm_chr}
if(!"other_allele" %in% names(pgs)){pgs$other_allele <- pgs$hm_inferOtherAllele}
pgs <- pgs %>% 
  # Remove positions not harmonized
  filter(!is.na(hm_pos)) %>%
  # Keep only chromosome of interest
  filter(chr_name == CHR) %>%
  # Format columns to numeric
  mutate(across(c(chr_name, hm_pos, effect_weight), as.numeric))  %>%
  # Filter to keep variants present in target genetic data
  inner_join(., variants, by = c("chr_name" = "chr", "hm_pos" = "pos")) 

# Format VCF, to keep only SNPs in PGS
gen <- gen[POS %in% pgs$hm_pos]

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
remove_var <- c()
for(each_var in duplicated_vars){
  
  # Get alleles for that variant form
  data_alleles <- c(gen$REF[each_var], gen$ALT[each_var])
  pgs_alleles <-  c(gen$effect_allele[each_var], gen$other_allele[each_var])
  
  # If no PGS allele is present in the data, flip data
  if(sum(data_alleles %in% pgs_alleles) == 0){
    data_alleles <- reverse_complement(data_alleles)
  }
  
  # Check if both alleles match, or if is missing
  if(!all(data_alleles %in% pgs_alleles)){
    remove_var <- c(remove_var, each_var)
  }
}

# Remove those from genetic data
gen <- as.data.frame(gen)
gen <- gen[-remove_var,]

############################################
## Apply PRS for each individual haplotype
############################################

# Print message before start
sample_count <- length(gen_header) - 10
cat(paste0(" Applying ", PGS_ID, " to ", sample_count," samples in reference data..."), sep="\n")

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

# Perform SNP weighting using parallel processing
haplo_matrix <- foreach(i = 1:nrow(gen), .combine='cbind', .packages=c("data.table")) %dopar% {
  
  # Get alleles from PGS and data
  pgs_alleles <-  c(gen$effect_allele[i], gen$other_allele[i])
  data_alleles <-  c(gen$REF[i], gen$ALT[i])
  
  # If none of the alleles in data match PGS, flip to reverse complement
  if(!all(pgs_alleles %in% data_alleles)){
    data_alleles <- reverse_complement(data_alleles)
  }
  
  # Extract possible alleles at SNP
  ref <- data_alleles[1]
  alt <- data_alleles[2]
  risk <- gen$effect_allele[i]
  weight <- as.numeric(gen$effect_weight[i])
  
  # Determine which gen allele is the risk allele (0=REF, 1=ALT)
  risk_code <- ifelse(risk == alt, 1,
                      ifelse(risk == ref, 0, NA))
  
  # Return error for SNPs with mismatched alleles
  if (is.na(risk_code)){
    
    break
    hap2_vals <- NA
    
  } else{
    
    # Get genotypes for all individuals
    genotypes <- unlist(gen[i, geno_cols])
    hap1 <- as.integer(gsub("\\|.*", "", unlist(genotypes)))
    hap2 <- as.integer(gsub(".*\\|", "", unlist(genotypes)))
    
    # If risk allele is present, add weight 
    # If the allele is missing, add an NA
    # If alternative allele, add 0
    hap1_vals <- ifelse(hap1 == risk_code, weight, ifelse(hap1 == ".", NA, 0))
    hap2_vals <- ifelse(hap2 == risk_code, weight, ifelse(hap2 == ".", NA, 0))
    
  }
  
  # Add weighted values to main matrix
  mat <- matrix(0, nrow = 2 * length(geno_cols), ncol = 1)
  mat[seq(1, by = 2, length.out = length(geno_cols))] <- hap1_vals
  mat[seq(2, by = 2, length.out = length(geno_cols))] <- hap2_vals
  
  return(mat)
}

# Cleanup
stopCluster(cl)

# Assign SNP and haplotype IDs to matrix
haplo_matrix <- as.matrix(haplo_matrix)
rownames(haplo_matrix) <- unlist(lapply(geno_cols, function(id) c(paste0(id, "_1"), paste0(id, "_2"))))
colnames(haplo_matrix) <- gen$ID


##########################
## Impute missing weights
##########################

# If there's any missing data, impute it using the mean of each panel

# Iterate panels
panels <- unique(map$Panel)
imputation_counter <- 0
for(each_panel in panels){

  # Identify subset of data in each panel
  panel_ids <- map$Sample[map$Panel == each_panel]
  panel_ids <- paste0(panel_ids, c("_1", "_2"))
  index_panel <- row.names(haplo_matrix) %in% panel_ids
    
  # Iterate SNPs within data & impute
  for(snp in 1:ncol(haplo_matrix)){
    
     # Identify missing data for that SNP within panel
    index_missing <- which( is.na(haplo_matrix[[snp]]) & index_panel )
    
    # If there is missing data, proceed with imputation
    if(length(index_missing) > 0){
      
      # Calculate the mean weight for that SNP within panel
      mean_snp <- mean(haplo_matrix[index_panel, snp], na.rm = T)
      
      # Impute missing data using mean weights within panel
      haplo_matrix[index_missing, snp] <- mean_snp
      
      # Update counter of data points imputed
      imputation_counter <- imputation_counter + length(index_missing)
    }
    
  }
}

cat(paste0(" ", imputation_counter, " missing values were imputed using mean weights for variant within each panel"), sep="\n")

################
## Save results
################

# Format & Save haplotype matrix 
cat(paste0(" Saving results..."), sep="\n")
haplo_matrix %>% 
  # Convert matrix to data frame
  as.data.frame(., stringsAsFactors = FALSE) %>%
  # Format sample names
  rownames_to_column(var="Haplotypes") %>%
  # Save weighted PGS dataframe
  write_tsv(., file=OUTPUT_FILE)


# Print final message
cat(paste0(" ", PGS_ID, " succesfully applied to reference data. File saved at ", OUTPUT_FILE), sep = "\n")

## End
