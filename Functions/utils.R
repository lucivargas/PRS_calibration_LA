# Script to define functions that will be used across the pipeline
# Author: Luciana Vargas

# Load libraries
suppressPackageStartupMessages(library(tidyverse))

# Function to flip alleles to reverse complement, given a list of sequences
reverse_complement <- function(sequence_list){
  
  # Copy list to new object
  new_list <- sequence_list
  
  for(i in 1:length(sequence_list)){
    
    # Extract each sequence from list
    each_seq <- sequence_list[i]
    
    # If sequence has more than one base, split bases
    if(nchar(each_seq) > 1){
      each_seq <- unlist(str_split(each_seq, pattern = ""))
    }
    
    # Copy sequence to new object
    new_seq <- each_seq
    
    # Iterate bases in sequence, and reverse complement them
    for(j in 1:length(new_seq)){
      new_seq[j] <- case_when(each_seq[j] == "G" ~ "C",
                              each_seq[j] == "C" ~ "G",
                              each_seq[j] == "T" ~ "A",
                              each_seq[j] == "A" ~ "T", .default = NA)
      
    }
    
    # Repaste sequence together, if previously split
    if(length(new_seq) > 1){
      new_seq <- paste(new_seq, collapse = "")
    }
    
    # Append result to sequence list
    new_list[i] <- new_seq
  }
  
  # Return rev. complement of original list
  return(new_list)
}



# Function to compare alelle mismatch between PGS and genentic data 
alelle_mismatch <- function(variant_position, pgs_alleles, data){
  
  # Get index of variant in data
  data_index <- which(data$POS == variant_position)
  
  # If more than one index, iterate indexes
  missing_list <- c()
  for(index in data_index){
    
    # Get alleles in the data
    data_alleles <- c(data$REF[index], data$ALT[index])
    
    # If no PGS allele is present in the data, flip data
    if(sum(data_alleles %in% pgs_alleles) == 0){
      data_alleles <- reverse_complement(data_alleles)
    }
    
    # Check if both alleles match, or if is missing
    is_missing <- ifelse(all(data_alleles %in% pgs_alleles), FALSE, TRUE)
    
    # Add to match list
    missing_list <- c(missing_list, is_missing)
  }
  
  # Find if BOTH of the alternative variants for this position are missing from PGS alleles
  missing <- all(missing_list)
  
  return(missing)
}

## End
