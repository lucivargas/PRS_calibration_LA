# Script to vizualize and extract list of GWAS hits for harmonization
# Author: Luciana Vargas
# Date: 11/06/2025

#######################
## Prepare environment
#######################

# Load libraries
suppressWarnings({
  suppressPackageStartupMessages({
    library(tidyverse)
    library(argparse)
  })
})

# Parse named arguments
parser <- ArgumentParser()
parser$add_argument("--gwas-out", required = TRUE, help = "Prefix to GWAS output")
parser$add_argument("--pheno-file", required = TRUE, help = "Path to phenotype file")
parser$add_argument("--pca-file", required = TRUE, help = "Path to PCA file")
parser$add_argument("--overlap-list", required = TRUE, help = "Path to list of overlapping SNPs")
parser$add_argument("--exclude-list", required = TRUE, help = "Path to list of excluded SNPs")
args <- parser$parse_args()

# Define variables
GWAS_OUT <- args$gwas_out
PHENO_FILE <- args$pheno
PCA_FILE <- args$pca_file
OVERLAP_LIST <- args$overlap_list
EXCLUDE_LIST <- args$exclude_list


#############################
## Import GWAS & PCA results
#############################

# GWAS_OUT <- "/scratch/alpine/ldebritovargas@xsede.org/PRS_calibration_LA/PGS002802/gen/gwas/postgwas_chr22_PGS002802"
# PHENO_FILE <- "/scratch/alpine/ldebritovargas@xsede.org/PRS_calibration_LA/PGS002802/gen/gwas/pheno_chr22_PGS002802.txt"
# PCA_FILE <- "/scratch/alpine/ldebritovargas@xsede.org/PRS_calibration_LA/PGS002802/gen/gwas/merged_chr22_PGS002802.eigenvec"
# OVERLAP_LIST <- "/scratch/alpine/ldebritovargas@xsede.org/PRS_calibration_LA/PGS002802//snp_lists/overlap_PGS002802_chr22_variant_list.txt"
# EXCLUDE_LIST <- "/scratch/alpine/ldebritovargas@xsede.org/PRS_calibration_LA/PGS002802//snp_lists/excluded_PGS002802_chr22_variant_list.txt"

# Define folder where results were saved
gwas_folder <- dirname(GWAS_OUT)

# Import results of GWAS on batch effect
GWAS_FILE <- paste0(GWAS_OUT, ".batch.glm.logistic.hybrid")
gwas <- read_tsv(GWAS_FILE, show_col_types = FALSE) %>%
  rename(CHR = "#CHROM", LOG_OR_SE="LOG(OR)_SE") %>%
  mutate(ID = paste0(CHR,":", POS))

# Import phenotype
pheno <- read_delim(PHENO_FILE, show_col_types = FALSE)

# Import PCA results
pca <- read_tsv(PCA_FILE, show_col_types = FALSE) %>%
  rename(IID="#IID")

# Import Freq results
freq_target <- read_tsv(paste0(GWAS_OUT, "_target.afreq"), show_col_types = FALSE) %>%
  rename(ALT_FREQS_TARGET="ALT_FREQS") %>% select(-OBS_CT) 
freq_ref <- read_tsv(paste0(GWAS_OUT, "_reference.afreq"), show_col_types = FALSE) %>%
  rename(ALT_FREQS_REF="ALT_FREQS") %>% select(-OBS_CT) 
freq <- freq_target %>%
  left_join(., freq_ref, by=c("#CHROM", "ID", "REF", "ALT")) %>%
  rename(CHR="#CHROM") %>%
  separate(col=ID, into=c(NA, "POS"), sep = ":", remove = F, extra = "drop") %>%
  mutate(POS=as.numeric(POS)) %>%
  mutate(ID = paste0(CHR, ":", POS))
rm(freq_target, freq_ref); gc()


###################
## Format results
###################

# Format GWAS
gwas <- gwas %>%
  # Filter out covariate results
  filter(TEST=="ADD") %>%
  # Join GWAS & freq 
  left_join(., freq, by=c("CHR", "POS","ID", "REF", "ALT"))

# Join pheno and pca
pca <- pca %>%
  left_join(., pheno, by="IID")


#######################
## Plot PCA by batch
#######################

# Create plot
pca %>%
  mutate(batch = as.factor(batch)) %>%
  ggplot(., aes(x=PC1, y=PC2, color=batch)) + 
  geom_point(alpha=.6) + theme_classic()

# Save plot
plot_file_pca <- paste0(dirname(GWAS_FILE), "/",gsub("postgwas_","pca_", 
                        gsub(".batch.glm.logistic.hybrid", ".png", basename(GWAS_FILE), perl = F)))
ggsave(filename=plot_file_pca,  width = 25, height = 15, units = "cm", dpi=300)


#######################
## Plot freq by batch
#######################

# Create plot
#gwas %>%
  #pivot_longer(cols=c("ALT_FREQS_TARGET", "ALT_FREQS_REF"), names_to = "batch", values_to = "ALT_FREQ") %>%
gwas %>%
  arrange(-log10(P)) %>%
  ggplot(., aes(x=ALT_FREQS_TARGET, y=ALT_FREQS_REF, color=-log10(P))) + 
  geom_point(alpha=.6) + scale_color_gradient(low="lightyellow", high="red") + 
  theme_classic()

# Save plot
plot_file_freq <- paste0(dirname(GWAS_FILE), "/",gsub("postgwas_","freq_", 
                                                     gsub(".batch.glm.logistic.hybrid", ".png", basename(GWAS_FILE), perl = F)))
ggsave(filename=plot_file_freq,  width = 25, height = 15, units = "cm", dpi=300)


#######################
## Plot GWAS results
#######################

# Add cumulative base-pair position for plotting
gwas <- gwas %>%
  arrange(CHR, POS) %>%
  mutate(POS_cum = POS)
offset <- 0
chr_order <- unique(gwas$CHR)
chr_order <- chr_order[order(chr_order)]
for (chr in chr_order) {
  idx <- which(gwas$CHR == chr)
  gwas$POS_cum[idx] <- gwas$POS[idx] + offset
  offset <- max(gwas$POS_cum[idx])
}

# Prepare position for X axis labels
axis_df <- gwas %>%
  group_by(CHR) %>%
  summarize(center = (min(POS_cum) + max(POS_cum)) / 2)

# Plot all p values
ggplot(gwas, aes(x = POS_cum, y = -log10(P), color = as.factor(CHR %% 2))) +
  geom_point(alpha = 0.6, size = 0.8) +
  scale_color_manual(values = c("gray20", "gray60")) +
  scale_x_continuous(label = axis_df$CHR, breaks = axis_df$center) +
  geom_hline(yintercept = -log10(5e-8), colour = "red", lty=2) +
  labs(x = "Chromosome", y = "-log10(P)", 
       title = "Manhattan Plot of GWAS on Batch Effect") +
  theme_minimal() +
  theme(legend.position = "none", 
        panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 1))
plot_file <- gsub(".batch.glm.logistic.hybrid", ".png", GWAS_FILE, perl = F)
ggsave(filename=plot_file,  width = 25, height = 15, units = "cm", dpi=300)



#################################
## Save list of significant vars
#################################

# Check number of variants that exceed P value threshold
#sum(gwas$P <= 0.05, na.rm = T)
#custom_p <- 0.05/nrow(gwas)
#sum(gwas$P <= custom_p, na.rm = T)
cat("Number of variants that are significantly associated with batch after adjusting for PC1 (significance level 5e-08):", sep="\n")
sum(gwas$P <= 5e-8, na.rm = T)


# Create list of variants that have significant batch effect (p <= 5e-8)
sig_vars <- gwas %>% 
  filter(P <= 5e-8) %>% select(ID) %>%
  mutate(ID = gsub("\\;.*", "", ID)) %>%
  separate(col="ID", into=c("CHR", "POS", "REF", "ALT")) %>%
  mutate(CHR = gsub("^chr", "", CHR)) %>%
  select("CHR", "POS")

# Save list
sig_vars_file <- gsub("\\.batch\\.glm\\.logistic\\.hybrid$", "_batch_variants.txt", GWAS_FILE)
write_tsv(sig_vars, file=sig_vars_file, col_names = F)


#####################################################
## Update SNP lists w/ overlapping and excluded vars
#####################################################

# Import lists
overlap <- read_tsv(OVERLAP_LIST, col_names=c("chr", "pos"), show_col_types = FALSE)
exclude <-  read_tsv(EXCLUDE_LIST, show_col_types = FALSE)

# Remove batch SNPs from overlap
sig_vars$ID <- paste0(sig_vars$CHR, ":", sig_vars$POS)
overlap <- overlap %>%
  mutate(ID = paste0(chr, ":", pos)) %>%
  filter( ! ID %in% sig_vars$ID ) %>%
  select(-ID)

# Add batch SNPs to exclude
add <- data.frame(chr = sig_vars$CHR, pos = sig_vars$POS,
                  effect_allele=NA, other_allele=NA,
                  AMBIGUOUS_SNP = FALSE, 
                  MISSING_POS_FROM_REF = FALSE, INCONSISTENT_ALLELE_FROM_REF = FALSE,
                  MISSING_POS_FROM_TARGET = FALSE, INCONSISTENT_ALLELE_FROM_TARGET = FALSE,
                  BATCH_EFFECT_AFTER_PC1 = TRUE)
exclude$BATCH_EFFECT_AFTER_PC1 <- FALSE # add new column for batch effect in exclude table
exclude <- rbind(exclude, add)


## End
