#!/usr/bin/env Rscript

################################################################################
# Deriving consensus Chromosome Y and Mitochondrial haplotypes using GigaMUGA
# reference genotype data

# Inputs:
# 1. CC/DO founder & F1 genotype data, including replicates.
# 2. GigaMUGA marker metadata from Karl Broman, Dan Gatti, and Belinda Cornes
# 3. Output from initial GigaMUGA marker QC

# Sam Widmayer
# samuel.widmayer@jax.org

# 20240214
################################################################################

library(dplyr)
library(tidyr)
library(stringr)
library(fst)
setwd('/projects/compsci/vmp/USERS/widmas/MUGA_reference_data')

### SETUP ###

# Set sample genotype paths
sample_geno_dir <- 'data/GigaMUGA/GigaMUGA_reference_genotypes'

# Set output paths
out_dir <- 'output/GigaMUGA/'

# Read Chr Y Genotypes
y_genos <- fst::read.fst(path = file.path(sample_geno_dir, 'gm_genos_chr_Y.fst'))

# Read Mt Genotypes
mt_genos <- fst::read.fst(path = file.path(sample_geno_dir, 'gm_genos_chr_M.fst'))

# Read marker metadata
gm_metadata <- read.csv(file = 'data/GigaMUGA/gm_uwisc_v4.csv')

# Read sample metadata
load(file = 'data/GigaMUGA/GigaMUGA_QC_Results.RData')
# 'founderSamples' = information for samples derived from CC/DO founders




### MAIN ###

# Removing no-calls and heterozygous sites from Chr Y genotype df
y_wh <- y_genos$sample_id %in% founderSamples$sample_id & (!y_genos$genotype %in% c("-","H"))
y_genos <- y_genos[y_wh,]

# Removing no-calls and heterozygous sites from Chr Mt genotype df
mt_wh <- mt_genos$sample_id %in% founderSamples$sample_id & (!mt_genos$genotype %in% c("-","H"))
mt_genos <- mt_genos[mt_wh,]

# Obtain founder letter code from sample id
y_genos$bg <- stringr::str_sub(string = y_genos$sample_id, 
                               start = 0, end = 2)

mt_genos$bg <- stringr::str_sub(string = mt_genos$sample_id, 
                               start = 0, end = 2)

# Generate letter codes for inbred strains
inbred_strains <- paste0(LETTERS[1:8],LETTERS[1:8])

#############
### CHR Y ###
#############

# Filter Chr Y genotypes to samples from inbred strains
y_genos %<>%
  dplyr::select(bg, marker, genotype) %>%
  dplyr::filter(bg %in% inbred_strains) %>%
  dplyr::arrange(marker)

consensus_y_genos <- y_genos %>%
  dplyr::arrange(bg) %>% 
  dplyr::distinct(bg, genotype, marker) %>%
  tidyr::pivot_wider(names_from = bg, values_from = genotype)

y_genos_consensus <- consensus_y_genos[rowSums(is.na(consensus_y_genos)) < 3,]
colnames(y_genos_consensus) <- c("marker",LETTERS[1:8])


##########
### MT ###
##########
# Filter Chr Mt genotypes to samples from inbred strains
mt_genos %<>%
  dplyr::select(bg, marker, genotype) %>%
  dplyr::filter(bg %in% inbred_strains) %>%
  dplyr::arrange(marker)

# Nest Chr Mt marker genotypes
mt_genos_nested <- mt_genos %>%
  dplyr::group_by(marker) %>%
  tidyr::nest()

# Identify segregating markers
seg <- lapply(mt_genos_nested$data, function(x) nrow(unique(x[,2]))) == 2
mt_genos <- mt_genos[which(mt_genos$marker %in% mt_genos_nested$marker[seg]),]

# Filter Chr Mt genotypes to segregating SNPs and cast unique alleles as consensus genotypes
mt_genos %<>%
  dplyr::distinct(bg, marker, genotype) %>%
  dplyr::arrange(bg) %>%
  tidyr::pivot_wider(names_from = bg, values_from = genotype)

mt_genos_consensus <- mt_genos[complete.cases(mt_genos),]
colnames(mt_genos_consensus) <- c("marker",LETTERS[1:8])


# Write consensus files
write.csv(y_genos_consensus, file = file.path(out_dir,'GigaMUGA_founder_consensus_genotypes_chrY.csv'))
write.csv(mt_genos_consensus, file = file.path(out_dir,'GigaMUGA_founder_consensus_genotypes_chrM.csv'))

