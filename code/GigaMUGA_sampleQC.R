#!/usr/bin/env Rscript
require(dplyr)
require(fst)
require(tidyr)
require(stringr)
require(parallel)
require(data.table)
require(purrr)
require(furrr)
require(magrittr)
require(vroom)

# Reading in source data
gm_metadata <- vroom::vroom("data/GigaMUGA/gm_uwisc_v2.csv", num_threads = detectCores())
control_genotype_files <- list.files("data/GigaMUGA/", pattern = "gm_genos_chr")

# Creating parallelization plan for furrr
plan(multisession, workers = 23)
make_chunks <- furrr:::make_chunks
make_chunks(n_x = 23, chunk_size = 1)

# Calculating frequency of each genotype for each marker
control_allele_freqs <- furrr::future_map(control_genotype_files,function(x){
                    chr_geno <- read.fst(paste0("data/GigaMUGA/",x))
                    control_allele_freqs <- chr_geno %>%
                      dplyr::group_by(marker, allele1) %>%
                      dplyr::count() %>% 
                      dplyr::ungroup() %>%
                      dplyr::group_by(marker) %>%
                      dplyr::mutate(allele1 = if_else(condition = allele1 == "-", 
                                                       true = "N", 
                                                       false = as.character(allele1)),
                                    freq = round(n/sum(n), 3),
                                    allele1 = as.factor(allele1))
                    return(control_allele_freqs)})
control_allele_freqs_df <- Reduce(dplyr::bind_rows, control_allele_freqs)

## Filtering to frequencies of missing genotypes
no.calls <- control_allele_freqs_df %>%
  dplyr::ungroup() %>%
  dplyr::filter(allele1 == "N") %>%
  tidyr::pivot_wider(names_from = allele1, 
                     values_from = n) %>%
  dplyr::left_join(., gm_metadata) %>%
  dplyr::select(marker, chr, bp_grcm39, freq) %>%
  dplyr::mutate(chr = as.factor(chr))

## Identifying markers with missing genotypes at a frequency higher than the 95th percentile of "N" frequencies across all markers
cutoff <- quantile(no.calls$freq, probs = seq(0,1,0.05))[[19]]
above.cutoff <- no.calls %>%
  dplyr::filter(freq > cutoff)





## Calculating the number of missing markers for each sample
n.calls.strains <- furrr::future_map(control_genotype_files,
                                     function(x){
                                       chr_geno <- read.fst(paste0("data/GigaMUGA/",x))
                                       sample_Ns <- chr_geno %>%
                                         dplyr::group_by(sample_id, allele1) %>%
                                         dplyr::count() %>% 
                                         dplyr::ungroup() %>%
                                         dplyr::group_by(sample_id) %>%
                                         dplyr::mutate(allele1 = if_else(condition = allele1 == "-", 
                                                                         true = "N", 
                                                                         false = as.character(allele1))) %>%
                                         dplyr::filter(allele1 == "N")
                                       return(sample_Ns)}) %>% 
  Reduce(dplyr::bind_rows, .)
  
n.calls.strains.df <- n.calls.strains %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(n.no.calls = sum(n)) %>%
  dplyr::select(n.no.calls, sample_id)
#   n.no.calls sample_id      
#        <int> <chr>          
# 1       8088 10570m4381     
# 2       8322 11531m8014     
# 3       8314 11666m45109    
# 4       8158 11679m4874     
# 5       8034 11888m7081     
# 6       8420 129P1/ReJm35858

save(control_allele_freqs_df, n.calls.strains.df, file = "data/GigaMUGA/Marker_QC.RData")
