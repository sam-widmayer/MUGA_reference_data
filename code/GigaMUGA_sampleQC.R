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
suppressMessages(make_chunks(n_x = 23, chunk_size = 1))

# Calculating frequency of each genotype for each marker
control_allele_freqs <- furrr::future_map(control_genotype_files,function(x){
                    chr_geno <- read.fst(paste0("data/GigaMUGA/",x))
                    control_allele_freqs <- chr_geno %>%
                      dplyr::group_by(marker, genotype) %>%
                      dplyr::count() %>% 
                      dplyr::ungroup() %>%
                      dplyr::group_by(marker) %>%
                      dplyr::mutate(genotype = if_else(condition = genotype == "-", 
                                                       true = "N", 
                                                       false = as.character(genotype)),
                                    freq = round(n/sum(n), 3),
                                    genotype = as.factor(genotype))
                    return(control_allele_freqs)})
control_allele_freqs_df <- Reduce(dplyr::bind_rows, control_allele_freqs)

## Filtering to markers with missing genotypes
no.calls <- control_allele_freqs_df %>%
  dplyr::ungroup() %>%
  dplyr::filter(genotype == "N") %>%
  tidyr::pivot_wider(names_from = genotype, 
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
                                         dplyr::group_by(sample_id, genotype) %>%
                                         dplyr::count() %>% 
                                         dplyr::ungroup() %>%
                                         dplyr::group_by(sample_id) %>%
                                         dplyr::mutate(genotype = if_else(condition = genotype == "-", 
                                                                         true = "N", 
                                                                         false = as.character(genotype))) %>%
                                         dplyr::filter(genotype == "N")
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

bad_sample_cutoff <- quantile(n.calls.strains.df$n.no.calls, probs = seq(0,1,0.05))[20]
high.n.samples <- n.calls.strains.df %>%
  dplyr::filter(n.no.calls > bad_sample_cutoff)



## Using regex searching of sample IDs to deduce the sex of each sample

###########################################################
## Key processes and expected outputs at each iteration: ##
###########################################################

#####################################################################
## 1) Extracting the a substring of X digits into the sample name. ##
#####################################################################

#####################################################################################
## 2) Assigning the predicted sex based on expected mouse nomenclature convention. ##
#####################################################################################


##############################################################################################
# 3) Inferring the strain background by removing the mouse id from the sample name when a sex is predicted. In certain cases, symbols had to be extracted prior to sex and background inference.
##############################################################################################


predicted.sexes <- n.calls.strains.df %>% 
  dplyr::select(sample_id) %>%
                # One character
  dplyr::mutate(mouse.id.1 = stringr::str_sub(sample_id, -1),
                predicted.sex.1 = dplyr::case_when(stringr::str_sub(mouse.id.1, 1, 1) %in% c("m","M") ~ "m",
                                                   stringr::str_sub(mouse.id.1, 1, 1) %in% c("f","F") ~ "f",
                                                 TRUE ~ "unknown"),
                # Two characters
                mouse.id.2 = stringr::str_sub(sample_id, -2),
                predicted.sex.2 = dplyr::case_when(stringr::str_sub(mouse.id.2, 1, 1) %in% c("m","M") ~ "m",
                                                   stringr::str_sub(mouse.id.2, 1, 1) %in% c("f","F") ~ "f",
                                                 TRUE ~ "unknown"),
                # Three characters
                mouse.id.3= stringr::str_sub(sample_id, -3),
                predicted.sex.3 = dplyr::case_when(stringr::str_sub(mouse.id.3, 1, 1) %in% c("m","M") ~ "m",
                                                 stringr::str_sub(mouse.id.3, 1, 1) %in% c("f","F") ~ "f",
                                                 TRUE ~ "unknown"),
                # Four characters
                mouse.id.4 = stringr::str_sub(sample_id, -4),
                predicted.sex.4 = dplyr::case_when(stringr::str_sub(mouse.id.4, 1, 1) %in% c("m","M") ~ "m",
                                                 stringr::str_sub(mouse.id.4, 1, 1) %in% c("f","F") ~ "f",
                                                 TRUE ~ "unknown"),
                # Five characters
                mouse.id.5 = stringr::str_sub(sample_id, -5),
                mouse.id.5 = stringr::str_replace(string = mouse.id.5,  ## a couple symbols in these ids mess up the regex search
                                                  pattern = "[:symbol:]", 
                                                  replacement = ""),
                predicted.sex.5 = dplyr::case_when(stringr::str_sub(mouse.id.5, 1, 1) %in% c("m","M") ~ "m",
                                                 stringr::str_sub(mouse.id.5, 1, 1) %in% c("f","F") ~ "f",
                                                 TRUE ~ "unknown"),
                # Six characters
                mouse.id.6 = stringr::str_sub(sample_id, -6),
                mouse.id.6 = stringr::str_replace(string = mouse.id.6,  ## a couple symbols in these ids mess up the regex search
                                                  pattern = "[:punct:]", 
                                                  replacement = ""),
                mouse.id.6 = stringr::str_replace(string = mouse.id.6,  ## a couple symbols in these ids mess up the regex search
                                                  pattern = "[:symbol:]", 
                                                  replacement = ""),
                predicted.sex.6 = dplyr::case_when(stringr::str_sub(mouse.id.6, 1, 1) %in% c("m","M") ~ "m",
                                                 stringr::str_sub(mouse.id.6, 1, 1) %in% c("f","F") ~ "f",
                                                 TRUE ~ "unknown"),
                # Seven characters
                mouse.id.7 = stringr::str_sub(sample_id, -7),
                mouse.id.7 = stringr::str_replace(string = mouse.id.7,  ## a couple symbols in these ids mess up the regex search
                                                  pattern = "[:punct:]", 
                                                  replacement = ""),
                mouse.id.7 = stringr::str_replace(string = mouse.id.7,  ## a couple symbols in these ids mess up the regex search
                                                  pattern = "[:symbol:]", 
                                                  replacement = ""),
                predicted.sex.7 = dplyr::case_when(stringr::str_sub(mouse.id.7, 1, 1) %in% c("m","M") ~ "m",
                                                 stringr::str_sub(mouse.id.7, 1, 1) %in% c("f","F") ~ "f",
                                                 TRUE ~ "unknown"),
                # Eight characters
                mouse.id.8 = stringr::str_sub(sample_id, -8),
                mouse.id.8 = stringr::str_replace(string = mouse.id.8,  ## a couple symbols in these ids mess up the regex search
                                                  pattern = "[:punct:]", 
                                                  replacement = ""),
                mouse.id.8 = stringr::str_replace(string = mouse.id.8,  ## a couple symbols in these ids mess up the regex search
                                                  pattern = "[:symbol:]", 
                                                  replacement = ""),
                predicted.sex.8 = dplyr::case_when(stringr::str_sub(mouse.id.8, 1, 1) %in% c("m","M") ~ "m",
                                                 stringr::str_sub(mouse.id.8, 1, 1) %in% c("f","F") ~ "f",
                                                 TRUE ~ "unknown"),
                predicted.sex.CAP = dplyr::case_when(stringr::str_detect(string = sample_id, pattern = "_M_") ~ "m",
                                                     stringr::str_detect(string = sample_id, pattern = "_F_") ~ "f",
                                                     TRUE ~ "unknown")) %>%
  dplyr::mutate(predicted.sex = dplyr::case_when(predicted.sex.1 == "m" ~ "m", 
                                                 predicted.sex.2 == "m" ~ "m", 
                                                 predicted.sex.3 == "m" ~ "m",
                                                 predicted.sex.4 == "m" ~ "m",
                                                 predicted.sex.5 == "m" ~ "m",
                                                 predicted.sex.6 == "m" ~ "m",
                                                 predicted.sex.7 == "m" ~ "m",
                                                 predicted.sex.8 == "m" ~ "m",
                                                 predicted.sex.CAP == "m" ~ "m",
                                                 
                                                 predicted.sex.1 == "f" ~ "f",
                                                 predicted.sex.2 == "f" ~ "f",
                                                 predicted.sex.3 == "f" ~ "f",
                                                 predicted.sex.4 == "f" ~ "f",
                                                 predicted.sex.5 == "f" ~ "f",
                                                 predicted.sex.6 == "f" ~ "f",
                                                 predicted.sex.7 == "f" ~ "f",
                                                 predicted.sex.8 == "f" ~ "f",
                                                 predicted.sex.CAP == "f" ~ "f",
                                                 TRUE ~ "unknown")) %>%
  dplyr::distinct(sample_id, predicted.sex)

## Reading in probe intensities
X.fst <- read.fst("data/GigaMUGA/gm_genos_chr_X.fst")
Y.fst <- read.fst("data/GigaMUGA/gm_genos_chr_Y.fst")
long_XY_intensities <- X.fst %>%
    dplyr::bind_rows(.,Y.fst) %>%
    dplyr::left_join(., gm_metadata)

save(control_allele_freqs_df, above.cutoff,
    n.calls.strains.df, high.n.samples,
    predicted.sexes, long_XY_intensities, file = "data/GigaMUGA/Marker_QC.RData")
