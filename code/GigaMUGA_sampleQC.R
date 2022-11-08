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
control_genotype_files <- list.files("data/GigaMUGA/GigaMUGA_reference_genotypes/", pattern = "gm_genos_chr")

# Creating parallelization plan for furrr
plan(multisession, workers = 23)
make_chunks <- furrr:::make_chunks
suppressMessages(make_chunks(n_x = 23, chunk_size = 1))

# Calculating frequency of each genotype for each marker
control_allele_freqs <- furrr::future_map(control_genotype_files,function(x){
                    chr_geno <- read.fst(paste0("data/GigaMUGA/GigaMUGA_reference_genotypes/",x))
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
                                       chr_geno <- read.fst(paste0("data/GigaMUGA/GigaMUGA_reference_genotypes/",x))
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
X.fst <- read.fst("data/GigaMUGA/GigaMUGA_reference_genotypes/gm_genos_chr_X.fst")
Y.fst <- read.fst("data/GigaMUGA/GigaMUGA_reference_genotypes/gm_genos_chr_Y.fst")
long_XY_intensities <- X.fst %>%
    dplyr::bind_rows(.,Y.fst) %>%
    dplyr::left_join(., gm_metadata)

## Flagging markers and samples based on previous QC steps
flagged_XY_intensities <- long_XY_intensities %>%
  dplyr::mutate(marker_flag = dplyr::if_else(condition = marker %in% above.cutoff$marker,
                                             true = "FLAG",
                                             false = "")) %>%
  dplyr::mutate(high_missing_sample = dplyr::if_else(condition = sample_id %in% high.n.samples$sample_id,
                                                     true = "FLAG",
                                                     false = ""))

# Input: Sex chromosome probe intensities for each marker with 1) marker metdata, 2) marker and sample flags, 3) background and sex predictions
Xchr.int <- flagged_XY_intensities %>%
  dplyr::ungroup() %>%
  dplyr::filter(marker_flag != "FLAG",
                chr == "X") %>%
  dplyr::mutate(x.chr.int = X + Y,
                genotype = dplyr::if_else(genotype == "-", "N", as.character(genotype))) %>%
  dplyr::filter(!is.na(x.chr.int)) %>%
  dplyr::group_by(sample_id, high_missing_sample) %>%
  dplyr::summarise(mean.x.chr.int = mean(x.chr.int))


Ychr.int <- flagged_XY_intensities %>%
  dplyr::ungroup() %>%
  dplyr::filter(chr == "Y") %>% # every marker on the Y chromosome was a flagged marker...so for the purposes of sexing I tried using the bad markers anyways
  dplyr::mutate(genotype = dplyr::if_else(genotype == "-", "N", as.character(genotype))) %>%
  dplyr::filter(!is.na(Y)) %>%
  dplyr::group_by(sample_id, high_missing_sample) %>%
  dplyr::summarise(mean.y.int = mean(Y))
# Expected output: Sample-averaged y-channel probe intensities for all chromosome Y markers. Note: replicated sample information collapses at this step. This is tolerable under the assumption that the samples with identical names are in fact duplicates of the same individual.



# Column binding the two intensities if the sample information matches
sex.chr.intensities <- Xchr.int %>%
  dplyr::full_join(., Ychr.int) %>%
  dplyr::left_join(., predicted.sexes) 
  # dplyr::filter(high_missing_sample == "")

# Clear visual clustering of samples motivated us to use a rough clustering method to quickly assign groups to samples based on X and Y chromsome probe intensities. K-means clustering is below supplying two clusters for each sex.
# Inputs: 
# 1) Sample-averaged summed x- and y-channel probe intensities for all chromosome X markers
# 2) Sample-averaged y-channel probe intensities for all chromosome Y markers
sex.chr.intensities.goodsamples <- sex.chr.intensities %>%
  dplyr::filter(high_missing_sample == "")
kmeans.x <- sex.chr.intensities.goodsamples %>%
  dplyr::ungroup() %>%
  dplyr::select(mean.x.chr.int) %>%
  dplyr::filter(!is.na(mean.x.chr.int)) %>%
  kmeans(., centers = 2)
kmeans.y <- sex.chr.intensities.goodsamples %>%
  dplyr::ungroup() %>%
  dplyr::select(mean.y.int) %>%
  dplyr::filter(!is.na(mean.y.int)) %>%
  kmeans(., centers = 2)

# Joining each sample's cluster assignment to the sample-averaged intensity metrics
sex.chr.k.means.x <- cbind(sex.chr.intensities.goodsamples %>% dplyr::filter(!is.na(mean.x.chr.int)), 
      kmeans.x$cluster)
colnames(sex.chr.k.means.x) <- c(colnames(sex.chr.k.means.x)[-6],"x.clust") 

sex.chr.k.means.y <- cbind(sex.chr.intensities.goodsamples %>% dplyr::filter(!is.na(mean.y.int)), 
      kmeans.y$cluster)
colnames(sex.chr.k.means.y) <- c(colnames(sex.chr.k.means.y)[-6],"y.clust")

# Generating a contingency table for how each cluster paired with each sex. 
sex.by.cluster.tab.x <- sex.chr.k.means.x %>%
  dplyr::group_by(predicted.sex, x.clust) %>% 
  dplyr::count() %>%
  dplyr::arrange(desc(n))
sex.by.cluster.tab.y <- sex.chr.k.means.y %>%
  dplyr::group_by(predicted.sex, y.clust) %>% 
  dplyr::count() %>%
  dplyr::arrange(desc(n))

# The most common clusters should be the two sexes, k-means doesn't always assign the same cluster name to the same sex. Therefore, the top clusters must be pulled out and assigned sexes dynamically.
top.clusters.x <- sex.by.cluster.tab.x[1:2,] %>%
  dplyr::ungroup() %>%
  dplyr::mutate(inferred.sex = predicted.sex) %>%
  dplyr::select(-n,-predicted.sex)
top.clusters.y <- sex.by.cluster.tab.y[1:2,] %>%
  dplyr::ungroup() %>%
  dplyr::mutate(inferred.sex = predicted.sex) %>%
  dplyr::select(-n,-predicted.sex)

# Samples are then recoded according to the k-means assigned sexes
reSexed.x <- sex.chr.k.means.x %>%
  dplyr::select(-predicted.sex) %>%
  dplyr::left_join(., top.clusters.x)
reSexed.y <- sex.chr.k.means.y %>%
  dplyr::select(-predicted.sex) %>%
  dplyr::left_join(., top.clusters.y)
reSexed_samples <- full_join(reSexed.x, reSexed.y) %>%
  dplyr::full_join(sex.chr.intensities.goodsamples %>%
                     dplyr::select(sample_id, predicted.sex),.) %>% 
  dplyr::select(-x.clust, -y.clust) %>%
  dplyr::distinct()

# Vector of founder strain names
founder_strains <- c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ",
                     "NZO/HILtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ")

reSexed_samples$rough_founder_pull <- strsplit(reSexed_samples$sample_id, split =  "_") %>%
  purrr::map(., function(x){return(x[[1]])}) %>% 
  unlist()

founderSamples <- reSexed_samples %>%
  dplyr::mutate(founder = dplyr::if_else(stringr::str_length(rough_founder_pull) == 2, true = "founder", false = "")) %>%
  dplyr::filter(founder == "founder") %>%
  dplyr::mutate(dam = str_sub(rough_founder_pull, start = 1, end = 1),
                sire = str_sub(rough_founder_pull, start = 2, end = 2),
                dam = dplyr::case_when(dam == "A" ~ founder_strains[1],
                                       dam == "B" ~ founder_strains[2],
                                       dam == "C" ~ founder_strains[3],
                                       dam == "D" ~ founder_strains[4],
                                       dam == "E" ~ founder_strains[5],
                                       dam == "F" ~ founder_strains[6],
                                       dam == "G" ~ founder_strains[7],
                                       dam == "H" ~ founder_strains[8]),
                sire = dplyr::case_when(sire == "A" ~ founder_strains[1],
                                        sire == "B" ~ founder_strains[2],
                                        sire == "C" ~ founder_strains[3],
                                        sire == "D" ~ founder_strains[4],
                                        sire == "E" ~ founder_strains[5],
                                        sire == "F" ~ founder_strains[6],
                                        sire == "G" ~ founder_strains[7],
                                        sire == "H" ~ founder_strains[8]),
                bg = dplyr::if_else(dam == sire, true = "INBRED", false = "CROSS")) %>%
  dplyr::select(sample_id, inferred.sex, dam, sire, bg, predicted.sex)


save(control_allele_freqs_df,
     n.calls.strains.df,
    founderSamples, file = "data/GigaMUGA/GigaMUGA_QC_Results.RData")
save(above.cutoff, high.n.samples, file = "data/GigaMUGA/GigaMUGA_BadSamples_BadMarkers.RData")
save(predicted.sexes, sex.chr.intensities, reSexed_samples, file = "data/GigaMUGA/GigaMUGA_SexCheck_Results.RData")
