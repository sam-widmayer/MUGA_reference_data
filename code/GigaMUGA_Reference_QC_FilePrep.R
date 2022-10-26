#!/usr/bin/env Rscript
require(dplyr)
require(tidyr)
require(stringr)
require(parallel)
require(data.table)
require(multidplyr)
require(purrr)
require(furrr)
require(magrittr)
require(tictoc)
require(furrr)
require(vroom)

## Reading in marker annotations fro Broman, Gatti, & Cornes analysis
gm_metadata <- vroom::vroom("../data/GigaMUGA/gm_uwisc_v2.csv", num_threads = detectCores())
gm_samples <- vroom::vroom("../data/GigaMUGA/Sample_Map.txt", num_threads = detectCores())

# Reading in genotypes
print("time to read in genotypes:")
tictoc::tic()
control_genotypes <- vroom::vroom("../data/GigaMUGA/GigaMUGA_control_genotypes.txt", 
                                  num_threads = detectCores(), progress = T)
tictoc::toc()

# Joining metadata to genotypes
print("time to join metadata to genotypes:")
tictoc::tic()
control_genotypes %<>%
  dplyr::rename(marker = `SNP Name`,
                sample_id = `Sample ID`,
                allele1 = `Allele1 - Forward`,
                allele2 = `Allele2 - Forward`) %>%
  dplyr::left_join(., gm_metadata)
tictoc::toc()

# Nesting genotypes by chromosome
print("time to nest genotypes by chromosome:")
tictoc::tic()
control_genotypes_nest_chr <- control_genotypes %>%
  dplyr::group_by(chr) %>%
  tidyr::nest()
tictoc::toc()



