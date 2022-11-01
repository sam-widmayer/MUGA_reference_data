#!/usr/bin/env Rscript
require(dplyr)
require(fst)
require(tidyr)
require(stringr)
require(parallel)
require(data.table)
require(multidplyr)
require(purrr)
require(furrr)
require(magrittr)
require(vroom)

args <- commandArgs(trailingOnly = TRUE)
callGeno <- function(x){
  
  # Compare genotypes from both strains
      if(c(x[[2]] == x[[3]])){
        #If they are the same, keep the first value
        predicted.geno <- x[[2]]
      } else {
        #If they are different, code as a het
        predicted.geno <- "H"
      }
      return(predicted.geno)
}
writeChrGenos <- function(x,y){
  x$genotype <- apply(x, 1, callGeno)
  fst::write.fst(x, path = paste0("data/GigaMUGA/gm_genos_chr_",y,".fst"))
}

# print(paste("Number of cores: ",detectCores()))

# # Reading in genotypes
# if(!file.exists("data/GigaMUGA/GigaMUGA_control_genotypes.fst")){
#   control_genotypes <- vroom::vroom("data/GigaMUGA/GigaMUGA_control_genotypes.txt",
#                                     num_threads = parallel::detectCores(),
#                                     progress = T)
#   fst::write.fst(control_genotypes, "data/GigaMUGA/GigaMUGA_control_genotypes.fst")
#   }

gm_metadata <- vroom::vroom("data/GigaMUGA/gm_uwisc_v2.csv",
                            progress = T)
print(paste("Writing genotype files for chromosome",args[1]))
control_genotypes_selected_chr <- read.fst("data/GigaMUGA/GigaMUGA_control_genotypes.fst") %>%
  dplyr::rename(marker = `SNP Name`,
                sample_id = `Sample ID`,
                allele1 = `Allele1 - Forward`,
                allele2 = `Allele2 - Forward`) %>%
  dplyr::select(marker, allele1, allele2, everything()) %>%
  dplyr::left_join(., gm_metadata) %>%
  dplyr::filter(chr == args[1])

noNIH <- control_genotypes_selected_chr %>%
  dplyr::mutate(flag = dplyr::if_else(stringr::str_detect(string = sample_id, pattern = "NIH_"),true = "T",false = "F")) %>%
  dplyr::filter(flag == "F") %>%
  dplyr::select(-flag)

writeChrGenos(x = noNIH, y = args[1])
# 
# future::plan(multisession, workers = parallel::detectCores())
# make_chunks <- furrr:::make_chunks
# if(detectCores() > length(control_genotypes_nest_chr$chr)){
#   make_chunks(n_x = length(control_genotypes_nest_chr$chr), n_workers = parallel::detectCores())
# } else {
#   make_chunks(n_x = length(control_genotypes_nest_chr$chr), chunk_size = parallel::detectCores()/8)
# }

# furrr::future_map2(control_genotypes_nest_chr$data,
#                    control_genotypes_nest_chr$chr,
#                    writeChrGenos)
