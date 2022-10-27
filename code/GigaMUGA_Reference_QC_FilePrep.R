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
require(vroom)


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
print(paste0("Writing Chromosome ",y)) 
 vroom::vroom_write(x, 
                    file = paste0("data/GigaMUGA/gm_genos_chr_",y,".csv"),
                    delim = ",",
                    num_threads = parallel::detectCores())
}

## Reading in marker annotations fro Broman, Gatti, & Cornes analysis
gm_metadata <- vroom::vroom("data/GigaMUGA/gm_uwisc_v2.csv", num_threads = parallel::detectCores())


# Reading in genotypes
print("time to read in genotypes:")
tictoc::tic()
control_genotypes <- vroom::vroom("data/GigaMUGA/GigaMUGA_control_genotypes.txt", num_threads = detectCores(), progress = T)
tictoc::toc()



# Joining metadata to genotypes
print("time to join metadata to genotypes:")
tictoc::tic()
control_genotypes %<>%
  dplyr::rename(marker = `SNP Name`,
                sample_id = `Sample ID`,
                allele1 = `Allele1 - Forward`,
                allele2 = `Allele2 - Forward`) %>%
  dplyr::select(marker, allele1, allele2, everything()) %>%
  dplyr::left_join(., gm_metadata)
tictoc::toc()



# Nesting genotypes by chromosome
print("time to nest genotypes by chromosome:")
tictoc::tic()
control_genotypes_nest_chr <- control_genotypes %>%
  dplyr::group_by(chr) %>%
  tidyr::nest()
tictoc::toc()


print("Writing chromosome-level genotype files")

future::plan(multicore)
furrr::future_map2(control_genotypes_nest_chr$data,
                   control_genotypes_nest_chr$chr,
                   writeChrGenos)
