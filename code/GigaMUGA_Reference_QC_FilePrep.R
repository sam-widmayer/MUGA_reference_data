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
  x$genotype <- apply(x, 1, callGeno)
  fst::write.fst(x, path = paste0("data/GigaMUGA/gm_genos_chr_",y,".fst"))
}

print(paste("Number of cores: ",detectCores()))

# Reading in genotypes
print("time to read in genotypes and write fst:")
tictoc::tic()
control_genotypes <- vroom::vroom("data/GigaMUGA/GigaMUGA_control_genotypes.txt",
                                  num_threads = parallel::detectCores(),
                                  progress = T)
fst::write.fst(control_genotypes, "data/GigaMUGA/GigaMUGA_control_genotypes.fst")
tictoc::toc()



tictoc::tic()
control_genotypes_nest_chr <- read.fst("data/GigaMUGA/GigaMUGA_control_genotypes.fst") %>%
  dplyr::rename(marker = `SNP Name`,
                sample_id = `Sample ID`,
                allele1 = `Allele1 - Forward`,
                allele2 = `Allele2 - Forward`) %>%
  dplyr::select(marker, allele1, allele2, everything()) %>%
  dplyr::group_by(chr) %>%
  tidyr::nest()
tictoc::tic()


print("Writing chromosome-level genotype files")

future::plan(multisession, workers = detectCores())
make_chunks <- furrr:::make_chunks
if(detectCores() > length(control_genotypes_nest_chr$chr)){
  make_chunks(n_x = length(control_genotypes_nest_chr$chr), n_workers = detectCores())
} else {
  make_chunks(n_x = length(control_genotypes_nest_chr$chr), chunk_size = detectCores()/8)
}

furrr::future_map2(control_genotypes_nest_chr$data,
                   control_genotypes_nest_chr$chr,
                   writeChrGenos)
