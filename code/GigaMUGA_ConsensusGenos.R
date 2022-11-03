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

###########################################
# Function used in the loop which forms founder consensus genotypes prior to background QC
# Inputs:
# mk = marker
# data = genotype data for each founder sample
# f = founder strain
# Outputs:
# data frame with 1 row with columns: marker; consensus genotype for founder strain f 
###########################################
removeExtremeInts <- function(mk, data, f){
    
    # Calculate summary statistics for probe intensities
    sd.x.int <- sd(data$X)
    sd.y.int <- sd(data$Y)
    mean.x.int <- mean(data$X)
    mean.y.int <- mean(data$Y)

    data %>%
      # With a prior expectation provided by an "N" call, determine whether probe intensities are unusual and 
      # flag samples that meet both criteria
      dplyr::mutate(x.ex = dplyr::if_else(X > (mean.x.int + sd.x.int) | X < (mean.x.int - sd.x.int), true = "EX", false = ""),
                    y.ex = dplyr::if_else(Y > (mean.y.int + sd.y.int) | Y < (mean.y.int - sd.y.int), true = "EX", false = ""),
                    flag = dplyr::if_else(genotype == "N" & (x.ex == "EX" | y.ex == "EX"), true = "FLAG", false = "")) %>%
      dplyr::filter(flag != "FLAG") %>%
      dplyr::mutate(marker = mk) %>%
      dplyr::distinct(marker,genotype) %>%
      dplyr::select(marker,genotype) %>%
      `colnames<-`(c("marker",f))
  }

findConsensusGenotypes <- function(mk, data){
  
  # Examples use mk = "UNC13666424"
  
  # Annotate sample genotype data with summed intensities as a heuristic and filter to sample genotypes assigned a real genotype
  call_filtered_data <- data %>%
    dplyr::mutate(sum_int = X+Y) %>%
    dplyr::filter(sample_genotype %in% c("A","C","T","G","H"))
  
  # Derive a lower bound by which to identify true N calls using the intensity data
  real_N_cutoff <- quantile(call_filtered_data$sum_int, probs = seq(0,1,0.05))[[2]]
  # Use this threshold to identify potential miscalls and remove them in order to recode samples without consensus calls
  filtered_data <- call_filtered_data %>%
    dplyr::mutate(test = if_else(sum_int < real_N_cutoff & is.na(consensus_genotype), true = "miscall", false = "")) %>%
    dplyr::filter(test != "miscall")
  
  # Determine the number of clusters to use for k-means clustering of intensity values;
  # If the previous steps have succeeded, there should be one or two real genotypes segregating and samples can be re-called by k-means
  n_genos <- length(unique(filtered_data$sample_genotype[which(filtered_data$sample_genotype %in% c("A","C","T","G","H"))]))
  ints_for_kmeans <- filtered_data %>%
    dplyr::select(X, Y)
  consensus_geno_clusters <- kmeans(ints_for_kmeans, centers = n_genos)$cluster
  
  # Joining each sample's cluster assignment to the sample intensity metrics
  genos.k.means <- filtered_data %>%
    dplyr::mutate(clust = as.factor(consensus_geno_clusters))
  
  # Create a join table to pair up clusters with genotypes
  gens_by_cluster_tab <- genos.k.means %>%
    dplyr::group_by(sample_genotype, consensus_genotype, clust) %>%
    dplyr::count() %>%
    arrange(-n)
# Example:
#   sample_genotype consensus_genotype clust     n
# 1 G               G                  2        32
# 2 A               A                  1        23
# 3 A               NA                 1         6
  
  # Filter the join table to hopefully have a 1:1 relationship between consensus genotypes and clusters
  top_clusters <- gens_by_cluster_tab %>%
    dplyr::filter(!is.na(consensus_genotype)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(sample_genotype, clust) %>%
    dplyr::rename(consensus_genotype = sample_genotype)
  # Example:
#   consensus_genotype clust
# 1 G                  2    
# 2 A                  1    

  # Samples are then recoded according to the k-means assigned genotypes and noted whether a consensus was re-assigned
  unknown_before <- filtered_data %>%
    dplyr::filter(is.na(consensus_genotype))
  recoded_consensus_genotypes <- genos.k.means %>%
    dplyr::select(-consensus_genotype) %>%
    dplyr::left_join(.,top_clusters, by = "clust") %>%
    dplyr::mutate(marker = mk) %>% 
    dplyr::mutate(recoded = dplyr::if_else(sample_id %in% unknown_before$sample_id, true = "RECODED", false = "")) %>%
    dplyr::select(-clust)
  # Example: recoded_consensus_genotypes[20:30,]
  # Note NOD samples are now recoded with the consensus genotype predicted from k-means across all founder samples
#      strain      sample              x_int y_int sample_genotype sum_int test  consensus_genotype marker      recoded  
#    <chr>       <chr>               <dbl> <dbl> <chr>             <dbl> <chr> <chr>              <chr>       <chr>    
#  1 129S1/SvImJ X129S1.SvImJf0827.1 0.716 0.083 A                 0.799 ""    A                  UNC13666424 ""       
#  2 129S1/SvImJ X129S1.SvImJm       0.618 0.079 A                 0.697 ""    A                  UNC13666424 ""       
#  3 129S1/SvImJ X129S1.SvImJm.1     0.599 0.084 A                 0.683 ""    A                  UNC13666424 ""       
#  4 129S1/SvImJ X129S1.SvImJm1314   0.69  0.097 A                 0.787 ""    A                  UNC13666424 ""       
#  5 NOD/ShiLtJ  NOD.ShiLtJf0713.1   0.567 0.07  A                 0.637 ""    A                  UNC13666424 "RECODED"
#  6 NOD/ShiLtJ  NOD.ShiLtJf0713.2   0.6   0.095 A                 0.695 ""    A                  UNC13666424 "RECODED"
#  7 NOD/ShiLtJ  NOD.ShiLtJf0713.3   0.883 0.091 A                 0.974 ""    A                  UNC13666424 "RECODED"
#  8 NOD/ShiLtJ  NOD.ShiLtJm0150     0.646 0.088 A                 0.734 ""    A                  UNC13666424 "RECODED"
#  9 NOD/ShiLtJ  NOD.ShiLtJm35324    0.538 0.077 A                 0.615 ""    A                  UNC13666424 "RECODED"
# 10 NOD/ShiLtJ  NOD.ShiLtJm39173    0.649 0.087 A                 0.736 ""    A                  UNC13666424 "RECODED"
# 11 NZO/HILtJ   NZO.HILtJf0588      0.041 0.728 G                 0.769 ""    G                  UNC13666424 ""       
  
  
  # Generate a wide table of consensus as expected for output file
  recoded <- recoded_consensus_genotypes %>%
    dplyr::distinct(marker, strain, consensus_genotype) %>% 
    tidyr::pivot_wider(names_from = strain, values_from = consensus_genotype)
#     marker      `A/J` `C57BL/6J` `129S1/SvImJ` `NOD/ShiLtJ` `NZO/HILtJ` `CAST/EiJ` `PWK/PhJ` `WSB/EiJ`
#   <chr>       <chr> <chr>      <chr>         <chr>        <chr>       <chr>      <chr>     <chr>    
# 1 UNC13666424 G     G          A             A            G           A          A         G     
  
  # Check for any founders that still segregate a genotype within the marker
  still_consensus <- recoded_consensus_genotypes %>%
    dplyr::group_by(strain, consensus_genotype) %>%
    dplyr::count() %>%
    dplyr::group_by(strain) %>%
    dplyr::count()
  
  # If there are still discrepancies (i.e. non-N genotypes segregating among samples), return the equivalent of the input data for downstream work
  if(2 %in% still_consensus$n){
    original_consensus <- suppressWarnings(data %>%
                                     dplyr::mutate(marker = mk) %>%
                                     dplyr::distinct(marker, strain, consensus_genotype) %>%
                                     tidyr::pivot_wider(names_from = strain, values_from = consensus_genotype))
    
    original_sample_calls <- suppressWarnings(data %>%
                                                dplyr::mutate(marker = mk) %>%
                                                dplyr::distinct(marker, sample, sample_genotype) %>%
                                                tidyr::pivot_wider(names_from = sample,
                                                                   values_from = sample_genotype))
    original_data <- data %>%
      dplyr::mutate(marker = mk,
                    recoded = "")
    
    return(list(original_consensus, original_sample_calls, original_data))
  } else{
    sample_new_consensus_calls <- recoded_consensus_genotypes %>%
                  dplyr::distinct(marker, sample, consensus_genotype) %>% 
                  tidyr::pivot_wider(names_from = sample, 
                                     values_from = consensus_genotype)
    return(list(recoded, sample_new_consensus_calls, recoded_consensus_genotypes))
  }
}

#args <- commandArgs(trailingOnly = TRUE)
load("data/GigaMUGA/GigaMUGA_QC_Results.RData")
load("data/GigaMUGA/bad_samples_markers.RData")
load("data/GigaMUGA/sex_check_results.RData")
gm_metadata <- vroom::vroom("data/GigaMUGA/gm_uwisc_v2.csv",
                            progress = T)

good_chr_genos <- fst::read.fst(paste0("data/GigaMUGA/gm_genos_chr_",args[1],".fst")) %>%
            dplyr::mutate(marker_flag = dplyr::if_else(condition = marker %in% above.cutoff$marker, 
                            true = "FLAG", 
                            false = ""),
                         genotype = dplyr::if_else(condition = genotype == "-", true = "N", false = genotype)) %>%
            dplyr::filter(marker_flag == "") %>%
            dplyr::select(-marker_flag)

founder_strains <- c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ",
                     "NZO/HILtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ")


for(f in founder_strains){
    print(paste("Generating Calls for",f))
  
  # Pulling the samples and genotypes for each CC/DO founder strain
  founder.geno.array <- founderSamples %>%
                dplyr::filter(dam == f,
                            sire == f) %>%
    # Attach all genotypes
    dplyr::left_join(.,good_chr_genos, by = "sample_id")


  # Count the number of unique allele calls for each marker
  founder.allele.counts <- founder.geno.array %>%
    dplyr::group_by(marker, genotype) %>%
    dplyr::count()
  
  # Collect markers which have identical genotypes across samples from the same founder
  complete_founder_genos <- founder.allele.counts %>%
    dplyr::filter(n == max(founder.allele.counts$n)) %>%
    dplyr::select(-n) %>%
    `colnames<-`(c("marker",f))

  # Collect markers where there is some disagreement
  incomplete_founder_genos <- founder.allele.counts %>%
    dplyr::filter(n != max(founder.allele.counts$n)) %>%
    dplyr::select(-n) %>%
    `colnames<-`(c("marker",f))

  # Filter sample genotypes to markers with genotype disagreement
  incomplete_founder_genos_samples <- founder.geno.array %>%
    dplyr::filter(marker %in% unique(incomplete_founder_genos$marker)) %>%
    dplyr::select(sample_id, genotype, marker) %>% 
    dplyr::arrange(marker)

  # Join intensities to discordant genotyped samples
  incomplete_founder_genos_ints_samples <- founder.geno.array %>%
    dplyr::right_join(.,incomplete_founder_genos_samples) %>%
    dplyr::left_join(., gm_metadata %>% dplyr::select(marker, chr))

  # Create nested list by marker of sample genotypes and respective intensities to be able to eliminate certain genotypes off the bat based on outlier intensity values *within* a founder background
  incomplete_founder_genos_ints_samples_nested <- incomplete_founder_genos_ints_samples %>%
    dplyr::group_by(marker) %>%
    tidyr::nest()
  
  # Remove samples with extreme intensity values to try to create better consensus for the founder strain
  incomplete_founder_consensus <- purrr::pmap(.l = list(incomplete_founder_genos_ints_samples_nested$marker,
                                                        incomplete_founder_genos_ints_samples_nested$data,
                                                        rep(f, length(incomplete_founder_genos_ints_samples_nested$data))),
                                              .f = removeExtremeInts) %>%
    Reduce(rbind,.)

  # Identify markers where there is now 1 genotype across samples after removing extreme intensities
  recaptured_tally <- incomplete_founder_consensus %>%
    dplyr::group_by(marker) %>%
    dplyr::count() %>%
    dplyr::filter(n == 1)

  # Re-attach these re-inferred genotype calls to the consensus that already exists
  founder_genos <- complete_founder_genos %>%
    dplyr::bind_rows(.,incomplete_founder_consensus %>% 
                       dplyr::filter(marker %in% recaptured_tally$marker))

  # Assign these calls to a founder object
  assign(paste0("Calls_",f), founder_genos)

}

founder_consensus_calls <- list(`Calls_A/J`, `Calls_C57BL/6J`, `Calls_129S1/SvImJ`, `Calls_NOD/ShiLtJ`,
                               `Calls_NZO/HILtJ`, `Calls_CAST/EiJ`, `Calls_PWK/PhJ`, `Calls_WSB/EiJ`) %>% 
                               Reduce(full_join,.) %>% dplyr::ungroup()

# Identify consensus calls that have complete data across founders
founder_consensus_complete <- founder_consensus_calls[complete.cases(founder_consensus_calls),]

# Identify consensus calls that have DO NOT have complete data across founders (i.e. maybe some sort of discrepancy among samples contributing to that founder)
founder_consensus_incomplete <- founder_consensus_calls[!complete.cases(founder_consensus_calls),]

inbredFounderSamples <- founderSamples %>% dplyr::filter(bg == "INBRED")

missing_consensus_calls_nested <- founder_consensus_incomplete %>%
  tidyr::pivot_longer(-marker, names_to = "strain", values_to = "genotype") %>%
  dplyr::rename(consensus_genotype = genotype) %>%
  dplyr::left_join(., good_chr_genos %>%
                    dplyr::filter(marker %in% founder_consensus_incomplete$marker,
                                    sample_id %in% inbredFounderSamples$sample_id) %>% 
                    dplyr::select(marker, sample_id, genotype, X, Y) %>%
                    dplyr::rename(sample_genotype = genotype)) %>% 
                    dplyr::ungroup() %>%
                    dplyr::group_by(marker) %>%
                    tidyr::nest()

findConsensusGenotypes(mk = missing_consensus_calls_nested$marker[[1]], data =  missing_consensus_calls_nested$data[[1]])
