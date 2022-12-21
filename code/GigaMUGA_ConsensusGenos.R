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
# Function to form F1 genotypes from founder consensus genotypes
# Inputs:
# x = Row of genotype values
# Outputs:
# Numeric vector of single letter genotypes
###########################################
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

###########################################
# Function to form F1 genotypes from founder consensus genotypes
# Inputs:
# x = Row of genotype values
# Outputs:
# Numeric vector of single letter genotypes
###########################################
callHemiGeno <- function(x){
  # Mitochondria comes from dam
  predicted.geno <- x[[2]]
  return(predicted.geno)
}

###########################################
# Function to form reference sample genotypes for comparison
# Inputs:
# Consensus genotypes for dam strain
# Consensus genotypes for sire strain

# Outputs:
# Numeric vector of recoded genotypes
###########################################
founder_background_QC <- function(dam, sire, GC_15 = F){
  
  # Extract strain names from genotype objects
  dam.df <- data.frame(dam)
  sire.df <- data.frame(sire)
  # dam_strain <- gsub(colnames(dam.df)[2], pattern = "X", replacement = "")
  # sire_strain <- gsub(colnames(sire.df)[2], pattern = "X", replacement = "")
  
  # Identify samples from supplied strains
  cross_bg <- paste0(colnames(dam)[2], colnames(sire)[2])
  
  if(GC_15 == F){
    maternalF1s <- founderSamples %>%
      dplyr::mutate(bg = stringr::str_sub(sample_id, start = 0, end = 2)) %>%
      dplyr::filter(bg == cross_bg)
    
  } else {
    maternalF1s <- founderSamples_GC_15_F %>%
      dplyr::mutate(bg = "EC") %>%
      dplyr::filter(bg == "EC")
  }
  
  print(paste0("Running QC: ", cross_bg))
  
  if(nrow(maternalF1s) == 0){
    # Some crosses don't exist in reference data, can't be QC'd
    return("No samples from this cross; skipping")
  } else {
    
    # Remove missing consensus genotypes for each strain
    mom <- dam.df[which(!dam.df[,2] == "-"),]
    dad <- sire.df[which(!sire.df[,2] == "-"),]
    
    # Identify informative markers common to both parental strains
    common_markers <- intersect(mom$marker, dad$marker)
    
    # Form the F1 hybrid
    cross <- dplyr::inner_join(mom, dad, "marker") %>%
      dplyr::filter(marker %in% common_markers)
    
    # Predict F1 genotypes from consensus genotypes for each strain
    predicted.genotypes <- apply(cross, 1, FUN = callGeno)
    cross$predicted_genotypes <- predicted.genotypes
  
    # Gather genotypes for F1 samples and remove bad markers
    sample_genos <- maternalF1s %>%
      dplyr::right_join(good_chr_genos,.,by = "sample_id") %>%
      dplyr::filter(marker %in% common_markers,
                    genotype != "N") %>%
      dplyr::distinct(marker, sample_id, genotype, .keep_all = TRUE)
    
    # Nest genotypes by sample and sex to prepare for calling Chr X genotypes
    nested_samples <- sample_genos %>%
      dplyr::group_by(sample_id, inferred.sex) %>%
      tidyr::nest()
    
    final_geno_comp_list <- list()
    for(i in 1:length(nested_samples$sample_id)){
      
      if(args[1] == "X" & nested_samples$inferred.sex[[i]] == "m"){
        print("Detected Male X Chromosome")
        # Male X Chrs
        new_predicted_genotypes <- apply(cross, 1, callHemiGeno)
        new_cross <- cross %>%
          dplyr::mutate(predicted_genotypes = new_predicted_genotypes)
          
          final_geno_comp_list[[i]] <- nested_samples$data[[i]] %>%
            dplyr::inner_join(., new_cross %>% 
                                dplyr::select(marker, predicted_genotypes), 
                              by = "marker") %>%
            dplyr::mutate(matching_genos = dplyr::if_else(genotype == predicted_genotypes,
                                                          true = "MATCH",
                                                          false = "NO MATCH"),
                          alt_chr = dplyr::case_when(chr == "M" ~ "M",
                                                     chr == "X" ~ "X",
                                                     chr == "Y" ~ "Y",
                                                     is.na(chr) ~ "Other",
                                                     TRUE ~ "Autosome"),
                          alt_chr = as.factor(alt_chr),
                          sample = nested_samples$sample_id[[i]],
                          inferred.sex = nested_samples$inferred.sex[[i]])
          } else {
        # Autosome
            print("Detected Female X Chromosome or Autosome")
            final_geno_comp_list[[i]] <- nested_samples$data[[i]] %>%
          dplyr::inner_join(., cross %>% 
                              dplyr::select(marker, predicted_genotypes), 
                            by = "marker") %>% 
          dplyr::mutate(matching_genos = dplyr::if_else(genotype == predicted_genotypes,
                                                        true = "MATCH",
                                                        false = "NO MATCH"),
                        alt_chr = dplyr::case_when(chr == "M" ~ "M",
                                                   chr == "X" ~ "X",
                                                   chr == "Y" ~ "Y",
                                                   is.na(chr) ~ "Other",
                                                   TRUE ~ "Autosome"),
                        alt_chr = as.factor(alt_chr),
                        sample = nested_samples$sample_id[[i]],
                        inferred.sex = nested_samples$inferred.sex[[i]])
        }
      }
    final_geno_comp <- Reduce(rbind,final_geno_comp_list)
    
    
    # Tabulate the percentage of concordant genotypes between theoretical and actual F1 samples
    geno_comp_summary <- final_geno_comp %>%
      dplyr::group_by(sample, inferred.sex, alt_chr, matching_genos) %>%
      dplyr::count() %>%
      tidyr::pivot_wider(names_from = matching_genos, values_from = n) %>%
      dplyr::mutate(bg = cross_bg,
                    chr = unique(good_chr_genos$chr))
    
    # Return genotypes and their summary statistics
    return(list(final_geno_comp,geno_comp_summary))
    }
}
# founder_background_QC(dam = dam_calls[[1]], 
#                       sire = sire_calls[[1]])


args <- commandArgs(trailingOnly = TRUE)
load("data/GigaMUGA/GigaMUGA_QC_Results.RData")
load("data/GigaMUGA/GigaMUGA_BadSamples_BadMarkers.RData")
load("data/GigaMUGA/GigaMUGA_SexCheck_Results.RData")
gm_metadata <- vroom::vroom("data/GigaMUGA/gm_uwisc_v2.csv",
                            progress = T)

# Sample Genotype Filepath
geno_path <- c("data/GigaMUGA/GigaMUGA_reference_genotypes")

# Consensus Genotype Filepath
consensus_geno_path <- c("data/GigaMUGA/GigaMUGA_founder_consensus_genotypes")

# Read in sample genotypes for a given chromosome and 
# filter out bad markers from sex check and marker QC.
# Using *ALL* reference samples, not just CC/DO founders and their F1s
good_chr_genos <- fst::read.fst(path = file.path(geno_path, 
                                                 paste0("gm_genos_chr_",args[1],".fst"))) %>%
  dplyr::mutate(marker_flag = dplyr::if_else(condition = marker %in% above.cutoff$marker,
                                             true = "FLAG",
                                             false = ""),
                genotype = dplyr::if_else(condition = genotype == "-",
                                          true = "N",
                                          false = genotype)) %>%
  dplyr::filter(marker_flag == "") %>%
  dplyr::select(-marker_flag)


# Read in consensus genotypes created from DO sample imputation
consensus_genotypes <- vroom::vroom(file = file.path(consensus_geno_path,
                                                     paste0("GigaMUGA_founder_consensus_imputed_genotypes_chr",
                                                            args[1],".csv")), delim = ",")
colnames(consensus_genotypes) <- c("marker",LETTERS[1:8])

founder_letters <- LETTERS[1:8]

# First form a list of dams that comprise each F1 cross type
parents <- data.frame(tidyr::expand_grid(founder_letters, founder_letters, .name_repair = "minimal")) %>%
  `colnames<-`(c("dams","sires"))

# Select dams
dams <- parents %>%
  dplyr::select(dams) %>%
  as.list()

# Select sires
sires <- parents %>%
  dplyr::select(sires) %>%
  as.list()

# Pull in consensus genotypes for each parental strain as a big list of objects
condenseFounderGenos <- function(x){
  # Identify founder in consensus file
  founder_consensus <- consensus_genotypes[,c("marker",x)]
  
  # Identify homozygous sites in consensus for the strain of interest
  fc_string <- strsplit(c(founder_consensus[,2])[[1]],"")
  homo_sites <- unlist(lapply(fc_string, function(x){x[1]==x[2]}))
  
  # Make founder strain genos the single letter genos
  single_letter_genos <- unlist(purrr::transpose(fc_string[which(homo_sites)])[[1]])
  founder_consensus[which(homo_sites),2] <- single_letter_genos
  
  return(founder_consensus)
}

dam_calls <- purrr::map(dams$dams,
                        condenseFounderGenos)
sire_calls <- purrr::map(sires$sires,
                         condenseFounderGenos)


# Compare the predicted genotypes (from consensus calls) to the actual genotypes of each sample in parallel
plan(multisession, workers = 16)
make_chunks <- furrr:::make_chunks
make_chunks(n_x = length(dam_calls), n_workers = 16)
bg_QC <- furrr::future_map2(dam_calls,
                            sire_calls, 
                            founder_background_QC,
                            .options = furrr_options(seed = TRUE))

# Keep outputs from the QC that are lists; if QC wasn't performed for a given background, the output was a character vector warning
founder_background_QC_tr <- bg_QC %>%
  purrr::keep(., is.list) %>%
  # Instead of having 64 elements of lists of two, have two lists of 64:
  # 1) All good genotypes from each cross with concordance values
  # 2) All concordance summaries for each cross
  purrr::transpose(.)

# Bind together all concordance summaries
founder_concordance_df <- Reduce(rbind, founder_background_QC_tr[[2]])

# If all markers for a given chromosome type were either concordant or discordant, NAs are returned
# This step assigns those NA values a 0
founder_concordance_df[is.na(founder_concordance_df)] <- 0

# Form concordance as a percentage
founder_concordance_df_2 <- founder_concordance_df %>%
  # dplyr::mutate(dam = gsub(dam, pattern = "[.]", replacement = "/"),
  #               sire = gsub(sire, pattern = "[.]", replacement = "/")) %>%
  dplyr::select(sample, inferred.sex, alt_chr, bg, chr, MATCH, `NO MATCH`)
founder_concordance_df_2$alt_chr <- factor(founder_concordance_df_2$alt_chr,
                                           levels = c("Autosome","X","Y","M","Other"))

# Performing background QC of GC_15 as the proper strain
founderSamples_GC_15_F <- founderSamples %>%
  dplyr::filter(sample_id == "GC_15_F")
GC_15_F_results <- founder_background_QC(dam = dam_calls[[33]],  # NZO
                                         sire = sire_calls[[3]], #129S1
                                         GC_15 = T) 


# Writing Concordance Files
if(dir.exists(paths = "data/GigaMUGA/GigaMUGA_founder_sample_concordance/")){
  fst::write.fst(founder_concordance_df_2, 
                 path = paste0("data/GigaMUGA/GigaMUGA_founder_sample_concordance/gm_founder_concordance_",args[1],".fst"))
} else {
  dir.create(path = "data/GigaMUGA/GigaMUGA_founder_sample_concordance/")
  fst::write.fst(founder_concordance_df_2, 
                 path = paste0("data/GigaMUGA/GigaMUGA_founder_sample_concordance/gm_founder_concordance_",args[1],".fst"))
}

# Write Results for Resexing of GC_15 to an NZO/129S1 mouse
if(dir.exists(paths = "data/GigaMUGA/GC_15_F_ReSex_Results")){
  write.csv(GC_15_F_results[[2]], file = paste0("data/GigaMUGA/GC_15_F_ReSex_Results/GC_15_F_",args[1],".csv"), row.names = F)
} else {
  dir.create(path = "data/GigaMUGA/GC_15_F_ReSex_Results/")
  write.csv(GC_15_F_results[[2]], file = paste0("data/GigaMUGA/GC_15_F_ReSex_Results/GC_15_F_",args[1],".csv"), row.names = F)
}












# From previous iteration in which consensus calls are generated from intensities
# 
# List of founder strains
# founder_strains <- c("A/J","C57BL/6J","129S1/SvImJ","NOD/ShiLtJ",
#                      "NZO/HILtJ","CAST/EiJ","PWK/PhJ","WSB/EiJ")
# 
# # 
# 
# # for(f in founder_strains){
# #   print(paste("Generating Calls for",f))
# #   
# #   # Pulling reference sample metadata for an individual CC/DO founder strain
# #   founder.geno.array <- founderSamples %>%
# #     dplyr::filter(dam == f,
# #                   sire == f) %>%
# #     # Attaching genotypes for good markers for each sample from that strain
# #     dplyr::left_join(.,good_chr_genos, by = "sample_id")
# #   
# #   
# #   # Count the number of unique allele calls for each marker
# #   founder.allele.counts <- founder.geno.array %>%
# #     dplyr::group_by(marker, genotype) %>%
# #     dplyr::count()
# #   
# #   # Collect markers which have identical genotypes across samples derived from the same founder strain
# #   complete_founder_genos <- founder.allele.counts %>%
# #     dplyr::filter(n == max(founder.allele.counts$n)) %>%
# #     dplyr::select(-n) %>%
# #     `colnames<-`(c("marker",f))
# #   
# #   # Collect markers where there is some disagreement
# #   incomplete_founder_genos <- founder.allele.counts %>%
# #     dplyr::filter(n != max(founder.allele.counts$n)) %>%
# #     dplyr::select(-n) %>%
# #     `colnames<-`(c("marker",f))
# #   
# #   # Filter sample genotype array to markers with genotype disagreement
# #   incomplete_founder_genos_samples <- founder.geno.array %>%
# #     dplyr::filter(marker %in% unique(incomplete_founder_genos$marker)) %>%
# #     dplyr::select(sample_id, genotype, marker) %>% 
# #     dplyr::arrange(marker)
# #   
# #   # Join intensities to discordant genotyped samples
# #   incomplete_founder_genos_ints_samples <- founder.geno.array %>%
# #     dplyr::right_join(.,incomplete_founder_genos_samples) %>%
# #     dplyr::left_join(., gm_metadata %>% dplyr::select(marker, chr))
# #   
# #   # Create nested list by marker of sample genotypes and respective intensities to be able to eliminate certain genotypes off the bat based on outlier intensity values *within* a founder background
# #   incomplete_founder_genos_ints_samples_nested <- incomplete_founder_genos_ints_samples %>%
# #     dplyr::group_by(marker) %>%
# #     tidyr::nest()
# #   
# #   # Remove samples with extreme intensity values to try to create better consensus for the founder strain
# #   incomplete_founder_consensus <- purrr::pmap(.l = list(incomplete_founder_genos_ints_samples_nested$marker,
# #                                                         incomplete_founder_genos_ints_samples_nested$data,
# #                                                         rep(f, length(incomplete_founder_genos_ints_samples_nested$data))),
# #                                               .f = removeExtremeInts) %>%
# #     Reduce(rbind,.)
# #   
# #   # Identify markers where there is now 1 genotype across samples after removing extreme intensities
# #   recaptured_tally <- incomplete_founder_consensus %>%
# #     dplyr::group_by(marker) %>%
# #     dplyr::count() %>%
# #     dplyr::filter(n == 1)
# #   
# #   # Re-attach these re-inferred genotype calls to the consensus that already exists
# #   founder_genos <- complete_founder_genos %>%
# #     dplyr::bind_rows(.,incomplete_founder_consensus %>% 
# #                        dplyr::filter(marker %in% recaptured_tally$marker))
# #   
# #   # Assign these calls to a founder object
# #   assign(paste0("Calls_",f), founder_genos)
# #   
# # }
# # 
# # founder_consensus_calls <- list(`Calls_A/J`, `Calls_C57BL/6J`, `Calls_129S1/SvImJ`, `Calls_NOD/ShiLtJ`,
# #                                 `Calls_NZO/HILtJ`, `Calls_CAST/EiJ`, `Calls_PWK/PhJ`, `Calls_WSB/EiJ`) %>% 
# #   Reduce(full_join,.) %>% dplyr::ungroup()
# # 
# # # Identify consensus calls that have complete data across founders
# # founder_consensus_complete <- founder_consensus_calls[complete.cases(founder_consensus_calls),]
# # 
# # # Identify consensus calls that have DO NOT have complete data across founders (i.e. maybe some sort of discrepancy among samples contributing to that founder)
# # founder_consensus_incomplete <- founder_consensus_calls[!complete.cases(founder_consensus_calls),]
# # 
# # inbredFounderSamples <- founderSamples %>% 
# #   dplyr::filter(bg == "INBRED") %>%
# #   dplyr::mutate(strain = dam)
# # 
# # missing_consensus_calls_nested <- founder_consensus_incomplete %>%
# #   tidyr::pivot_longer(-marker, names_to = "strain", values_to = "genotype") %>%
# #   dplyr::rename(consensus_genotype = genotype) %>%
# #   dplyr::left_join(., good_chr_genos %>%
# #                      dplyr::select(marker, sample_id, genotype, X, Y) %>%
# #                      dplyr::filter(marker %in% founder_consensus_incomplete$marker,
# #                                    sample_id %in% inbredFounderSamples$sample_id) %>% 
# #                      dplyr::rename(sample_genotype = genotype) %>%
# #                      dplyr::left_join(., inbredFounderSamples)) %>% 
# #   dplyr::ungroup() %>%
# #   dplyr::group_by(marker) %>%
# #   tidyr::nest()
# # 
# # plan(multisession, workers = 16)
# # make_chunks <- furrr:::make_chunks
# # make_chunks(n_x = length(missing_consensus_calls_nested$marker), n_workers = 16)
# # founder_consensus_incomplete_recoded <- suppressWarnings(furrr::future_map2(missing_consensus_calls_nested$marker,
# #                                                                             missing_consensus_calls_nested$data,
# #                                                                             findConsensusGenotypes, 
# #                                                                             .options = furrr_options(seed = TRUE)))
# # founder_consensus_incomplete_recoded_tr <- purrr::transpose(founder_consensus_incomplete_recoded)
# # new_consensus_founders <- Reduce(dplyr::bind_rows,founder_consensus_incomplete_recoded_tr[[1]])
# # updated_founder_sample_genotypes <- Reduce(dplyr::bind_rows,founder_consensus_incomplete_recoded_tr[[2]])
# # 
# # # Update original consensus calls with re-assigned consensus calls
# # clean_founder_consensus_genotypes <- founder_consensus_complete %>%
# #   dplyr::bind_rows(.,new_consensus_founders)
# # 
# # original_founder_sample_genotypes <- good_chr_genos %>%
# #   dplyr::filter(sample_id %in% colnames(updated_founder_sample_genotypes)[-1],
# #                 !marker %in% updated_founder_sample_genotypes$marker) %>%
# #   dplyr::select(sample_id, marker, genotype) %>%
# #   tidyr::pivot_wider(names_from = sample_id, values_from = genotype)
# # 
# # founder_sample_genotypes <- dplyr::bind_rows(updated_founder_sample_genotypes,original_founder_sample_genotypes)
# # 
# # 
# # # founder_strains_2 <- gsub(founder_strains, pattern = "[.]", replacement =  "/")
# # for(f in founder_strains){
# #   clean_founder_con <- clean_founder_consensus_genotypes %>%
# #     dplyr::ungroup() %>%
# #     dplyr::select(marker, f)
# #   assign(paste0("Clean_Calls_",f), clean_founder_con)
# # }
# # 


# for(i in 1:length(dams$dams)){
#   
#   # Identify founder in consensus file
#   founder_consensus <- consensus_genotypes[,c("marker",dams$dams[i])]
#   
#   # Identify homozygous sites in consensus for the strain of interest
#   fc_string <- strsplit(c(founder_consensus[,2])[[1]],"")
#   homo_sites <- unlist(lapply(fc_string, function(x){x[1]==x[2]}))
#   
#   # Make founder strain genos the single letter genos
#   single_letter_genos <- unlist(purrr::transpose(fc_string[which(homo_sites)])[[1]])
#   founder_consensus[which(homo_sites),2] <- single_letter_genos
#   
#   dam_calls[[i]] <- founder_consensus
# }
# 
# 
# sire_calls <- list()
# for(i in 1:length(sires$sires)){
#   sire_calls[[i]] <- get(ls(pattern = paste0("Clean_Calls_",sires$sires[i])))
# }


# ###########################################
# # Function to re-derive founder consensus genotypes for founders that are missing consensus calls due to one or a few bad calls
# # Inputs:
# # mk = marker name
# # data = data frame with founder strain samples, intensities, individual sample genotype calls , and existing consensus genotype calls
# 
# # Outputs:
# # 1) Consensus genotype calls for each founder strain
# # 2) Wide genotype table for all founder samples
# # 3) Data frame with recoded consensus genotypes for each samples with intensity data and strain names, as well as a flag for whether the sample contributed to a recoded consensus genotype
# ###########################################
# findConsensusGenotypes <- function(mk, data){
#   
#   # Examples use mk = "UNC13666424"
#   
#   # Annotate sample genotype data with summed intensities as a heuristic and filter to sample genotypes assigned a real genotype
#   call_filtered_data <- data %>%
#     dplyr::mutate(sum_int = X+Y) %>%
#     dplyr::filter(sample_genotype %in% c("A","C","T","G","H"))
#   
#   # Derive a lower bound by which to identify true N calls using the intensity data
#   real_N_cutoff <- quantile(call_filtered_data$sum_int, probs = seq(0,1,0.05))[[2]]
#   # Use this threshold to identify potential miscalls and remove them in order to recode samples without consensus calls
#   filtered_data <- call_filtered_data %>%
#     dplyr::mutate(test = if_else(sum_int < real_N_cutoff & is.na(consensus_genotype), true = "miscall", false = "")) %>%
#     dplyr::filter(test != "miscall")
#   
#   # Determine the number of clusters to use for k-means clustering of intensity values;
#   # If the previous steps have succeeded, there should be one or two real genotypes segregating and samples can be re-called by k-means
#   n_genos <- length(unique(filtered_data$sample_genotype[which(filtered_data$sample_genotype %in% c("A","C","T","G","H"))]))
#   ints_for_kmeans <- filtered_data %>%
#     dplyr::select(X, Y)
#   consensus_geno_clusters <- kmeans(ints_for_kmeans, centers = n_genos)$cluster
#   
#   # Joining each sample's cluster assignment to the sample intensity metrics
#   genos.k.means <- filtered_data %>%
#     dplyr::mutate(clust = as.factor(consensus_geno_clusters))
#   
#   # Create a join table to pair up clusters with genotypes
#   gens_by_cluster_tab <- genos.k.means %>%
#     dplyr::group_by(sample_genotype, consensus_genotype, clust) %>%
#     dplyr::count() %>%
#     arrange(-n)
#   # Example:
#   #   sample_genotype consensus_genotype clust     n
#   # 1 G               G                  2        32
#   # 2 A               A                  1        23
#   # 3 A               NA                 1         6
#   
#   # Filter the join table to hopefully have a 1:1 relationship between consensus genotypes and clusters
#   top_clusters <- gens_by_cluster_tab %>%
#     dplyr::filter(!is.na(consensus_genotype)) %>%
#     dplyr::ungroup() %>%
#     dplyr::distinct(sample_genotype, clust) %>%
#     dplyr::rename(consensus_genotype = sample_genotype)
#   # Example:
#   #   consensus_genotype clust
#   # 1 G                  2    
#   # 2 A                  1    
#   
#   # Samples are then recoded according to the k-means assigned genotypes and noted whether a consensus was re-assigned
#   unknown_before <- filtered_data %>%
#     dplyr::filter(is.na(consensus_genotype))
#   recoded_consensus_genotypes <- genos.k.means %>%
#     dplyr::select(-consensus_genotype) %>%
#     dplyr::left_join(.,top_clusters, by = "clust") %>%
#     dplyr::mutate(marker = mk) %>% 
#     dplyr::mutate(recoded = dplyr::if_else(sample_id %in% unknown_before$sample_id, true = "RECODED", false = "")) %>%
#     dplyr::select(-clust)
#   # Example: recoded_consensus_genotypes[20:30,]
#   # Note NOD samples are now recoded with the consensus genotype predicted from k-means across all founder samples
#   #      strain      sample              x_int y_int sample_genotype sum_int test  consensus_genotype marker      recoded  
#   #    <chr>       <chr>               <dbl> <dbl> <chr>             <dbl> <chr> <chr>              <chr>       <chr>    
#   #  1 129S1/SvImJ X129S1.SvImJf0827.1 0.716 0.083 A                 0.799 ""    A                  UNC13666424 ""       
#   #  2 129S1/SvImJ X129S1.SvImJm       0.618 0.079 A                 0.697 ""    A                  UNC13666424 ""       
#   #  3 129S1/SvImJ X129S1.SvImJm.1     0.599 0.084 A                 0.683 ""    A                  UNC13666424 ""       
#   #  4 129S1/SvImJ X129S1.SvImJm1314   0.69  0.097 A                 0.787 ""    A                  UNC13666424 ""       
#   #  5 NOD/ShiLtJ  NOD.ShiLtJf0713.1   0.567 0.07  A                 0.637 ""    A                  UNC13666424 "RECODED"
#   #  6 NOD/ShiLtJ  NOD.ShiLtJf0713.2   0.6   0.095 A                 0.695 ""    A                  UNC13666424 "RECODED"
#   #  7 NOD/ShiLtJ  NOD.ShiLtJf0713.3   0.883 0.091 A                 0.974 ""    A                  UNC13666424 "RECODED"
#   #  8 NOD/ShiLtJ  NOD.ShiLtJm0150     0.646 0.088 A                 0.734 ""    A                  UNC13666424 "RECODED"
#   #  9 NOD/ShiLtJ  NOD.ShiLtJm35324    0.538 0.077 A                 0.615 ""    A                  UNC13666424 "RECODED"
#   # 10 NOD/ShiLtJ  NOD.ShiLtJm39173    0.649 0.087 A                 0.736 ""    A                  UNC13666424 "RECODED"
#   # 11 NZO/HILtJ   NZO.HILtJf0588      0.041 0.728 G                 0.769 ""    G                  UNC13666424 ""       
#   
#   
#   # Generate a wide table of consensus as expected for output file
#   recoded <- recoded_consensus_genotypes %>%
#     dplyr::distinct(marker, strain, consensus_genotype) %>% 
#     tidyr::pivot_wider(names_from = strain, values_from = consensus_genotype)
#   #     marker      `A/J` `C57BL/6J` `129S1/SvImJ` `NOD/ShiLtJ` `NZO/HILtJ` `CAST/EiJ` `PWK/PhJ` `WSB/EiJ`
#   #   <chr>       <chr> <chr>      <chr>         <chr>        <chr>       <chr>      <chr>     <chr>    
#   # 1 UNC13666424 G     G          A             A            G           A          A         G     
#   
#   # Check for any founders that still segregate a genotype within the marker
#   still_consensus <- recoded_consensus_genotypes %>%
#     dplyr::group_by(strain, consensus_genotype) %>%
#     dplyr::count() %>%
#     dplyr::group_by(strain) %>%
#     dplyr::count()
#   
#   # If there are still discrepancies (i.e. non-N genotypes segregating among samples), return the equivalent of the input data for downstream work
#   if(2 %in% still_consensus$n){
#     original_consensus <- suppressWarnings(data %>%
#                                              dplyr::mutate(marker = mk) %>%
#                                              dplyr::distinct(marker, strain, consensus_genotype) %>%
#                                              tidyr::pivot_wider(names_from = strain, values_from = consensus_genotype))
#     
#     original_sample_calls <- suppressWarnings(data %>%
#                                                 dplyr::mutate(marker = mk) %>%
#                                                 dplyr::distinct(marker, sample_id, sample_genotype) %>%
#                                                 tidyr::pivot_wider(names_from = sample_id,
#                                                                    values_from = sample_genotype))
#     original_data <- data %>%
#       dplyr::mutate(marker = mk,
#                     recoded = "")
#     
#     return(list(original_consensus, original_sample_calls, original_data))
#   } else{
#     sample_new_consensus_calls <- recoded_consensus_genotypes %>%
#       dplyr::distinct(marker, sample_id, consensus_genotype) %>% 
#       tidyr::pivot_wider(names_from = sample_id, 
#                          values_from = consensus_genotype)
#     return(list(recoded, sample_new_consensus_calls, recoded_consensus_genotypes))
#   }
# }



# # Form a list of mitochondrial markers
# alt_chr_M <- good_chr_genos %>%
#   dplyr::filter(chr %in% c("M")) %>%
#   dplyr::distinct(marker)
# # Filter the artificial F1 hybrid genotype calls down to mitochondrial markers
# mito_cross <- alt_chr_M %>% 
#   dplyr::left_join(., cross)
# # Assign maternal mitochondrial genotype to predicted cross instead of hets where strains have different mitotypes
# mito_cross$predicted_genotypes <- apply(mito_cross, 1, FUN = callHemiGeno)
# cross_mitorecoded <- cross %>%
#   dplyr::filter(!marker %in% alt_chr_M$marker) %>%
#   dplyr::bind_rows(.,mito_cross)

# Create a list of mitochondrial genotypes to be supplied to the callXGeno function
# cross_mitorecoded_list <- list()
# for(i in 1:length(nested_samples$data)){
#   cross_mitorecoded_list[[i]] <- cross_mitorecoded
# }

# # Form a list of Chr X markers
# alt_chr_X <- good_chr_genos %>%
#   dplyr::filter(chr %in% c("X")) %>%
#   dplyr::distinct(marker)

# Create a list of X genotypes to be supplied to the callXGeno function
# alt_chr_X_list <- list()
# for(i in 1:length(nested_samples$data)){
#   alt_chr_X_list[[i]] <- alt_chr_X
# }
# Recode X genotypes to match expectations for male and female samples and form a dataframe of sample genotypes (and concordance status)
# final_geno_comp <- purrr::pmap(.l = list(nested_samples$sample_id,
#                                          nested_samples$inferred.sex,
#                                          nested_samples$data,
#                                          cross_mitorecoded_list,
#                                          alt_chr_X_list), 
#                                .f = callXGeno) %>% 
# Reduce(rbind,.)


# 
# founder_background_QCGC_15_F <- function(dam, sire){
#   
#   # Extract strain names from genotype objects\
#   dam.df <- data.frame(dam)
#   sire.df <- data.frame(sire)
#   dam_strain <- gsub(colnames(dam.df)[2], pattern = "X", replacement = "")
#   sire_strain <- gsub(colnames(sire.df)[2], pattern = "X", replacement = "")
#   
#   if(dam_strain == sire_strain){
#     print(paste0("Running QC: ", dam_strain))  
#   } else {
#     print(paste0("Running QC: (", dam_strain, "x", sire_strain, ")F1"))  
#   }
#   
#   # Identify samples from supplied strains
#   maternalF1s <- founderSamples_GC_15_F %>%
#     dplyr::mutate(dam = gsub(dam, pattern = "/", replacement = "."),
#                   sire = gsub(sire, pattern = "/", replacement = ".")) %>%
#     dplyr::filter(dam == dam_strain,
#                   sire == sire_strain) %>%
#     dplyr::mutate(sire = as.factor(sire),
#                   dam = as.factor(dam))
#   
#   if(nrow(maternalF1s) == 0){
#     # Some crosses don't exist in reference data, can't be QC'd
#     return("No samples from this cross; skipping")
#   } else {
#     
#     # Remove consensus genotypes for each strains that are N's
#     mom <- dam.df[which(!dam.df[,2] %in% c("H","N")),]
#     dad <- sire.df[which(!sire.df[,2] %in% c("H","N")),]
#     
#     # Form a hypothetical F1 hybrid by combining genotypes from markers that exist in both strains
#     if(dam_strain == sire_strain){
#       cross <- dplyr::inner_join(mom[complete.cases(mom),],dad[complete.cases(dad),], "marker")
#     } else {
#       cross <- dplyr::inner_join(mom[complete.cases(mom),],dad[complete.cases(dad),])
#     }
#     
#     # Predict F1 genotypes from consensus genotypes for each strain
#     predicted.genotypes <- apply(cross, 1, FUN = callGeno)
#     cross$predicted_genotypes <- predicted.genotypes
#     
#     # Gather genotypes for F1 samples and remove bad markers
#     sample_genos <- maternalF1s %>%
#       dplyr::right_join(good_chr_genos,.,by = "sample_id") %>%
#       dplyr::filter(marker %in% mom$marker, 
#                     genotype != "N") %>%
#       dplyr::distinct(marker, sample_id, genotype, .keep_all = TRUE)
#     
#     # Nest genotypes by sample and sex to prepare for calling Chr X genotypes
#     nested_samples <- sample_genos %>%
#       dplyr::group_by(sample_id, inferred.sex) %>%
#       tidyr::nest()
#     
#     final_geno_comp_list <- list()
#     for(i in 1:length(nested_samples$sample_id)){
#       
#       if(args[1] == "M"){
#         # Mt
#         print("Detected Mt")
#         new_predicted_genotypes <- apply(cross, 1, callHemiGeno)
#         new_cross <- cross %>%
#           dplyr::mutate(predicted_genotypes = new_predicted_genotypes)
#         
#         final_geno_comp_list[[i]] <- nested_samples$data[[i]] %>%
#           dplyr::inner_join(., new_cross %>% 
#                               dplyr::select(marker, predicted_genotypes), 
#                             by = "marker") %>%
#           dplyr::mutate(matching_genos = dplyr::if_else(genotype == predicted_genotypes,
#                                                         true = "MATCH",
#                                                         false = "NO MATCH"),
#                         alt_chr = dplyr::case_when(chr == "M" ~ "M",
#                                                    chr == "X" ~ "X",
#                                                    chr == "Y" ~ "Y",
#                                                    is.na(chr) ~ "Other",
#                                                    TRUE ~ "Autosome"),
#                         alt_chr = as.factor(alt_chr),
#                         sample = nested_samples$sample_id[[i]],
#                         inferred.sex = nested_samples$inferred.sex[[i]])
#       } else if(args[1] == "X" & nested_samples$inferred.sex[[i]] == "m"){
#         print("Detected Male X Chromosome")
#         # Male X Chrs
#         new_predicted_genotypes <- apply(cross, 1, callHemiGeno)
#         new_cross <- cross %>%
#           dplyr::mutate(predicted_genotypes = new_predicted_genotypes)
#         
#         final_geno_comp_list[[i]] <- nested_samples$data[[i]] %>%
#           dplyr::inner_join(., new_cross %>% 
#                               dplyr::select(marker, predicted_genotypes), 
#                             by = "marker") %>%
#           dplyr::mutate(matching_genos = dplyr::if_else(genotype == predicted_genotypes,
#                                                         true = "MATCH",
#                                                         false = "NO MATCH"),
#                         alt_chr = dplyr::case_when(chr == "M" ~ "M",
#                                                    chr == "X" ~ "X",
#                                                    chr == "Y" ~ "Y",
#                                                    is.na(chr) ~ "Other",
#                                                    TRUE ~ "Autosome"),
#                         alt_chr = as.factor(alt_chr),
#                         sample = nested_samples$sample_id[[i]],
#                         inferred.sex = nested_samples$inferred.sex[[i]])
#       } else {
#         # Autosome
#         print("Detected Female X Chromosome or Autosome")
#         final_geno_comp_list[[i]] <- nested_samples$data[[i]] %>%
#           dplyr::inner_join(., cross %>% 
#                               dplyr::select(marker, predicted_genotypes), 
#                             by = "marker") %>% 
#           dplyr::mutate(matching_genos = dplyr::if_else(genotype == predicted_genotypes,
#                                                         true = "MATCH",
#                                                         false = "NO MATCH"),
#                         alt_chr = dplyr::case_when(chr == "M" ~ "M",
#                                                    chr == "X" ~ "X",
#                                                    chr == "Y" ~ "Y",
#                                                    is.na(chr) ~ "Other",
#                                                    TRUE ~ "Autosome"),
#                         alt_chr = as.factor(alt_chr),
#                         sample = nested_samples$sample_id[[i]],
#                         inferred.sex = nested_samples$inferred.sex[[i]])
#       }
#     }
#     final_geno_comp <- Reduce(rbind,final_geno_comp_list)
#     
#     
#     # Tabulate the percentage of concordant genotypes between theoretical and actual F1 samples
#     geno_comp_summary <- final_geno_comp %>%
#       dplyr::group_by(sample, inferred.sex, alt_chr, matching_genos) %>%
#       dplyr::count() %>%
#       tidyr::pivot_wider(names_from = matching_genos, values_from = n) %>%
#       dplyr::mutate(dam = dam_strain, 
#                     sire = sire_strain,
#                     chr = unique(good_chr_genos$chr))
#     
#     # Return genotypes and their summary statistics
#     return(list(final_geno_comp,geno_comp_summary))
#   }
# }



###########################################
# Function used in the loop which forms founder consensus genotypes prior to background QC
# Inputs:
# mk = marker
# data = genotype data for each founder sample
# f = founder strain
# Outputs:
# data frame with 1 row with columns: marker; consensus genotype for founder strain f 
# ###########################################
# removeExtremeInts <- function(mk, data, f){
#   
#   # Calculate summary statistics for probe intensities
#   sd.x.int <- sd(data$X)
#   sd.y.int <- sd(data$Y)
#   mean.x.int <- mean(data$X)
#   mean.y.int <- mean(data$Y)
#   
#   data %>%
#     # With a prior expectation provided by an "N" call, determine whether probe intensities are unusual and 
#     # flag samples that meet both criteria
#     dplyr::mutate(x.ex = dplyr::if_else(X > (mean.x.int + sd.x.int) | X < (mean.x.int - sd.x.int), true = "EX", false = ""),
#                   y.ex = dplyr::if_else(Y > (mean.y.int + sd.y.int) | Y < (mean.y.int - sd.y.int), true = "EX", false = ""),
#                   flag = dplyr::if_else(genotype == "N" & (x.ex == "EX" | y.ex == "EX"), true = "FLAG", false = "")) %>%
#     dplyr::filter(flag != "FLAG") %>%
#     dplyr::mutate(marker = mk) %>%
#     dplyr::distinct(marker,genotype) %>%
#     dplyr::select(marker,genotype) %>%
#     `colnames<-`(c("marker",f))
# }