#!/usr/bin/env Rscript
library(dplyr)
library(magrittr)
library(stringr)
library(vroom)
library(tidyr)
library(qtl2)
library(parallel)
library(fst)
library(furrr)


# Read in GM reference metadata
GM_founder_metadata <- vroom::vroom(file = "output/GigaMUGA/GigaMUGA_founder_metadata.csv")
# GM_reference_genos <- vroom::vroom(file = "output/GigaMUGA/GigaMUGA_founder_consensus_genotypes.csv.gz")
# colnames(GM_reference_genos) <- c("marker", LETTERS[1:8])
# GM_reference_sample_genos <- vroom::vroom(file = "output/GigaMUGA/GigaMUGA_genotypes.csv.gz")
# GM_reference_sample_genos %<>%
#   dplyr::select(marker, GM_founder_metadata$sample)
GM_marker_metadata <- vroom::vroom(file = "data/GigaMUGA/gm_uwisc_v2.csv")
GM_original_reference_geno_path <- "data/GigaMUGA/GigaMUGA_reference_genotypes/"
GM_original_reference_genos <- list.files(GM_original_reference_geno_path)


# Read in DO dataset
DO_lifespan_genos <- vroom::vroom(file = "data/GigaMUGA/DO_churchill_lifespan_1/Churchill-churchill2-GigaMUGA_geno.csv", 
                                  num_threads = parallel::detectCores())
DO_lifespan_covar <- vroom::vroom(file = "data/GigaMUGA/DO_churchill_lifespan_1/Churchill-churchill2-GigaMUGA_covar.csv")
DO_lifespan_foundergenos <- vroom::vroom(file = "data/GigaMUGA/DO_churchill_lifespan_1/Churchill-churchill2-GigaMUGA_foundergeno.csv")
DO_lifespan_genoprobs <- readRDS(file = "data/GigaMUGA/DO_churchill_lifespan_1/churchill2__GigaMUGA_genoprobs_36state.rds")


# Count no calls among DO founders in the GigaMUGA reference data
GM_reference_no_calls_list <- list()
for(g in GM_original_reference_genos){
  cat(paste0("Reading in ",g))
  raw_geno <- fst::read.fst(path = paste0(GM_original_reference_geno_path,g))
  raw_geno %<>% 
    dplyr::select(marker, genotype, sample_id) %>%
    dplyr::filter(sample_id %in% GM_founder_metadata$sample)
    
  n_samples <- length(unique(raw_geno$sample_id))
  
  wide_raw_geno <- raw_geno %>%
    tidyr::pivot_wider(names_from = sample_id, values_from = genotype)
  
  no_calls_raw_genos <- apply(wide_raw_geno[,-1], 
                              1, 
                              function(x) sum(stringr::str_count(as.character(x),pattern = "-"))/n_samples)
  GM_reference_no_calls <- data.frame(wide_raw_geno$marker, no_calls_raw_genos)
  colnames(GM_reference_no_calls) <- c("marker","GM_no_calls")
  GM_reference_no_calls <- dplyr::left_join(GM_reference_no_calls, GM_marker_metadata)
  GM_reference_no_calls_list[[which(g == GM_original_reference_genos)]] <- GM_reference_no_calls
  }
GM_reference_no_calls <- Reduce(dplyr::bind_rows, GM_reference_no_calls_list)
redundant_markers <- GM_marker_metadata[which(!GM_marker_metadata$marker %in% GM_reference_no_calls$marker),]


# Count no calls among DO samples in lifespan data
DO_lifespan_no_calls <- apply(DO_lifespan_genos[,-1], 
                              1, 
                              function(x) (ncol(DO_lifespan_genos[,-1])-table(is.na(x))[[1]])/ncol(DO_lifespan_genos[,-1]))
DO_lifespan_no_calls <- data.frame(DO_lifespan_genos$marker, DO_lifespan_no_calls)
colnames(DO_lifespan_no_calls) <- c("marker","DO_no_call_rate")


# Pair up the two no call datasets
no_call_comp <- dplyr::left_join(GM_reference_no_calls, DO_lifespan_no_calls) %>%
  dplyr::select(marker, chr, bp_grcm39, snp, GM_no_calls, DO_no_call_rate)


# Plotting no calls in DO samples against DO founder reference samples
ggplot(data = no_call_comp[which(no_call_comp$DO_no_call_rate != 0),], mapping = aes(x = GM_no_calls, 
                                          y = DO_no_call_rate)) + 
  theme_bw() + 
  geom_point(alpha = 0.2) +
  # geom_smooth(method = "lm") + 
  theme(panel.grid = element_blank()) + 
  labs(x = "No call rate in GigaMUGA reference samples (%)", 
       y = "No call rate in DO lifespan study (%)")

# Finding the inferred diplotypes at each marker for all DO samples
GM_reference_no_calls$chr <- as.factor(GM_reference_no_calls$chr)
GM_reference_no_calls$chr <- factor(GM_reference_no_calls$chr,
                                    levels = c(paste(seq(1:19)),"X","Y","M"))
GM_reference_no_calls_nested <- GM_reference_no_calls %>%
  dplyr::arrange(chr) %>%
  dplyr::group_by(chr) %>%
  tidyr::nest()

# findMarkerDiplotypes(GM_reference = GM_reference_no_calls_nested$data[[1]], 
#                      DO_genoprobs = DO_lifespan_genoprobs[[1]])


findMarkerDiplotypes <- function(GM_reference, DO_genoprobs){

  max_marker_diplos <- list()
  for(i in 1:length(GM_reference$marker)){
    if(!GM_reference$marker[i] %in% names(DO_genoprobs[1,1,])){
      max_marker_diplos[[i]] <- paste0(GM_reference$marker[i]," not in genoprobs")
    } else {
      DO_marker_diplotypes <- data.frame(DO_genoprobs[,,GM_reference$marker[i]])
      rownames(DO_marker_diplotypes) <- unlist(purrr::transpose(strsplit(rownames(DO_marker_diplotypes),"_"))[[7]])
      DO_marker_diplotypes$sample <- as.factor(rownames(DO_marker_diplotypes))
      rownames(DO_marker_diplotypes) <- NULL
      DO_samples_top_diplotype <- DO_marker_diplotypes %>%
        tidyr::pivot_longer(-sample, names_to = "diplotype", values_to = "prob") %>%
        dplyr::group_by(sample) %>%
        dplyr::filter(prob == max(prob, na.rm=TRUE)) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(sample) %>%
        dplyr::mutate(marker = GM_reference$marker[i],
                      GM_no_call_rate = GM_reference$GM_no_calls[i],
                      chr = unique(GM_reference$chr))
      max_marker_diplos[[i]] <- DO_samples_top_diplotype
    }
  }
  max_marker_diplos_df <- purrr::keep(max_marker_diplos, is.data.frame) %>%
    Reduce(dplyr::bind_rows,.)
  return(max_marker_diplos_df)
}
future::plan(multisession, workers = 16)
make_chunks <- furrr:::make_chunks
make_chunks(n_x = length(GM_reference_no_calls_list[1:20]), n_workers = 16)
max_marker_diplos_df <- furrr::future_map2(.x = GM_reference_no_calls_list[1:20], 
                                           .y = DO_lifespan_genoprobs,
                                           .f = findMarkerDiplotypes, 
                                           .options = furrr_options(seed = TRUE))

save(max_marker_diplos_df, file = "data/GigaMUGA/DO_marker_diplotypes.RData")






