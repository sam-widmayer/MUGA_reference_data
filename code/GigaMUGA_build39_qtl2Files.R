# load required packages
library(data.table)
library(qtl2convert)
library(vroom)
library(dplyr)
library(tidyr)
writeMaps <- function(map, map_type){
  if(map_type == "genetic"){
    print(paste0("Writing genetic map for chromosome ",unique(map$chr)))
    qtl2convert::write2csv(df = map,
                           filename = paste0("output/GigaMUGA/GM_gmap",unique(map$chr),".csv"), 
                           comment = paste0("Genetic map for GigaMUGA chromosome ", unique(map$chr)),
                           overwrite = T)
  } else if(map_type == "physical"){
    print(paste0("Writing physical map for chromosome ",unique(map$chr)))
    qtl2convert::write2csv(df = map,
                           filename = paste0("output/GigaMUGA/GM_pmap",unique(map$chr),".csv"), 
                           comment = paste0("Physical map for GigaMUGA chromosome ", unique(map$chr)),
                           overwrite = T)
  } else {
    print("Invalid map type; skipping")
  }
}
writeFounderGenos <- function(founder_geno, founder_geno_chr){
  recodedFounderGeno <- founder_geno %>%
    tidyr::pivot_longer(-c(marker,`A`,`B`)) %>%
    dplyr::mutate(value = dplyr::case_when(value == `A` ~ "A",
                                           value == `B` ~ "B")) %>%
    tidyr::pivot_wider(names_from = name, values_from = value) %>%
    dplyr::select(-`A`, -`B`)
  
  recodedFounderGeno[,2:9][is.na(recodedFounderGeno[,2:9])] <- "-"
  
  colnames(recodedFounderGeno) <- c("marker",LETTERS[1:8])
  print(paste0("Writing founder genotypes for chromosome ", founder_geno_chr))
  qtl2convert::write2csv(df = recodedFounderGeno,
                         filename = paste0("output/GigaMUGA/GM_foundergeno",founder_geno_chr,".csv"), 
                         comment = paste0("Founder genotypes (A/B/-) for GigaMUGA chromosome ", unique(map$chr)),
                         overwrite = T)
  
}

# read in GM marker metadata
gm_ann <- vroom::vroom("data/GigaMUGA/gm_uwisc_v2.csv")

# read in QC'd founder consensus genotypes
founder_genos_path <- "output/GigaMUGA/GigaMUGA_founder_consensus_genotypes.csv.gz"
system(paste("gunzip", founder_genos_path))
founder_genos <- vroom::vroom("output/GigaMUGA/GigaMUGA_founder_consensus_genotypes.csv")

# Restrict annotation files to markers with HQ consensus genotypes
# and convert SNP genotypes to allele codes
allele_codes <- gm_ann[which(gm_ann$marker %in% founder_genos$marker),] %>%
  dplyr::select(marker, chr,  snp) %>%
  tidyr::separate(snp, into = c("A","B"), sep = -1)

# Write build 39 allele codes
write.csv(allele_codes, "output/GigaMUGA/GM_allele_codes_GRCm39.csv")
qtl2convert::write2csv(df = allele_codes, 
                       filename = "output/GigaMUGA/GM_allelecodes.csv", 
                       comment = "Allele codes for GigaMUGA",
                       overwrite = T)

# Filter to genetic map with HQ consensus genos
gmap <- gm_ann[which(gm_ann$marker %in% founder_genos$marker),] %>%
  dplyr::select(marker, chr,  cM_cox) %>%
  dplyr::rename(pos = cM_cox) %>%
  dplyr::mutate(chr_dummy = chr) %>%
  dplyr::group_by(chr_dummy) %>%
  tidyr::nest()

# Filter to physical map with HQ consensus genos
pmap <- gm_ann[which(gm_ann$marker %in% founder_genos$marker),] %>%
  dplyr::select(marker, chr,  bp_grcm39) %>%
  dplyr::rename(pos = bp_grcm39) %>%
  dplyr::mutate(chr_dummy = chr) %>%
  dplyr::group_by(chr_dummy) %>%
  tidyr::nest()

# Specify map types for each loop
map_type_genetic <- rep("genetic",length(gmap$data))
map_type_physical <- rep("physical",length(pmap$data))

# Write build 39 genetic maps
purrr::map2(.x = gmap$data, 
            .y = map_type_genetic,
            .f = writeMaps)

# Write build 39 phyiscal maps
purrr::map2(.x = pmap$data, 
            .y = map_type_physical,
            .f = writeMaps)

# Attach allele codes to CC/DO founder genotypes, 
# order by genetic position, and
# nest by chromosome
founderGenos_codeJoined <- founder_genos %>%
  dplyr::left_join(., allele_codes) %>%
  dplyr::right_join(gm_ann %>%
                      dplyr::arrange(chr, cM_cox) %>% 
                      dplyr::select(marker, chr, cM_cox),.) %>%
  dplyr::select(-cM_cox) %>%
  dplyr::group_by(chr) %>%
  tidyr::nest()

# Recode founder genotypes in (A/B/-) format and write in qtl2 style
purrr::map2(.x = founderGenos_codeJoined$data, 
            .y = founderGenos_codeJoined$chr, 
            .f = writeFounderGenos)

