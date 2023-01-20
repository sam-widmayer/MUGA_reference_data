# load required packages
library(data.table)
library(qtl2convert)
library(vroom)
library(dplyr)
library(tidyr)

# read in GM marker metadata
gm_ann <- vroom::vroom("data/GigaMUGA/gm_uwisc_v2.csv")

# read in "old" GM allele codes
allele_codes <- vroom::vroom("/projects/compsci/vmp/USERS/widmas/lcGBS_wf/data/GM/GM_allelecodes.csv", skip = 3)

# read in QC'd founder consensus genotypes
founder_genos_path <- "output/GigaMUGA/GigaMUGA_founder_consensus_genotypes.csv.gz"
system(paste("gunzip", founder_genos_path))
founder_genos <- vroom::vroom("output/GigaMUGA/GigaMUGA_founder_consensus_genotypes.csv")

# Restrict annotation files to markers with HQ consensus genotypes
# and convert SNP genotypes to allele codes
allele_codes <- gm_ann[which(gm_ann$marker %in% founder_genos$marker),] %>%
  dplyr::select(marker, chr,  snp) %>%
  tidyr::separate(snp, into = c("A","B"), sep = -1)

write.csv(allele_codes, "output/GigaMUGA/GM_allele_codes_GRCm39.csv")
qtl2convert::write2csv(df = allele_codes, 
                       filename = "output/GigaMUGA/GM_allele_codes_GRCm39.csv", 
                       comment = "Allele codes for GigaMUGA",
                       overwrite = T)
