#.................................
# Dependencies
#.................................
devtools::install_github("nickbrazeau/rplasmodium")
library(rplasmodium)
library(vcfR)
library(tidyverse)
library(stringr)


#--------------------------------------------------------------
# Metadata files
#--------------------------------------------------------------

# monoclonal samples
mc <- readRDS("data/monoclonals.rds")
# drug res sites
drugres <- rplasmodium::pf_3d7_PutDrugRxSites
drugres$start <- drugres$start + 1 # take 1-based, vcfs are 1-based

#.................................
# make map file
#.................................
maplm <- function(data){
  lm(cM ~ pos,
     data = data)
}
mp <- rplasmodium::recomb_map_pf3d7
mapmodels <- mp %>%
  group_by(chr) %>%
  nest() %>%
  dplyr::mutate(model = purrr::map(data, maplm),
                chr_fct = rplasmodium::factor_chrom(chr),
                CHROM = chr
                ) %>%
  select(-c("chr_fct", "chr"))


#.................................
# alter drug res
#.................................
drugres <- drugres %>%
  dplyr::filter(chr != "Pf_M76611") %>%  # no recombo in mtdna, no EHH
  dplyr::mutate(chrom_fct = rplasmodium::factor_chrom(chr),
                seqname = paste0("chr", as.character(chrom_fct)),
                start = start - 5*1e4, # already made 1-based
                end = end + 5*1e4)
# make sure we are still on the genomic map
min(drugres$start)
aggregate(drugres$end, list(factor(drugres$chr)), max)
rplasmodium::chromsizes_3d7()
