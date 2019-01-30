#.................................
# Dependencies
#.................................
devtools::install_github("mrc-ide/MIPanalyzer")
devtools::install_github("nickbrazeau/rplasmodium")
library(tidyverse)
library(MIPanalyzer)
library(rplasmodium)
library(vcfR)
library(rehh)


#.................................
# imports
#.................................
mc <- readRDS("analyses/data/monoclonal_samples.rds")
mipbi <- readRDS("analyses/data/biallelic_processed.rds")
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
                chr_orig = chr,
                CHROM = paste0("chr", as.character(chr_fct))) %>%
  select(-c("chr_fct", "chr"))


#.................................
# alter drug res
#.................................
drugres <- drugres %>%
  dplyr::filter(chr != "Pf_M76611") %>%  # no recombo in mtdna
  dplyr::mutate(chrom_fct = rplasmodium::factor_chrom(chr),
                seqname = paste0("chr", as.character(chrom_fct)),
                start = start - 5*1e4, # already made 1-based
                end = end + 5*1e4)
# make sure we are still on the genomic map
# min(drugres$start)
# aggregate(drugres$end, list(factor(drugres$chr)), max)
# rplasmodium::chromsizes_3d7()

#.................................
# make vcfR
#.................................
mipbivcfR <- MIPanalyzerbi2vcfR(input = mipbi, cutoff = 0.1)
# subset to monoclonal samples
mc$name_short <- stringr::str_split_fixed(mc$name, pattern = "-", n=2)[,1]
mipbivcfR <- mipbivcfR[,colnames(mipbivcfR@gt)[2:ncol(mipbivcfR@gt)] %in% mc$name_short]
#!!!! LOSING 7 samples from OJ. Ask him if he can use mipmapper sample names
write.vcf(mipbivcfR, file = "~/Documents/MountPoints/mountedMeshnick/Projects/mip_ehh/analyses/data/mipbi_bigbarcode.vcf.gz")


