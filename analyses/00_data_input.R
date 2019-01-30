#.................................
# Dependencies
#.................................
devtools::install_github("nickbrazeau/rplasmodium")
library(rplasmodium)
library(tidyverse)
library(stringr)

#.................................
# imports
#.................................
mc <- readRDS("analyses/data/monoclonal_samples.rds")
mc$name_short <- stringr::str_split_fixed(mc$name, pattern = "-", n=2)[,1] #!!!! need to fix this with OJ

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
  dplyr::filter(chr != "Pf_M76611") %>%  # no recombo in mtdna, no EHH
  dplyr::mutate(chrom_fct = rplasmodium::factor_chrom(chr),
                seqname = paste0("chr", as.character(chrom_fct)),
                start = start - 5*1e4, # already made 1-based
                end = end + 5*1e4)
# make sure we are still on the genomic map
# min(drugres$start)
# aggregate(drugres$end, list(factor(drugres$chr)), max)
# rplasmodium::chromsizes_3d7()






# IMPORTANT... only need to run this once per go but make a better system call
# #.................................
# # make vcfR
# #.................................
# mipbivcfR <- MIPanalyzerbi2vcfR(input = mipbi, cutoff = 0.1)
# # subset to monoclonal samples
#
# mipbivcfR <- mipbivcfR[,colnames(mipbivcfR@gt) %in% c("FORMAT", mc$name_short)]
# #!!!! LOSING 7 samples from OJ. Ask him if he can use mipmapper sample names
#
#
# #.................................
# # write out to be compatible with pf3d7
# #.................................
# liftover <- tibble(chrom = rplasmodium::chromnames(genome = "pf3d7")[1:14],
#        chr = paste0("chr", seq(1:14)))
# CHROM <- left_join(tibble(chr = vcfR::getCHROM(mipbivcfR)), liftover)
# mipbivcfR@fix[,1] <- unlist(CHROM[,2])
# write.vcf(mipbivcfR, file = "~/Documents/MountPoints/mountedMeshnick/Projects/mip_bigbarcode_ehh/analyses/data/mipbi_bigbarcode.vcf.gz")
#
# system("bash ~/Documents/MountPoints/mountedMeshnick/Projects/mip_bigbarcode_ehh/analyses/data/polarize.sh")




mipbivcfR <- vcfR::read.vcfR("~/Documents/MountPoints/mountedMeshnick/Projects/mip_bigbarcode_ehh/analyses/data/polarized_mipbi_bigbarcode.vcf.gz")


liftover <- tibble(chrom = rplasmodium::chromnames(genome = "pf3d7")[1:14],
        chr = paste0("chr", seq(1:14)))
CHROM <- left_join(tibble(chrom = vcfR::getCHROM(mipbivcfR)), liftover)
mipbivcfR@fix[,1] <- unlist(CHROM[,2])
