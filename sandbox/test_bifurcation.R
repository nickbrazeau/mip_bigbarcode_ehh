source("R/01-EHH_tools.R")
source("R/02-spiderplot.R")
source("analyses/00_metadata.R")
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
library(rehh)

vcf = vcfR::read.vcfR(file = "~/Desktop/playvcf/crt_ehh_play.polarized.vcf")
playdf <- drugres %>%
  dplyr::filter(name == "crt")

pmap <- getpmap.polarize(vcf, chrommapdf = mapmodels )
thap <-  vcfRmanip::vcfR2thap(vcf)

outdir <- "~/Desktop/test_spiderplot/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

write.table(x = thap, file = paste0(outdir, "testspiderplot.thap"),
              sep = " ", quote = F, row.names =F, col.names = F)

pmap <- pmap %>%
  dplyr::mutate(pos = as.numeric( stringr::str_split_fixed(snpname, "_", n=4)[,4] ))
write.table(x = pmap, file = paste0(outdir, "testspiderplot.inp"),
              sep = " ", quote = F, row.names =F, col.names = F)

haplohh <- rehh::data2haplohh(hap_file="~/Desktop/test_spiderplot/testspiderplot.thap",
                              map_file="~/Desktop/test_spiderplot/testspiderplot.inp",
                              haplotype.in.columns=TRUE,
                              recode.allele = T,
                              min_perc_geno.hap = 0,
                              min_perc_geno.snp = 0,
                              min_maf = 0)




spiderplot(     hh = haplohh,
                pmapobj = pmap,
                focal = 400653,
                nucleotides = T,
                left = 22,
                right = 25,
                max.haps = 2,
                palette = "RdBu",
                reverse = FALSE,
                relabel = NULL     )






