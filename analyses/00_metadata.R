#.................................
# Dependencies
#.................................
devtools::install_github("nickbrazeau/rplasmodium")
library(rplasmodium)
library(vcfR)
library(tidyverse)
library(stringr)
library(GenomicRanges)

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
  dplyr::select(-c("chr_fct", "chr"))


#.................................
# alter drug res
#.................................
drugres <- drugres %>%
  dplyr::filter(chr != "Pf_M76611")  # no recombo in mtdna, no EHH


## Make final drug res table
## after talking to group 100 kb was selected
drugres <- drugres %>%
  dplyr::mutate(chrom_fct = rplasmodium::factor_chrom(chr),
                seqname = paste0("chr", as.character(chrom_fct)),
                start = start - 5e4, # already made 1-based
                # TODO fix this
                end = end  + 5e5)
# remember for the vcfR2SubsetChromPos function to work the chromosome name in the VCF must be in the `seqname` field. To protect against multiple chrom names floating around
drugres$seqname <- drugres$chr # updated for our work

# make sure we are still on the genomic map
drugres$start[ min(drugres$start) < 0 ] <- 0
drugends <- aggregate(drugres$end, list(factor(drugres$chr)), max)
colnames(drugends)[1] <- "CHROM"
chromends <- tibble(chr = names( rplasmodium::chromsizes_3d7()),
                    chromend =  rplasmodium::chromsizes_3d7())
( drugends <- chromends %>%
    dplyr::rename(CHROM = chr) %>%
    dplyr::left_join(drugends, .) %>%
  dplyr::mutate(offend = x > chromend) )
if(any(drugends$offend)){
  warning("Have mapped beyone the end of the chromosome; automatic fix")

  drugres <- drugres %>%
    dplyr::left_join(x=., y=chromends, by =  "chr") %>%
    dplyr::mutate(end = ifelse(end > chromend, chromend, end)) %>%
    dplyr::select(-c("chromend", "x"))
}

# but false for crt and dhps so OK (on chrom 7 and 8)




