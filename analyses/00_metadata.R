#.................................
# Dependencies
#.................................
devtools::install_github("andrewparkermorgan/rplasmodium", ref = "bc0801326ddf7ab99f961a8b9d58416baad9a258") # until PR is accepted
library(rplasmodium)
library(vcfR)
library(tidyverse)
library(stringr)
library(fuzzyjoin)

#--------------------------------------------------------------
# Metadata files
#--------------------------------------------------------------

#.................................
# Monoclonal Samples
#.................................
mc <- readRDS("data/monoclonal_samples.rds")


#.................................
# Recombination Map File
#.................................
rmap <- rplasmodium::recomb_map_pf3d7
# join ends of chromosomes for boundaries
mapmodels <- rmap %>% dplyr::bind_rows(tibble::tibble(
          chr = names(rplasmodium::chromsizes()),
           pos = rplasmodium::chromsizes(),
           marker = "end",
           cM  = NA)) %>%
  dplyr::mutate(chr = rplasmodium::factor_chrom(chr)) %>%
  dplyr::filter(! chr %in% c("M", "A")) # no recombination in mtdna or apicoplast, no EHH

mapmodels <- mapmodels %>%
  dplyr::group_by(chr) %>%
  dplyr::mutate(lag_cM = dplyr::lag(cM),
                lag_pos = dplyr::lag(pos),

                lag_pos = ifelse(pos == min(pos), 0, lag_pos), # set up boundaries for min
                lag_cM = ifelse(pos == min(pos), 0, lag_cM) # set up boundaries for min
  ) %>%
  dplyr::mutate(cM_slope = (lag_cM - cM)/ (lag_pos - pos),
                cM_intercept = lag_cM ) %>%
  dplyr::group_by(chr) %>%
  dplyr::mutate(
               # cM_slope = ifelse(pos == min(pos), dplyr::lead(cM_slope), cM_slope), # this created to fast of a recombination rate, will interpolate instead
                cM_slope = ifelse(pos == max(pos), dplyr::lag(cM_slope), cM_slope), # fix boundaries
                cM_intercept = ifelse(pos == min(pos), 0, cM_intercept), # fix boundaries
                cM_intercept = ifelse(pos == max(pos), lag_cM, cM_intercept) # fix boundaries
                ) %>%
  # clean up
  dplyr::select(c("chr", "lag_pos", "pos", "cM_slope", "cM_intercept")) %>%
  dplyr::rename(start = lag_pos,
                end = pos) %>%
  dplyr::mutate(start = start + 1) # was doing inference on a 0-based, inclusive boundary system. Now add 1 to make it non-overlapping/1-based to match with VCF





#.................................
# Drug Res Sites
#.................................

drugres <- rplasmodium::pf_3d7_PutDrugRxSites
# TODO temporary addition until PR request is accepted

temp <- tibble(chr = "Pf3D7_01_v3",
               start = 265207,
               end = 269173,
               gene_symbol = "pfatp6",
               score = ".",
               strand = "-",
               gene_id = "PF3D7_0106300")
drugres <- dplyr::bind_rows(drugres, temp)


drugres$start <- drugres$start + 1 # take to 1-based, vcfs are 1-based (in bed format)



drugres <- drugres %>%
  dplyr::filter(chr != "Pf_M76611")  # no recombo in mtdna, no EHH


## Make final drug res table
# after talking to group, we decided for 100kb
drugres <- drugres %>%
  dplyr::mutate(chrom_fct = rplasmodium::factor_chrom(chr),
                seqname = paste0("chr", as.character(chrom_fct)),
                start = start - 2e5, # already made 1-based
                end = end  + 2e5)
# remember for the vcfR2SubsetChromPos function to work the chromosome name in the VCF must be in the `seqname` field. To protect against multiple chrom names floating around
drugres$seqname <- drugres$chr # updated for our work

# make sure we are still on the genomic map

if(any(drugres$start < 0)){
  warning("Have mapped beyone the start of the chromosome; automatic fix to 0 (beginning of chromosome)")
  drugres$start[ which( drugres$start < 0 ) ] <- 0
}

drugends <- drugres %>%
  dplyr::select(c("chr", "gene_symbol", "end"))
colnames(drugends)[1] <- "CHROM"
chromends <- tibble(chr = names( rplasmodium::chromsizes_3d7()),
                    chromend =  rplasmodium::chromsizes_3d7())
( drugends <- chromends %>%
    dplyr::rename(CHROM = chr) %>%
    dplyr::left_join(drugends, .) %>%
  dplyr::mutate(offend = end > chromend) )
if(any(drugends$offend)){
  warning("Have mapped beyone the end of the chromosome; automatic fix to end of chromosome")

  drugres <- drugres %>%
    dplyr::left_join(x=., y=chromends, by =  "chr") %>%
    dplyr::mutate(end = ifelse(end > chromend, chromend, end)) %>%
    dplyr::select(-c("chromend"))
}

# rename drug res gene_id to geneid for compatibility with vcfRmanip
drugres <- drugres %>%
  dplyr::rename(geneid = gene_id,
                name = gene_symbol) # for legacy reasons

# clean up
rm(drugends, temp)

