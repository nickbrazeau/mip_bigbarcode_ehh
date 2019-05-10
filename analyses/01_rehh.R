
#-------------------------------------------------------------------------------------------------------
# Purpose of this script is to perform EHH and iHS based on the ancestral and derived allele
# First we will examine loci as a whole without respect to population structure
# Second we will examine loci relationships with respect to structure (compare pop script)
# In this basic script, we keep all samples and all sites. Later we will filter samples that have missing data
# and do a more targeted approach by focusing on putative drug sites that we found
#-------------------------------------------------------------------------------------------------------

#.................................
# imports
#.................................
source("analyses/00_metadata.R")
source("R/01-EHH_tools.R")
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
library(rehh)


mipbb_dr_panel_vcfR <- vcfR::read.vcfR(file = "data/derived/vcfs/polarized_mipbi_drugres_bigbarcode_panels.vcf.bgz")
mipbb_dr_panel_vcfR@gt[mipbb_dr_panel_vcfR@gt == "./.:.:."] <- NA # import issue. See "issue#1" in vcfdo https://github.com/IDEELResearch/vcfdo/issues/1

# A437G is "properly" coded with the ancestral allele being the alternative "C" base that confers drug sensitivity
mipbb_dr_panel_vcfR@fix[mipbb_dr_panel_vcfR@fix[,"CHROM"] == "Pf3D7_08_v3" & mipbb_dr_panel_vcfR@fix[,"POS"] == "549685",]

#.................................
# datawrangle
#.................................
#..............
# set up drug res sites
#..............
drugreslist <- split(drugres, factor(1:nrow(drugres)))

#..............
# remove sites that lack polarization
#..............

nopolar <- getpmap.polarize(mipbb_dr_panel_vcfR) %>%
  dplyr::filter(is.na(anc)) %>%
  dplyr::mutate(
    CHROM = stringr::str_split_fixed(snpname, "_(?!.*_)", n=2)[,1],
    POS = stringr::str_split_fixed(snpname, "_(?!.*_)", n=2)[,2],
    POS = as.numeric(POS)
  ) %>%
  # needs to be in bed format
  dplyr::rename(chr = CHROM) %>%
  dplyr::mutate(
    start = POS,
    end = POS,
    gene_symbol = ".",
    score = ".",
    strand = ".",
    geneid = snpname,
    seqname = chr
  ) %>%
  dplyr::select(c("chr", "start", "end", "gene_symbol", "score", "geneid", "seqname"))

mipbb_dr_panel_vcfR_ancpolar <- vcfRmanip::vcffilter_ChromPos(vcfRobject = mipbb_dr_panel_vcfR, chromposbed = nopolar)


#..............
# make drug-res gene objects
#..............
drugres_ret_full <- drugres %>%
  select(-c("strand", "chrom_fct", "seqname"))
drugres_ret_full$vcfRobj <- purrr::map(drugreslist, vcfRmanip::vcfR2SubsetChromPos, vcfRobject = mipbb_dr_panel_vcfR_ancpolar)
drugres_ret_full$nvar <- unlist(purrr::map(drugres_ret_full$vcfRobj, function(x){ return(nrow(x@gt)) }))

drugres_ret_sub <- drugres_ret_full %>%
  filter(nvar > 10)


#-------------------------------------------------------------------------------------------------------
# LOCI EHH
#-------------------------------------------------------------------------------------------------------

#.................................
# make map and haplotype files
#.................................
drugres_ret_sub$pmap <- purrr::map(drugres_ret_sub$vcfRobj, getpmap.polarize)
drugres_ret_sub$thap <- purrr::map(drugres_ret_sub$vcfRobj, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/ihs/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(drugres_ret_sub)){

  file = paste0(outdir, drugres_ret_sub$name[i])

  write.table(x = drugres_ret_sub$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =F, col.names = F)

  write.table(x = drugres_ret_sub$pmap[[i]], file = paste0(file, ".inp"),
              sep = " ", quote = F, row.names =F, col.names = F)

}


#.................................
# pull in the hapobj
#.................................
rehhfile <- tibble(
  thap_files = list.files(outdir, pattern = ".thap", full.names = T),
  pmap_files = list.files(outdir, pattern = ".inp", full.names = T),
  name = stringr::str_replace(basename(thap_files), ".thap", ""))

rehhfile$haplohh <- purrr::map2(rehhfile$thap_files, rehhfile$pmap_files,
                                ~rehh::data2haplohh(hap_file=.x, map_file=.y,
                                                    haplotype.in.columns=TRUE,
                                                    recode.allele = T,
                                                    min_perc_geno.hap=0,
                                                    min_perc_geno.snp = 0,
                                                    min_maf = 0)) # assuming Bob's filters already valid

#.................................
# combine drug resistance information with the haplohh obj
#.................................
drugres_ret_sub <- rehhfile %>%
  dplyr::select(c("name", "haplohh")) %>%
  dplyr::left_join(x = drugres_ret_sub, y=., by = "name")

#.................................
# call ihs calcs
#.................................
drugres_ret_sub$scanhh <- purrr::map(drugres_ret_sub$haplohh, rehh::scan_hh,
                                     limhaplo = 2,
                                     limehh = 0.05,
                                     discard_integration_at_border = F)

# looks like NAs are being produced from
# https://github.com/cran/rehh/blob/master/src/hh_utils.c
# line 412 - 426
#   if (discard_integration_at_border && ((y_axis[0] > threshold) || (y_axis[n - 1] > threshold))) {  // If the EHH or EHHS is larger than the minimum value at either end of the chromosome, ...
#   return (UNDEFND);                                                           // ... then do not compute the integral, and quit
# Note, can turn this off by setting `discard_integration_at_border = F`
# Although this setting makes absolute sense in whole genome setting, I have borders at every loci
# by design with MIPs (e.g. really multiplexed targeted amplicon sequencing). This means that I will have hard
# cutoffs for every vcfRobject. As such, I am going to ignore boundaries under the assumption that because the MIP density is all constant
# this gave all regions a "fair chance" to have that same border.
# It will certainly inflate samples/regions integration with



drugres_ret_sub$marker <- purrr::map(drugres_ret_sub$scanhh, function(x){
  maxmarker <- which(x$iES_Sabeti_et_al_2007 == max(x$iES_Sabeti_et_al_2007, na.rm = T))
  maxmarker <- ifelse(purrr::is_empty(maxmarker), NA, maxmarker)
  return(maxmarker)
})

# Throw away loci that don't have a clear marker
drugres_ret_sub <- drugres_ret_sub %>%
  dplyr::filter(!is.na(marker))

#.................................
# call ehh based on ihs value
#.................................
drugres_ret_sub$ehh <- purrr::map2(drugres_ret_sub$haplohh, drugres_ret_sub$marker, ~calc_ehh(
  haplohh = .x,
  mrk = .y,
  limhaplo = 2,
  limehh = 0.05,
  discard_integration_at_border = F,
  plotehh = F))


#.................................
# make plots for ehh continous
#.................................
ehhplotdf <- tibble(geneid = drugres_ret_sub$geneid,
                    marker = unlist( drugres_ret_sub$marker ))


ehhplotdf$pos <- map(drugres_ret_sub$haplohh, "position")
ehhplotdf$ehh <- map(drugres_ret_sub$ehh, "ehh")
ehhplotdf$ehh <- map(ehhplotdf$ehh, function(x){ return(as.data.frame(t(x))) } )

ehhplotdf$ehhdf <- map2(ehhplotdf$pos, ehhplotdf$ehh, ~data.frame(pos = .x,
                                                                  ancestral = .y[,1],
                                                                  derived = .y[,2]))

ehhplotdf$ehhdf <- map(ehhplotdf$ehhdf, function(x){
  ret <- gather(x, key = "allele", value = "ehh", 2:3)
  return(ret)
})


plotehh <- function(ehhdf, marker, geneid){
  plotObj <- ggplot() +
    geom_line(data=ehhdf, aes(x=pos, y=ehh, colour = factor(allele), group = factor(allele))) +
    scale_color_manual("Allele Status", values = c("#4575b4", "#d73027")) +
    ggtitle( paste("EHH for Marker:", marker, "on gene:", geneid) ) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}


ehhplotdf$ehhplot <-  pmap(ehhplotdf[,c("ehhdf", "marker", "geneid")], plotehh)



save(drugres_ret_full,
     drugres_ret_sub,
     ehhplotdf,
     file = "data/derived/01-drugres_obj_rehh_countrylevel.rda")

