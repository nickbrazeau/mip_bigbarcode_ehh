
#-------------------------------------------------------------------------------------------------------
# Purpose of this script is to perform EHH and iHS based on the ancestral and derived allele
# note this is a different question than cross-population comparisons
#-------------------------------------------------------------------------------------------------------

#.................................
# imports
#.................................
source("R/00_functions_temp_if_PR.R")
source("analyses/00_data_input.R")
devtools::install_github("IDEELResearch/vcfR2plsmdmcalc")
library(vcfR2plsmdmcalc)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)


#.................................
# datawrangle
#.................................

drugreslist <- split(drugres, factor(1:nrow(drugres)))

drugres_ret_full <- drugres %>%
  select(-c("strand", "chrom_fct", "seqname"))
drugres_ret_full$vcfRobj <- purrr::map(drugreslist, vcfRmanip::vcfR2SubsetChromPos, vcfRobject = mipbivcfR)
drugres_ret_full$nvar <- unlist(purrr::map(drugres_ret_full$vcfRobj, function(x){ return(nrow(x@gt)) }))

drugres_ret_sub <- drugres_ret_full %>%
  filter(nvar > 20)

#.................................
# make map and haplotype files
#.................................

getpmap.polarize <- function(x, chrommapdf = mapmodels){
  chrompos <- vcfR::getFIX(x)[, c("CHROM", "POS")] %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(pos = as.numeric(POS))

  predmp <- inner_join(chrommapdf, chrompos) %>%
    mutate(cMpred = map2_dbl(model, pos, ~predict(.x, newdata = data.frame(pos = .y))))

  fixtidy <- vcfR::extract_info_tidy(x)
  refalt <- tibble::as_tibble(getFIX(x)[,c("REF", "ALT")])
  alleles <- dplyr::bind_cols(refalt, fixtidy)

  ancderiv <- tibble::tibble(
    anc = alleles$AA,
    der = ifelse(alleles$REF == alleles$AA,
                 alleles$ALT, alleles$REF)

  )
  # APM specific masks. Shouldn't be an issue in genic regions
  ancderiv$anc[fixtidy$AA %in% c("N", "X")] <- NA
  ancderiv$der[fixtidy$AA %in% c("N", "X")] <- NA

  ret <- tibble::tibble(snpname = paste0(predmp$chr_orig, "_", predmp$pos),
                        chrom = predmp$CHROM,
                        pos = predmp$cMpred,
                        anc = ancderiv$anc,
                        der = ancderiv$der
  )

  return(ret)
}

drugres_ret_sub$pmap <- purrr::map(drugres_ret_sub$vcfRobj, getpmap.polarize)
drugres_ret_sub$thap <- purrr::map(drugres_ret_sub$vcfRobj, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/ihs/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(drugres_ret_sub)){

  file = paste0(outdir, drugres_ret_sub$name[i])

  write.table(x = drugres_ret_sub$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =T, col.names = F)

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

rehhfile$haplohh <- purrr::map2(rehhfile$thap_files, rehhfile$pmap_files,  ~data2haplohh(hap_file=.x, map_file=.y,
                                                                                         haplotype.in.columns=TRUE, recode.allele = T,
                                                                                         min_perc_geno.hap=0,
                                                                                         min_perc_geno.snp=0)) # I don't want it's on the fly drops


#.................................
# combine
#.................................

drugres_ret_sub <- rehhfile %>%
  dplyr::select(c("name", "haplohh")) %>%
  dplyr::left_join(x = drugres_ret_sub, y=., by = "name")

#.................................
# calculate ihh
#.................................
drugres_ret_sub$scanhh <- purrr::map(drugres_ret_sub$haplohh, scan_hh,
                                     limhaplo = 2,
                                     limehh = 0.00)

drugres_ret_sub$marker <- purrr::map(drugres_ret_sub$scanhh, function(x){
  maxmarker <- which(x$iES_Sabeti_et_al_2007 == max(x$iES_Sabeti_et_al_2007, na.rm = T))
  return(maxmarker)
})

# !!!!! temporary -- something odd going on with mdr1
drugres_ret_sub$marker[[2]] <- 19




