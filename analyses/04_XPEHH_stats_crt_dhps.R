#-------------------------------------------------------------------------------------------------------
#
# We are going to subset to just look at the
# crt loci 76T and the dhps 540 and 581 loci for
# cross population statistics
#
#-------------------------------------------------------------------------------------------------------
#.................................
# imports
#.................................
source("R/01-EHH_tools.R")
source("R/02-spiderplot.R")
source("R/03-crude_iES.R")
source("analyses/00_metadata.R")
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
library(rehh)
set.seed(42)


load("data/derived/03-analyze_crt_dhps_drc.rda")

#.................................
# subset relevant tasks
#.................................
crtdhps_xpehh <- crtdhps_sub %>%
  dplyr::select(c("region", "name", "mut_name", "marker", "haplohh")) %>%
  dplyr::filter(mut_name %in% c("K76T", "K540E", "A581G")) %>%
  dplyr::full_join(., ., by = c("name", "mut_name", "marker")) %>% # make cross-compare
  dplyr::filter(region.x != region.y)

# make map df for pmap
crtdhps_xpehh_map <- crtdhps_xpehh %>%
  dplyr::select("haplohh.x", "haplohh.y", "marker") %>%
  dplyr::mutate(discard_integration_at_border = F) %>%
  dplyr::rename(haplohh.pop1 = haplohh.x,
                haplohh.pop2 = haplohh.y)

crtdhps_xpehh$crudexpehh <- unlist( purrr::pmap(crtdhps_xpehh_map, unnormalized_xpehh_derivedallele_logratio) )

# store results
crtdhps_xpehh_dist <- crtdhps_xpehh %>%
  dplyr::select(c("region.x", "region.y", "crudexpehh", "mut_name")) %>%
  dplyr::mutate(crudexpehh = round(crudexpehh, 2)) %>%
  group_by(mut_name) %>%
  nest()
# coerce to dist matrix for better viz
crtdhps_xpehh_dist$dist <- map(crtdhps_xpehh_dist$data, pairwisedistdf_to_squaredist)
# write out
for(i in 1:nrow(crtdhps_xpehh_dist)){
  readr::write_csv(x = as.data.frame( as.matrix( crtdhps_xpehh_dist$dist[[i]] ) ),
                   path = paste0("results/final_tables/xpehh_loci_", crtdhps_xpehh_dist$mut_name[[i]], "_full_sampleset.csv"))

}

#-------------------------------------------------------------------------------------------------------------------------------------
# Down-Sampling Section
#-------------------------------------------------------------------------------------------------------------------------------------
#.................................
# Find Downsampling n
#.................................
crtdhps_sub %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    nsmpls = mean(nsmpls)
  )

# DRC1.2 has 65 while DRC2.2 has 42. Ghana, Uganda, Zambia have 15, 11, 9 respetively
# 12 seems resonable
crtdhps_xpehh_dwnsmpl <- crtdhps_sub %>%
  dplyr::select(c("region", "name", "mut_name", "marker", "vcfRobj"))

# find 12 samples from DRC-E and DRC-W
drce <-  crtdhps_xpehh_dwnsmpl[crtdhps_xpehh_dwnsmpl$region == "DRC1.2" & crtdhps_xpehh_dwnsmpl$mut_name == "M74I", "vcfRobj"] # could have chosen any DRC1.2 to extract samples
drce <- drce[[1]][[1]]
drce <- sample(x = colnames(drce@gt)[2:ncol(drce@gt)], size = 12)

drcw <-  crtdhps_xpehh_dwnsmpl[crtdhps_xpehh_dwnsmpl$region == "DRC2.2" & crtdhps_xpehh_dwnsmpl$mut_name == "M74I", "vcfRobj"] # could have chosen any DRC1.2 to extract samples
drcw <- drcw[[1]][[1]]
drcw <- sample(x = colnames(drcw@gt)[2:ncol(drcw@gt)], size = 12)


#.................................
# make new vcfR obj
#.................................


mknewvcfR <- function(vcfRobj, region){
  switch(region,
         DRC1.2={
           vcfRobj_new <- vcfRmanip::select_samples(vcfRobj, drce)
           return(vcfRobj_new)
         },
         DRC2.2 = {
           vcfRobj_new <- vcfRmanip::select_samples(vcfRobj, drcw)
           return(vcfRobj_new)
         },
         {
           return(vcfRobj)
         })
}

crtdhps_xpehh_dwnsmpl$new_vcfRobj <- pmap(crtdhps_xpehh_dwnsmpl[,c("vcfRobj", "region")], mknewvcfR)

#.................................
# make NEW map and haplotype files
#.................................
crtdhps_xpehh_dwnsmpl$pmap <- purrr::map(crtdhps_xpehh_dwnsmpl$new_vcfRobj, getpmap.polarize)
crtdhps_xpehh_dwnsmpl$thap <- purrr::map(crtdhps_xpehh_dwnsmpl$new_vcfRobj, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/downsample_ehhcrosspop/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(crtdhps_xpehh_dwnsmpl)){

  file = paste0(outdir, crtdhps_xpehh_dwnsmpl$mut_name[i], "-", crtdhps_xpehh_dwnsmpl$region[i])

  write.table(x = crtdhps_xpehh_dwnsmpl$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =F, col.names = F)

  write.table(x = crtdhps_xpehh_dwnsmpl$pmap[[i]], file = paste0(file, ".inp"),
              sep = " ", quote = F, row.names =F, col.names = F)

}


#.................................
# pull in the hapobj
#.................................
rehhfile <- tibble(
  thap_files = list.files(outdir, pattern = ".thap", full.names = T),
  pmap_files = list.files(outdir, pattern = ".inp", full.names = T),
  mut_name = stringr::str_split_fixed(stringr::str_replace(basename(thap_files), ".thap", ""), "-", n=2)[,1],
  region = stringr::str_split_fixed(stringr::str_replace(basename(thap_files), ".thap", ""), "-", n=2)[,2]
)

rehhfile$haplohh <- purrr::map2(rehhfile$thap_files, rehhfile$pmap_files,  ~rehh::data2haplohh(hap_file=.x, map_file=.y,
                                                                                               haplotype.in.columns=TRUE,
                                                                                               recode.allele = T,
                                                                                               min_perc_geno.hap = 0,
                                                                                               min_perc_geno.snp = 0,
                                                                                               min_maf = 0)) # Assuming Bob's filters upstream good enough
#.................................
# pull it all together
#.................................

crtdhps_xpehh_dwnsmpl <- rehhfile %>%
  dplyr::select(c("mut_name", "region", "haplohh")) %>%
  dplyr::left_join(x = crtdhps_xpehh_dwnsmpl, y=., by = c("mut_name", "region"))


#.................................
# subset relevant tasks
#.................................
crtdhps_xpehh_dwnsmpl <- crtdhps_xpehh_dwnsmpl %>%
  dplyr::select(c("region", "name", "mut_name", "marker", "haplohh")) %>%
  dplyr::filter(mut_name %in% c("K76T", "K540E", "A581G")) %>%
  dplyr::full_join(., ., by = c("name", "mut_name", "marker")) %>% # make cross-compare
  dplyr::filter(region.x != region.y)

# make map df for pmap
crtdhps_xpehh_dwnsmpl_map <- crtdhps_xpehh_dwnsmpl %>%
  dplyr::select("haplohh.x", "haplohh.y", "marker") %>%
  dplyr::mutate(discard_integration_at_border = F) %>%
  dplyr::rename(haplohh.pop1 = haplohh.x,
                haplohh.pop2 = haplohh.y)

crtdhps_xpehh_dwnsmpl$crudexpehh <- unlist( purrr::pmap(crtdhps_xpehh_dwnsmpl_map, unnormalized_xpehh_derivedallele_logratio) )

# store results
crtdhps_xpehh_dwnsmpl_dist <- crtdhps_xpehh_dwnsmpl %>%
  dplyr::select(c("region.x", "region.y", "crudexpehh", "mut_name")) %>%
  dplyr::mutate(markerxpehh = round(crudexpehh, 2)) %>%
  group_by(mut_name) %>%
  nest()
# coerce to dist matrix for better viz
crtdhps_xpehh_dwnsmpl_dist$dist <- map(crtdhps_xpehh_dwnsmpl_dist$data, pairwisedistdf_to_squaredist)
# write out
for(i in 1:nrow(crtdhps_xpehh_dwnsmpl_dist)){
  readr::write_csv(x = as.data.frame( as.matrix( crtdhps_xpehh_dwnsmpl_dist$dist[[i]] ) ),
                   path = paste0("results/final_tables/xpehh_loci_", crtdhps_xpehh_dwnsmpl_dist$mut_name[[i]], "_full_sampleset.csv"))

}

