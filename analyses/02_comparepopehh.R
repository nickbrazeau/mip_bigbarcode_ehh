
#-------------------------------------------------------------------------------------------------------
# Purpose of this script is to perform EHH controlling for population structure
#-------------------------------------------------------------------------------------------------------

#.................................
# imports
#.................................
source("analyses/00_metadata.R")
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
load("data/derived/drug_res_nopopstructure.rda")

mipbb_dr_panel_vcfR <- vcfR::read.vcfR(file = "data/polarized_mipbi_bigbarcode.vcf.bgz")


#.................................
# need drug and population specific sites and combine
#.................................
regions <- mc %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(region = kmc2) %>%
  group_by(region) %>%
  nest()
regions$smpls <- purrr::map(regions$data, "id")

regionlong <- do.call("rbind", lapply(1:nrow(drugres), function(x){
  return(regions)
})) %>%
  dplyr::arrange(region)

drugregions <- as_tibble(cbind.data.frame(regionlong, drugres))

#.................................
# subset and expand vcf
#.................................

subsetlocismpls <- function(smpls, seqname, start, end, geneid){

  chromposbed = tibble(seqname = seqname, start = start, end = end, geneid = geneid)
  ret <- vcfRmanip::vcfR2SubsetChromPos(vcfRobject = mipbb_dr_panel_vcfR, # global
                                        chromposbed = chromposbed) # subset loci
  ret <- ret[, colnames(ret@gt) %in% c("FORMAT", smpls)] # subset samples
  return(ret)
}

drugregions$vcfRobj <- purrr::pmap(drugregions[,c("smpls", "seqname", "start", "end", "geneid")], subsetlocismpls)


#.................................
# remove bad sites
#.................................
drugregions$nvar <- unlist(purrr::map(drugregions$vcfRobj, function(x){ return(nrow(x@gt)) }))
drugregions$nsmpls <- unlist(purrr::map(drugregions$vcfRobj, function(x){ return(ncol(x@gt)-1) }))
drugregions_sub <- drugregions %>%
  filter(nvar > 10) %>%
  filter(nsmpls > 5)


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
  refalt <- tibble::as_tibble(vcfR::getFIX(x)[,c("REF", "ALT")])
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

drugregions_sub$pmap <- purrr::map(drugregions_sub$vcfRobj, getpmap.polarize)
drugregions_sub$thap <- purrr::map(drugregions_sub$vcfRobj, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/ehhcrosspop/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(drugregions_sub)){

  file = paste0(outdir, drugregions_sub$name[i], "-", drugregions_sub$region[i])

  write.table(x = drugregions_sub$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =T, col.names = F)

  write.table(x = drugregions_sub$pmap[[i]], file = paste0(file, ".inp"),
              sep = " ", quote = F, row.names =F, col.names = F)

}


#.................................
# pull in the hapobj
#.................................
rehhfile <- tibble(
  thap_files = list.files(outdir, pattern = ".thap", full.names = T),
  pmap_files = list.files(outdir, pattern = ".inp", full.names = T),
  name = stringr::str_split_fixed(stringr::str_replace(basename(thap_files), ".thap", ""), "-", n=2)[,1],
  region = stringr::str_split_fixed(stringr::str_replace(basename(thap_files), ".thap", ""), "-", n=2)[,2]
)

rehhfile$haplohh <- purrr::map2(rehhfile$thap_files, rehhfile$pmap_files,  ~rehh::data2haplohh(hap_file=.x, map_file=.y,
                                                                                         haplotype.in.columns=TRUE,
                                                                                         recode.allele = T,
                                                                                         min_perc_geno.hap = 50,
                                                                                         min_perc_geno.snp = 50,
                                                                                         min_maf = 0)) # turn off their cutoffs


#.................................
# pull it all together
#.................................

drugregions_sub <- rehhfile %>%
  dplyr::select(c("name", "region", "haplohh")) %>%
  dplyr::left_join(x = drugregions_sub, y=., by = c("name", "region"))


#.................................
# call ihs calcs
#.................................
drugregions_sub$scanhh <- purrr::map(drugregions_sub$haplohh, rehh::scan_hh,
                                     limhaplo = 2,
                                     limehh = 0.05)



# manually set for now
markerdf <- tibble::tibble(
  name = c("kelch", "crt", "mdr1", "dhps", "dhfr", "pfabcI3", "pfpare"),
  mutation = c("kelch-misc", "crt-N75E", "mdr1-N86Y", "dhps-K540E", "dhfr-misc",  "pfabcI3-misc", "pfpare-?-strong"),
  marker = c(15, 11, 14, 16, 5, 6, 48)
)
drugregions_sub <- left_join(drugregions_sub, markerdf, by = "name")


#.................................
# call ehh based on ihs value
#.................................
drugregions_sub$ehh <- map2(drugregions_sub$haplohh, drugregions_sub$marker, ~rehh::calc_ehh(
                            haplohh = .x,
                            mrk = .y,
                            limhaplo = 2,
                            limehh = 0.05,
                            plotehh = F)
)



drugregions_sub$pos <- map(drugregions_sub$haplohh, "position")
drugregions_sub$ehh <- map(drugregions_sub$ehh, "ehh")
drugregions_sub$ehh <- map(drugregions_sub$ehh, function(x){ return(as.data.frame(t(x))) } )

drugregions_sub$ehhdf <- map2(drugregions_sub$pos, drugregions_sub$ehh, ~data.frame(pos = .x,
                                                                                    ancestral = .y[,1],
                                                                                    derived = .y[,2]))
drugregions_sub$ehhdf <- map(drugregions_sub$ehhdf, function(x){
  ret <- gather(x, key = "allele", value = "ehh", 2:3)
  return(ret)
})

#.................................
# make plots for ehh continous
#.................................
plotehh <- function(ehhdf, region, mutation){
  plotObj <- ggplot() +
    geom_line(data=ehhdf, aes(x=pos, y=ehh, colour = factor(allele), group = factor(allele))) +
    scale_color_manual("Allele Status", values = c("#4575b4", "#d73027")) +
    ggtitle(paste0(mutation, " for region ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

drugregions_sub$ehhplot <-  pmap(drugregions_sub[,c("region", "mutation", "ehhdf")], plotehh)

#.................................
# make haplotype plots
#.................................
drugregions_sub$hapdf <- map(drugregions_sub$thap, function(x){
  x <- t(x)
  x <- as.data.frame(cbind(paste0("smpl", seq(1:nrow(x))), x))
  colnames(x) <- c("smpl", paste0("loci", seq(1:(ncol(x)-1))))
  ret <- gather(x, key = "loci", value = "allele", 2:ncol(x))
  return(ret)
})


plothap <- function(hapdf, region, name){
  plotObj <- ggplot() +
    geom_tile(data=hapdf, aes(x=loci, y=smpl, fill = factor(allele))) +
    scale_fill_manual("Allele", values = c("#d73027", "#fee090", "#66bd63", "#4575b4")) +
    ggtitle(paste0("Gene ", name, " for region ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.text = element_blank())
}
drugregions_sub$happlot <-  pmap(drugregions_sub[,c("region", "name", "hapdf")], plothap)



# TODO make bifurcation plots (see rehh)
# TODO in above plots -- highlight with grey bar the focal snp

#example haplohh object (280 haplotypes, 1424 SNPs)
#see ?haplohh_cgu_bta12 for details
data(haplohh_cgu_bta12)
#plotting bifurcation diagram for both ancestral and derived allele
#from the focal SNP at position 456
#which display a strong signal of selection
layout(matrix(1:2,2,1))
#ancestral allele
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc=456,all_foc=1,
                    nmrk_l=20,nmrk_r=20)
#derived allele
bifurcation.diagram(haplohh_cgu_bta12,mrk_foc=456,all_foc=2,
                    nmrk_l=20,nmrk_r=20)
##
dev.off()



#.................................
# call plots
#.................................
library(gridExtra)

final <- drugregions_sub %>%
  group_by(name) %>%
  nest()
final$ehhplot <- map(final$data, "ehhplot")
final$happlot <- map(final$data, "happlot")

map2(final$ehhplot, final$happlot, function(x, y){
  n <- length(x)
  nCol <- floor(sqrt(n))
  do.call("grid.arrange", c(x, ncol=nCol))

  m <- length(y)
  mCol <- floor(sqrt(m))
  do.call("grid.arrange", c(y, ncol=mCol))

})


#.................................
# write out
#.................................
save(drugregions_sub, file = "data/derived/drugregions_sub.rda")
