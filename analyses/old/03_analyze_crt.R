#-------------------------------------------------------------------------------------------------------
# After discussing with JJJ, JAB & SRM
# have decided to subset to just the DRC with
# K-means of 2 cluster for just the DRC
#
# We are going to subset to just look at the
# crt loci 76T and the dhpa 540 and 581 loci
#
#-------------------------------------------------------------------------------------------------------
#.................................
# imports
#.................................
source("R/01-EHH_tools.R")
source("R/02-spiderplot.R")
source("analyses/00_metadata.R")
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
set.seed(42)

mipbb_dr_panel_vcfR <- vcfR::read.vcfR(file = "data/polarized_mipbi_drugres_bigbarcode_panels.vcf.bgz")

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
  dplyr::arrange(region) # make a copy of the region data for each drug res site (not memory efficient but fits within long tidy framework)

drugregions <- as_tibble(cbind.data.frame(regionlong, drugres))




#.................................
# need drug and population specific sites and combine
# NOTE, based on group discussions have elected to keep the
# DRC as one large population instead of considering substructure
#.................................
regions <- mc %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::rename(region = country) %>%
  group_by(region) %>%
  nest()
regions$smpls <- purrr::map(regions$data, "id")

# make crt table
crt <- as_tibble(cbind.data.frame(regions, drugres[drugres$name == "crt", ]))

#.................................
# subset vcf to crt
#.................................

subsetlocismpls <- function(smpls, seqname, start, end, geneid){

  chromposbed = tibble(seqname = seqname, start = start, end = end, geneid = geneid)
  ret <- vcfRmanip::vcfR2SubsetChromPos(vcfRobject = mipbb_dr_panel_vcfR, # global
                                        chromposbed = chromposbed) # subset loci
  ret <- ret[, colnames(ret@gt) %in% c("FORMAT", smpls)] # subset samples
  return(ret)
}

crt$vcfRobj <- purrr::pmap(crt[,c("smpls", "seqname", "start", "end", "geneid")], subsetlocismpls)


#.................................
# remove countries with no data
#.................................
crt$nsmpls <- unlist(purrr::map(crt$vcfRobj, function(x){ return(ncol(x@gt)-1) }))
crt <- crt %>%
  filter(nsmpls > 10)

#.................................
# Downsample countries to minimum n
# This way there is no bias in terms of
# how many samples could be derived v. ancestral
#.................................
dwnsmplmin <- min(crt$nsmpls)

Down_Sample_Vcf <- function(vcfRobj, smplmin){
  smpls <- ncol(vcfRobj@gt) - 1
  rsmpls <- sort(sample(x = 1:smpls, size = smplmin))
  rsmpls <- rsmpls + 1 # add one to offset for format column
  vcfRobj@gt <- vcfRobj@gt[, c(1, rsmpls)] # must include Format (column 1)
  return(vcfRobj)
}

crt$vcfRobj_subset <- purrr::map(crt$vcfRobj, Down_Sample_Vcf, smplmin = 14)



#.................................
# make map and haplotype files for crt
#.................................
crt$pmap <- purrr::map(crt$vcfRobj_subset, getpmap.polarize)
crt$thap <- purrr::map(crt$vcfRobj_subset, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/crtehhcrosspop/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(crt)){

  file = paste0(outdir, crt$name[i], "-", crt$region[i])

  write.table(x = crt$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =F, col.names = F)

  write.table(x = crt$pmap[[i]], file = paste0(file, ".inp"),
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
                                                                                               min_perc_geno.hap = 0,
                                                                                               min_perc_geno.snp = 1, # if entire SNP is missing information, it is due to the fact that we weren't able to polarize it (otherwise would have been dropped upstream); need to drop those as that information is useless
                                                                                               min_maf = 0))
#.................................
# pull it all together
#.................................

crt <- rehhfile %>%
  dplyr::select(c("name", "region", "haplohh")) %>%
  dplyr::left_join(x = crt, y=., by = c("name", "region"))


#.................................
# call ihs calcs
#.................................
crt$scanhh <- purrr::map(crt$haplohh, rehh::scan_hh,
                         limhaplo = 2,
                         limehh = 0.05)

# Find the marker with the largest integrated EHHS statistic by Sabeti (unnormalized)
crt$marker <- purrr::map(crt$scanhh, function(x){
  maxmarker <- which(x$iES_Sabeti_et_al_2007 == max(x$iES_Sabeti_et_al_2007, na.rm = T))
  maxmarker <- ifelse(purrr::is_empty(maxmarker), NA, maxmarker)
  return(maxmarker)
})

if(length(unique(crt$marker)) != 1){
  warning("Multiple markers by region are maximized by iES")
}

#.................................
# call ehh based on ihs value
# here by hand because of interest snp
#.................................
# chr7	403621 crt-N75E	mark
N75Emark <- which(crt$haplohh[[1]]@snp.name == "Pf3D7_07_v3_403620")

crt$ehh_N75E_mark <- map(crt$haplohh, ~rehh::calc_ehh(
  haplohh = .x,
  mrk = N75Emark,
  limhaplo = 2,
  limehh = 0.05,
  plotehh = F)
)


#.................................
# make crt ehh plots
#.................................
crt$pos <- map(crt$haplohh, "position")
# N75E MARK
crt$N75E_ehhdf <- map(crt$ehh_N75E_mark, "ehh")
crt$N75E_ehhdf <- map(crt$N75E_ehhdf, function(x){ return(as.data.frame(t(x))) } )
crt$N75E_ehhdf <- map2(crt$pos, crt$N75E_ehhdf, ~data.frame(pos = .x,
                                                            ancestral = .y[,1],
                                                            derived = .y[,2]))

crt$N75E_ehhdf <- map(crt$N75E_ehhdf, function(x){
  ret <- gather(x, key = "allele", value = "ehh", 2:3)
  return(ret)
})



#.................................
# make plots for ehh continous
#.................................
plotehh <- function(ehhdf, region, marker){
  plotObj <- ggplot() +
    geom_line(data=ehhdf, aes(x=pos, y=ehh, colour = factor(allele), group = factor(allele))) +
    scale_color_manual("Allele Status", values = c("#4575b4", "#d73027")) +
    ggtitle(paste0("This is marker: ", marker, "\n for region ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

# make plots
crt$n75e_marker_plot <- mapply(plotehh,
                               ehhdf =  crt$N75E_ehhdf,
                               region = crt$region,
                               marker = crt$haplohh[[1]]@snp.name[ N75Emark ],
                               SIMPLIFY = F)


#.................................
# make bifurcation plot
#.................................
focal = crt$haplohh[[1]]@position[ N75Emark ] # 11 is where

map2(crt$haplohh, crt$pmap, ~spiderplot(hh = .x,
                                        pmapobj = .y,
                                        focal = focal,
                                        nucleotides = T,
                                        left = 10,
                                        right = 10,
                                        max.haps = 2,
                                        palette = "RdBu",
                                        reverse = FALSE,
                                        relabel = NULL))
spiderplot(hh = crt$haplohh[[2]],
           pmapobj = crt$pmap[[2]],
           focal = focal,
           nucleotides = T,
           left = 10,
           right = 10,
           max.haps = 2,
           palette = "RdBu",
           reverse = FALSE,
           relabel = NULL)




#.................................
# make haplotype plots
#.................................
crt$hapdf <- map(crt$thap, function(x){
  x <- t(x)
  x <- as.data.frame(cbind(paste0("smpl", seq(1:nrow(x))), x))
  colnames(x) <- c("smpl", paste0("loci", seq(1:(ncol(x)-1))))
  ret <- gather(x, key = "loci", value = "allele", 2:ncol(x))
  ret <- ret %>%
    dplyr::mutate(loci = factor(loci,
                                levels = c(paste0("loci", seq(1:(ncol(x)-1)))), # perserve order
                                ordered = T)
    )
  return(ret)
})


plothap <- function(hapdf, region, name){
  plotObj <- ggplot() +
    geom_tile(data=hapdf, aes(x=factor(loci), y=smpl, fill = factor(allele))) +
    scale_fill_manual("Allele", values = c("#d73027", "#fee090", "#66bd63", "#4575b4")) +
    ggtitle(paste0("Gene ", name, " for region ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_blank()
          )
}
crt$happlot <-  pmap(crt[,c("region", "name", "hapdf")], plothap)


