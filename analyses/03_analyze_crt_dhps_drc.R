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
# subset to DRC and crt/dhps
#.................................
drugregions <- drugregions %>%
  dplyr::filter(grepl("DRC", region)) %>%
  dplyr::filter(name %in% c("crt", "dhps")) %>%
  dplyr::arrange(., name)


#.................................
# subset vcf to crt and dhps
#.................................

subsetlocismpls <- function(smpls, seqname, start, end, geneid){

  chromposbed = tibble(seqname = seqname, start = start, end = end, geneid = geneid)
  ret <- vcfRmanip::vcfR2SubsetChromPos(vcfRobject = mipbb_dr_panel_vcfR, # global
                                        chromposbed = chromposbed) # subset loci
  ret <- ret[, colnames(ret@gt) %in% c("FORMAT", smpls)] # subset samples
  return(ret)
}

drugregions$vcfRobj <- purrr::pmap(drugregions[,c("smpls", "seqname", "start", "end", "geneid")],
                                   subsetlocismpls)


#.................................
# check data n
#.................................
drugregions$nsmpls <- unlist(purrr::map(drugregions$vcfRobj, function(x){ return(ncol(x@gt)-1) }))
# looks reasonably spread out


#.................................
# make map and haplotype files for crt
#.................................
drugregions$pmap <- purrr::map(drugregions$vcfRobj, getpmap.polarize)
drugregions$thap <- purrr::map(drugregions$vcfRobj, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/final_dhps_crt_select_drc/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(drugregions)){

  file = paste0(outdir, drugregions$name[i], "-", drugregions$region[i])

  write.table(x = drugregions$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =F, col.names = F)

  write.table(x = drugregions$pmap[[i]], file = paste0(file, ".inp"),
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

rehhfile$haplohh <- purrr::map2(rehhfile$thap_files, rehhfile$pmap_files,
                                ~rehh::data2haplohh(hap_file=.x, map_file=.y,
                                                    haplotype.in.columns = TRUE,
                                                    recode.allele = T,
                                                    min_perc_geno.hap = 0,
                                                    min_perc_geno.snp = 1, # if entire SNP is missing information, it is due to the fact that we weren't able to polarize it (otherwise would have been dropped upstream); need to drop those as that information is useless
                                                    min_maf = 0))



#.................................
# pull it all together
#.................................
drugregions <- rehhfile %>%
  dplyr::select(c("name", "region", "haplohh")) %>%
  dplyr::left_join(x = drugregions, y=., by = c("name", "region"))


#.................................
# call ihs calcs
#.................................
drugregions$scanhh <- purrr::map(drugregions$haplohh, rehh::scan_hh,
                         limhaplo = 2,
                         limehh = 0.05)

#.................................
# manually set markers
#.................................

# chr7	403621 crt-K76T	mark
K76Tmark <- which(drugregions$haplohh[[1]]@snp.name == "Pf3D7_07_v3_403623")

# chr8 549993 dhps-K540E
K540Emark <- which(drugregions$haplohh[[3]]@snp.name == "Pf3D7_08_v3_549993")

# chr8 550116 dhps-A581G
A581Gmark <- which(drugregions$haplohh[[3]]@snp.name == "Pf3D7_08_v3_550117") # TODO why is this one based off


mrk <- tibble::tibble(name   = c("crt",     "dhps",    "dhps"),
                      marker = c(K76Tmark, K540Emark, A581Gmark),
                      markername = c("K76T", "K540E", "A581G"))

drugregions <- left_join(mrk, drugregions, by = c("name"))

#.................................
# call ehh based on the marker that we picked
#.................................
drugregions$ehh <- map2(drugregions$haplohh, drugregions$marker, ~rehh::calc_ehh(
  haplohh = .x,
  mrk = .y,
  limhaplo = 2,
  limehh = 0.05,
  plotehh = F)
)

#.................................
# make plots for ehh continous by region
#.................................
crossehhplotdf <- tibble(geneid = drugregions$geneid,
                         marker = unlist( drugregions$marker ),
                         markername = drugregions$markername,
                         region = drugregions$region,
                         name = drugregions$name )



crossehhplotdf$pos <- map(drugregions$haplohh, "position")
crossehhplotdf$ehh <- map(drugregions$ehh, "ehh")
crossehhplotdf$ehh <- map(crossehhplotdf$ehh, function(x){ return(as.data.frame(t(x))) } )

crossehhplotdf$ehhdf <- map2(crossehhplotdf$pos, crossehhplotdf$ehh, ~data.frame(pos = .x,
                                                                                 ancestral = .y[,1],
                                                                                 derived = .y[,2]))
crossehhplotdf$ehhdf <- map(crossehhplotdf$ehhdf, function(x){
  ret <- gather(x, key = "allele", value = "ehh", 2:3)
  return(ret)
})

#.................................
# make plots for ehh continous
#.................................
plotehh <- function( region, markername, geneid, ehhdf){
  plotObj <- ggplot() +
    geom_line(data=ehhdf, aes(x=pos, y=ehh, colour = factor(allele), group = factor(allele))) +
    scale_color_manual("Allele Status", values = c("#4575b4", "#d73027")) +
    ggtitle(paste0("This is marker: ", markername, " on gene ", geneid, " for region ", region)) +
    labs(x = "position (cM)", y = "EHH Statistic") +
    scale_color_manual("Allele", labels = c("Ancestral", "Derived"), values = c("#2166ac", "#b2182b")) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

crossehhplotdf$ehhplot <-  pmap(crossehhplotdf[,c("region", "markername", "geneid", "ehhdf")], plotehh)

#.................................
# make haplotype plots
#.................................
crossehhplotdf$hapdf <- map(drugregions$thap, function(x){
  x <- t(x)
  x <- as.data.frame(cbind(paste0("smpl", seq(1:nrow(x))), x))
  colnames(x) <- c("smpl", paste0("loci", seq(1:(ncol(x)-1))))
  ret <- gather(x, key = "loci", value = "allele", 2:ncol(x))
  ret <- ret %>%
    dplyr::mutate(loci = factor(loci,
                                levels =  paste0("loci", seq(1:(ncol(x)-1))),
                                ordered = T))

  return(ret)
})


plothap <- function(hapdf, region, marker, markername, geneid){
  plotObj <- ggplot() +
    geom_tile(data=hapdf, aes(x=loci, y=smpl, fill = factor(allele))) +
    scale_fill_manual("Allele", values = c("#d73027", "#fee090", "#66bd63", "#4575b4")) +
    ggtitle(paste0("This is loci: ",  marker, " via ", markername, " on gene ", geneid, " for region ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 5),
          axis.text.y = element_blank()
    )
}
crossehhplotdf$happlot <-  pmap(crossehhplotdf[,c("hapdf", "region", "marker", "markername", "geneid")], plothap)


crossehhplotdf$ehhplot[[1]]


#.................................
# make bifurcation plots
#.................................
bifurmap <- drugregions %>%
  dplyr::select(c("haplohh", "pmap", "marker")) %>%
  dplyr::rename(hh = haplohh,
                pmapobj = pmap,
                focal = marker) %>%
  dplyr::mutate(nucleotides = T,
                nucleotidetrim = 5,
                left = 8,
                right = 8,
                max.haps = 2,
                palette = "RdBu",
                reverse = FALSE,
                relabel = NULL)


crossehhplotdf$spdrplot <- pmap(bifurmap, spiderplot)

save(crossehhplotdf, drugregions,
        file = "results/03-analyze_crt_dhps_drc.rda")

