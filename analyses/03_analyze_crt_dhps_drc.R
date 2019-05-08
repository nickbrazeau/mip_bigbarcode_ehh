#-------------------------------------------------------------------------------------------------------
# After discussing with JJJ, JAB & SRM
# have decided to subset to DHPS and CRT and will considerd the
# DRC as a split population using the K-means of 2 cluster
#
# We are going to subset to just look at the
# crt loci 76T and the dhps 540 and 581 loci
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
library(rehh)
set.seed(42)


load("data/derived/02-drugres_obj_rehh_subpopulationlevel.rda")
mipDRpanel <- readRDS("data/raw_snps_filtered/dr_monoclonal.rds")


#.................................
# subset to data to crt/dhps
#.................................
crtdhps_sub <- drugregions_sub %>%
  dplyr::filter(name %in% c("crt", "dhps")) %>%
  dplyr::select(-c("marker", "ehh"))


#.................................
# find putiative markers markers
#.................................
# This is the supplementary Table that Bob V put together
# Details the prevalence of putative drug resistance loci

chromliftover <- data.frame(CHROM = rplasmodium::chromnames()) %>%
  dplyr::filter(! CHROM %in% c("Pf_M76611", "Pf3D7_API_v3")) %>%
  dplyr::mutate(
    chrom = paste0("chr", 1:14))


crtdhps_drugsites  <-read_csv("data/Supplemental Table 2_ Drug resistance MIP panel prevalences - Sheet1.csv") %>%
  dplyr::filter(gene %in% c("crt", "dhps")) %>%
  dplyr::mutate(
    overall_perc = stringr::str_split_fixed(overall, "\\(", n=2)[,2],
    prev = stringr::str_extract(overall_perc,  "\\d+\\.*\\d*") #https://stackoverflow.com/questions/19252663/extracting-decimal-numbers-from-a-string
  ) %>%
  dplyr::filter(prev > 10) %>% # subset to putative sites with reasonable support
  dplyr::left_join(x=., y=chromliftover, by = "chrom") %>%
  dplyr::mutate(snpname = paste0(CHROM, "_", pos))

# find markers for crt
crtsites <- crtdhps_drugsites %>%
  dplyr::filter(gene == "crt") %>%
  dplyr::mutate(snppresent = snpname %in% crtdhps_sub$haplohh[[ which(crtdhps_sub$name == "crt")[1] ]]@snp.name)

markpos <- sapply(crtsites$snpname, function(x){
  if(x %in% crtdhps_sub$haplohh[[ which(crtdhps_sub$name == "crt")[1] ]]@snp.name ){
    marker <- which(crtdhps_sub$haplohh[[ which(crtdhps_sub$name == "crt")[1] ]]@snp.name == x)
    cM_Pos <- crtdhps_sub$haplohh[[ which(crtdhps_sub$name == "crt")[1] ]]@position[marker]
    return(cbind.data.frame(marker = marker, cM_Pos = cM_Pos))
  } else {
    return(NA)
  }
}) %>% do.call("rbind.data.frame", .)

crtsites$marker <- markpos$marker
crtsites$cM_Pos <- markpos$cM_Pos


# confirm missing data for crt from multiallelic sites
mipDRpanel$loci %>%
  dplyr::filter(CHROM %in% crtsites$chrom[is.na(crtsites$marker)]) %>%
  dplyr::filter(POS %in% crtsites$pos[is.na(crtsites$marker)])

# find markers for dhps
dhpssites <- crtdhps_drugsites %>%
  dplyr::filter(gene == "dhps") %>%
  dplyr::mutate(snppresent = snpname %in% crtdhps_sub$haplohh[[ which(crtdhps_sub$name == "dhps")[1] ]]@snp.name)

markpos <- sapply(dhpssites$snpname, function(x){
  if(x %in% crtdhps_sub$haplohh[[ which(crtdhps_sub$name == "dhps")[1] ]]@snp.name ){
    marker <- which(crtdhps_sub$haplohh[[ which(crtdhps_sub$name == "dhps")[1] ]]@snp.name == x)
    cM_Pos <- crtdhps_sub$haplohh[[ which(crtdhps_sub$name == "dhps")[1] ]]@position[marker]
    return(cbind.data.frame(marker = marker, cM_Pos = cM_Pos))
  } else {
    return(NA)
  }
}) %>% do.call("rbind.data.frame", .)

dhpssites$marker <- markpos$marker
dhpssites$cM_Pos <- markpos$cM_Pos



mrk <- rbind(crtsites, dhpssites) %>%
  dplyr::filter(snppresent) %>%
  dplyr::select(c("gene", "mut_name", "marker", "cM_Pos")) %>%
  dplyr::rename(name = gene)


#......................
# Add makers to df
#......................
crtdhps_sub <- dplyr::left_join(x=crtdhps_sub, y=mrk, by = "name")

#.................................
# call ehh based on the marker that we picked
#.................................
crtdhps_sub$ehh <- map2(crtdhps_sub$haplohh, crtdhps_sub$marker, ~rehh::calc_ehh(
  haplohh = .x,
  mrk = .y,
  limhaplo = 2,
  limehh = 0.05,
  discard_integration_at_border = T,
  plotehh = F)
)

#.................................
# make plots for ehh continous by region
#.................................
crossehhplotdf <- tibble(geneid = crtdhps_sub$geneid,
                         marker = unlist( crtdhps_sub$marker ),
                         mutation = crtdhps_sub$mut_name,
                         region = crtdhps_sub$region,
                         name = crtdhps_sub$name )



crossehhplotdf$pos <- map(crtdhps_sub$haplohh, "position")
crossehhplotdf$ehh <- map(crtdhps_sub$ehh, "ehh")
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
plotehh <- function( region, mutation, name, ehhdf){
  plotObj <- ggplot() +
    geom_line(data=ehhdf, aes(x=pos, y=ehh, colour = factor(allele), group = factor(allele))) +
    scale_color_manual("Allele Status", values = c("#4575b4", "#d73027")) +
    ggtitle(paste0("This is mutation: ", mutation, " on gene ", name, " for region ", region)) +
    labs(x = "position (cM)", y = "EHH Statistic") +
    scale_color_manual("Allele", labels = c("Ancestral", "Derived"), values = c("#2166ac", "#b2182b")) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

crossehhplotdf$ehhplot <-  pmap(crossehhplotdf[,c("region", "mutation", "name", "ehhdf")], plotehh)

#.................................
# make haplotype plots
#.................................
crossehhplotdf$hapdf <- map2(crtdhps_sub$thap, crtdhps_sub$pmap, function(hapsdf = .x, posdf = .y){
  hapsdf <- t(hapsdf)
  hapsdf <- as.data.frame(cbind(paste0("smpl", seq(1:nrow(hapsdf))), hapsdf))
  colnames(hapsdf) <- c("smpl", posdf$snpname)

  ret <- gather(hapsdf, key = "loci", value = "allele", 2:ncol(hapsdf))
  ret <- ret %>%
    dplyr::mutate(loci = factor(loci,
                                levels =  unique(loci),
                                ordered = T))

  return(ret)
})


plothap <- function(hapdf, region, name){
  plotObj <- ggplot() +
    geom_tile(data=hapdf, aes(x=loci, y=smpl, fill = factor(allele))) +
    scale_fill_manual("Allele", values = c("#d73027", "#fee090", "#66bd63", "#4575b4")) +
    ggtitle(paste0("Haplotype for Gene:", name, " in ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 5),
          axis.text.y = element_blank()
    )
}
crossehhplotdf$happlot <-  pmap(crossehhplotdf[,c("hapdf", "region", "name")], plothap)


#.................................
# write out
#.................................
save(crossehhplotdf, crtdhps_sub,
        file = "data/derived/03-analyze_crt_dhps_drc.rda")

