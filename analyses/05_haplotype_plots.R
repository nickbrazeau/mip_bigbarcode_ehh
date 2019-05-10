#-------------------------------------------------------------------------------------------------------
# Purpose of this script is to make haplotype plots for
# crt loci 76T and the dhps 437 540 and 581 loci
# with all samples, not just those with no missing data
#-------------------------------------------------------------------------------------------------------
source("R/01-EHH_tools.R")
source("analyses/00_metadata.R")
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
library(rehh)
load("data/derived/02-drugres_obj_rehh_subpopulationlevel.rda")
mipDRpanel <- readRDS("data/raw_snps_filtered/dr_monoclonal.rds")
#-------------------------------------------------------------------------------------------------------
# subset to data to just what we need to start. will make haplotype plots from all data now
##------------------------------------------------------------------------------------------------------
haplotypes_sub <- drugregions_sub %>%
  dplyr::select(-c("vcfRobj_new", "pmap", "thap", "haplohh", "scanhh", "nvar", "marker")) %>% # remove marker so I can pair it with snp name later
  dplyr::filter(name %in% c("dhps", "crt"))

#-------------------------------------------------------------------------------------------------------
# make map and haplotype files from all smpls
#------------------------------------------------------------------------------------------------------
haplotypes_sub$pmap <- purrr::map(haplotypes_sub$vcfRobj, getpmap.polarize)
haplotypes_sub$thap <- purrr::map(haplotypes_sub$vcfRobj, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/haplotypeplots/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(haplotypes_sub)){

  file = paste0(outdir, haplotypes_sub$name[i], "-", haplotypes_sub$region[i])

  write.table(x = haplotypes_sub$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =F, col.names = F)

  write.table(x = haplotypes_sub$pmap[[i]], file = paste0(file, ".inp"),
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
                                                                                               min_perc_geno.snp = 0,
                                                                                               min_maf = 0)) # Assuming Bob's filters upstream good enough

#.................................
# pull it all together
#.................................
haplotypes_sub <- rehhfile %>%
  dplyr::select(c("name", "region", "haplohh")) %>%
  dplyr::left_join(x = haplotypes_sub, y=., by = c("name", "region"))



#-------------------------------------------------------------------------------------------------------
# Find putative drug resistance markers that we will make haplotype plots for
#------------------------------------------------------------------------------------------------------
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
  dplyr::mutate(snppresent = snpname %in% haplotypes_sub$haplohh[[ which(haplotypes_sub$name == "crt")[1] ]]@snp.name)

markpos <- sapply(crtsites$snpname, function(x){
  if(x %in% haplotypes_sub$haplohh[[ which(haplotypes_sub$name == "crt")[1] ]]@snp.name ){
    marker <- which(haplotypes_sub$haplohh[[ which(haplotypes_sub$name == "crt")[1] ]]@snp.name == x)
    cM_Pos <- haplotypes_sub$haplohh[[ which(haplotypes_sub$name == "crt")[1] ]]@position[marker]
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
  dplyr::mutate(snppresent = snpname %in% haplotypes_sub$haplohh[[ which(haplotypes_sub$name == "dhps")[1] ]]@snp.name)

markpos <- sapply(dhpssites$snpname, function(x){
  if(x %in% haplotypes_sub$haplohh[[ which(haplotypes_sub$name == "dhps")[1] ]]@snp.name ){
    marker <- which(haplotypes_sub$haplohh[[ which(haplotypes_sub$name == "dhps")[1] ]]@snp.name == x)
    cM_Pos <- haplotypes_sub$haplohh[[ which(haplotypes_sub$name == "dhps")[1] ]]@position[marker]
    return(cbind.data.frame(marker = marker, cM_Pos = cM_Pos))
  } else {
    return(NA)
  }
}) %>% do.call("rbind.data.frame", .)

dhpssites$marker <- markpos$marker
dhpssites$cM_Pos <- markpos$cM_Pos

# confirm missing data for dhps from multiallelic sites
mipDRpanel$loci %>%
  dplyr::filter(CHROM %in% dhpssites$chrom[is.na(dhpssites$marker)]) %>%
  dplyr::filter(POS %in% dhpssites$pos[is.na(dhpssites$marker)])

mrk <- rbind(crtsites, dhpssites) %>%
  dplyr::filter(snppresent) %>%
  dplyr::select(c("gene", "mut_name", "marker", "cM_Pos")) %>%
  dplyr::rename(name = gene)

#......................
# Add makers to df
#......................
haplotypes_sub <- dplyr::left_join(x=haplotypes_sub, y=mrk, by = "name")

#-------------------------------------------------------------------------------------------------------
# Make Haplotype Plots
#------------------------------------------------------------------------------------------------------
# subset to putative loci

haplotypes_sub <- haplotypes_sub %>%
  dplyr::filter(mut_name %in% c("K76T", "G437A", "K540E", "A581G"))



haplotypes_sub$hapdf <- pmap(haplotypes_sub[, c("thap", "pmap", "haplohh", "marker")],
                             function(thap, pmap, haplohh, marker){
                               # have set this up so that we can facet by the original snp we identified
                               # want to be able to look at haplotypes based on ancestral versus derived site

                               thap <- t(thap)
                               thap <- as.data.frame(cbind(paste0("smpl", seq(1:nrow(thap))), thap))
                               colnames(thap) <- c("smpl", pmap$snpname)

                               focalsnpname <-haplohh@snp.name[marker]

                               allelepol <- thap %>%
                                 dplyr::select(c("smpl", focalsnpname)) %>%
                                 dplyr::mutate(snpname = focalsnpname,
                                               focalsnpused = focalsnpname) %>%
                                 dplyr::rename(allele = focalsnpname) %>%
                                 dplyr::left_join(x=., y=pmap, by = "snpname") %>%
                                 dplyr::mutate(
                                   allelepol = ifelse(allele == anc, "A", ifelse(allele == der, "D", NA))
                                 ) %>%
                                 dplyr::select(c("smpl", "allelepol", "focalsnpused"))

                               ret <- thap %>%
                                 dplyr::left_join(x=., y=allelepol, by = "smpl") %>%
                                 tidyr::gather(., key = "snpname", value = "allele", 2:(ncol(.)-2)) %>%
                                 dplyr::rename(loci = snpname)

                               return(ret)
                             })


plothap <- function(hapdf, region, mut_name){
  asterickmarks <- hapdf %>%
    dplyr::select(c("smpl", "allelepol", "focalsnpused")) %>%
    dplyr::filter(!is.na(allelepol)) %>%  # these sites are immediately dropped by rehh, so exclude
    dplyr::group_by(allelepol) %>%
    dplyr::summarise(
      height = mean(as.numeric(factor(smpl))),
      loci = unique(focalsnpused)
    ) %>%
    dplyr::mutate(ast = "*",
                  height = ifelse(height == 1, 0.5, height))

  plotObj <- hapdf %>%
    dplyr::mutate(allele = forcats::fct_explicit_na(allele)) %>%
    dplyr::filter(!is.na(allelepol)) %>%  # these sites are immediately dropped by rehh, so exclude
    ggplot() +
    geom_tile(aes(x=loci, y=smpl, fill = allele)) +
    geom_text(data = asterickmarks, aes(x=loci, y=height, label = ast), size = 8, color = "#000000") +
    facet_grid(allelepol ~ ., scale = "free_y", space = "free_y") +
    scale_fill_manual("Allele", values = c("#d73027", "#fee090", "#66bd63", "#4575b4", "#f0f0f0")) +
    ggtitle(paste0("Haplotype Plot for ", mut_name, " in ", region)) +
    ylab("Samples") + xlab("Loci") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.title = element_text(face = "bold", hjust = 0.5, size = 8.5),
          # axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 3),
          axis.text.x = element_blank(),
          axis.text.y = element_blank()
    )
}

haplotypes_sub$happlot <-  pmap(haplotypes_sub[,c("hapdf", "region", "mut_name")], plothap)




#..........................
# save it out
#..........................

saveRDS(object = haplotypes_sub, file = "data/derived/05-haplotype_plots_allsmpls.rds")


