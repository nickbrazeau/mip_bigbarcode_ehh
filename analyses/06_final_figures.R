library(tidyverse)
library(gridExtra)
source("~/Documents/GitHub/mip_bigbarcode_ehh/R/01-EHH_tools.R")
source("~/Documents/GitHub/mip_bigbarcode_ehh/R/02-spiderplot.R")


#------------------------------------------------------------------------------------------------------------------------------------------------------------
# EHH and Bifurcation plots
#------------------------------------------------------------------------------------------------------------------------------------------------------------
load("~/Documents/GitHub/mip_bigbarcode_ehh/data/derived/03-analyze_crt_dhps_drc.rda")

#--------------------------------------------------------------------------------------
# Remake haplohh to have bifurcation plots that can use pmap with bp instead of cM
#--------------------------------------------------------------------------------------
bifurmap <- crtdhps_sub %>%
  dplyr::select(-c("ehh", "pmap", "thap", "haplohh", "scanhh")) %>%
  dplyr::filter(mut_name %in% c("K76T", "G437A", "K540E", "A581G"))

bifurmap$pmap <- purrr::map(bifurmap$vcfRobj_new, getpmap.polarize.nucbp)
bifurmap$thap <- purrr::map(bifurmap$vcfRobj_new, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/bifurmap/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(bifurmap)){

  file = paste0(outdir, bifurmap$name[i], "-", bifurmap$region[i])

  write.table(x = bifurmap$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =F, col.names = F)

  write.table(x = bifurmap$pmap[[i]], file = paste0(file, ".inp"),
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
bifurmap <- rehhfile %>%
  dplyr::select(c("name", "region", "haplohh")) %>%
  dplyr::left_join(x = bifurmap, y=., by = c("name", "region"))



#--------------------------------------------------------------------------------------
# Make Spider Plots
#--------------------------------------------------------------------------------------
#..................................................
#                     CRT
#..................................................
#.................................
# loci crt K76T DRC 1.2
#.................................
bifur_crt_drc1.2 <- bifurmap %>%
  dplyr::filter(region == "DRC1.2" & mut_name == "K76T")



bifur_crt_drc1.2_spiderplot <- spiderplot(hh = bifur_crt_drc1.2$haplohh[[1]],
           focal = 403625,
           pmapobj = bifur_crt_drc1.2$pmap[[1]],
           nucleotides = T,
           left = 5, # marker here is 13
           right = 11, # max here would be 66-13
           max.haps = 2,
           palette = "RdBu",
           reverse = FALSE,
           relabel = NULL) +
  xlab("Position (kbp)") + ggtitle("K76T Haplotype Bifurcation in DRC-East") +
  scale_x_genome(scale = 1000) +
  theme(
       strip.text = element_text(family = "Arial", face = "bold", size = 10),
       axis.text = element_text(family = "Arial", size = 9, angle = 90),
       axis.title.x = element_text(family = "Arial", face = "bold", size = 10),
       plot.title = element_text(family = "Arial", face = "bold", size = 12, hjust = 0.5, vjust=0.5),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank(),
       panel.border = element_rect(colour = "black", fill=NA))



#.................................
# loci crt K76T DRC 2.2
#.................................
bifur_crt_drc2.2 <- bifurmap %>%
  dplyr::filter(region == "DRC2.2" & mut_name == "K76T")


cols <- RColorBrewer::brewer.pal(11, "RdBu")
cols <- cols[ c(8:11, 4:1) ]

bifur_crt_drc2.2_spiderplot <- spiderplot(hh = bifur_crt_drc2.2$haplohh[[1]],
           focal = 403625,
           pmapobj = bifur_crt_drc2.2$pmap[[1]],
           nucleotides = T,
           left = 5, # marker here is 13
           right = 11, # max here would be 66-13
           max.haps = 2,
           palette = "RdBu",
           reverse = FALSE,
           relabel = NULL) +
  xlab("Position (kbp)") + ggtitle("K76T Haplotype Bifurcation in DRC-West") +
  scale_x_genome(scale = 1000) +
  scale_colour_manual(values = cols) +
  theme(
    strip.text = element_text(family = "Arial", face = "bold", size = 10),
    axis.text = element_text(family = "Arial", size = 9, angle = 90),
    axis.title.x = element_text(family = "Arial", face = "bold", size = 10),
    plot.title = element_text(family = "Arial", face = "bold", size = 12, hjust = 0.5, vjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA))


#..................................................
#                     DHPS
#..................................................

#.................................
# loci dhps K540E DRC 1.2
#.................................
bifur_dhps_drc1.2 <- bifurmap %>%
  dplyr::filter(region == "DRC1.2" & mut_name == "K540E")

cols <- RColorBrewer::brewer.pal(11, "RdBu")
cols <- cols[ c(7:10, 5:1) ]


bifur_dhps_drc1.2_spiderplot <- spiderplot(hh = bifur_dhps_drc1.2$haplohh[[1]],
                                           focal = 549993,
                                           pmapobj = bifur_dhps_drc1.2$pmap[[1]],
                                           nucleotides = T,
                                           left = 6, # marker here is 13
                                           right = 6, # max here would be 66-13
                                           max.haps = 2,
                                           palette = "RdBu",
                                           reverse = FALSE,
                                           relabel = NULL) +
  xlab("Position (kbp)") + ggtitle("K540E Haplotype Bifurcation in DRC-East") +
  scale_x_genome(scale = 1000) +
  scale_colour_manual(values = cols) +
  theme(
    strip.text = element_text(family = "Arial", face = "bold", size = 10),
    axis.text = element_text(family = "Arial", size = 9, angle = 90),
    axis.title.x = element_text(family = "Arial", face = "bold", size = 10),
    plot.title = element_text(family = "Arial", face = "bold", size = 12, hjust = 0.5, vjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA))



#.................................
# loci dhps K540E DRC 2.2
#.................................
bifur_dhps_drc2.2 <- bifurmap %>%
  dplyr::filter(region == "DRC2.2" & mut_name == "K540E")

cols <- RColorBrewer::brewer.pal(11, "RdBu")
cols <- cols[ c(7:10, 3:1) ]

bifur_dhps_drc2.2_spiderplot <- spiderplot(hh = bifur_dhps_drc2.2$haplohh[[1]],
           focal = 549993,
           pmapobj = bifur_dhps_drc2.2$pmap[[1]],
           nucleotides = T,
           left = 6, # marker here is 13
           right = 6, # max here would be 66-13
           max.haps = 2,
           palette = "RdBu",
           reverse = FALSE,
           relabel = NULL) +
  xlab("Position (kbp)") + ggtitle("K540E Haplotype Bifurcation in DRC-West") +
  scale_x_genome(scale = 1000) +
  scale_colour_manual(values = cols) +
  theme(
    strip.text = element_text(family = "Arial", face = "bold", size = 10),
    axis.text = element_text(family = "Arial", size = 9, angle = 90),
    axis.title.x = element_text(family = "Arial", face = "bold", size = 10),
    plot.title = element_text(family = "Arial", face = "bold", size = 12, hjust = 0.5, vjust=0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA))







#--------------------------------------------------------------------------------------
# Bring it home for DRC crt
#--------------------------------------------------------------------------------------
drc1.2_K76T_ehhdecay <- crossehhplotdf$ehhplot[[2]] + ggtitle("EHH Decay Plot for K76T in DRC-East")
drc2.2_K76T_ehhdecay <- crossehhplotdf$ehhplot[[8]] + ggtitle("EHH Decay Plot for K76T in DRC-West")



jpeg(file = "results/final_figures/K76T_DRC_EvW_EHH_Bifur.jpg", width = 11, height = 8, units = "in", res = 800)
cowplot::plot_grid(drc1.2_K76T_ehhdecay, drc2.2_K76T_ehhdecay,
                   bifur_crt_drc1.2_spiderplot, bifur_crt_drc2.2_spiderplot,
                   labels = c("A", "B", "C", "D"), ncol = 2)
graphics.off()

#--------------------------------------------------------------------------------------
# Bring it home for DRC dhps
#--------------------------------------------------------------------------------------
drc1.2_K540E_ehhdecay <- crossehhplotdf$ehhplot[[5]] + ggtitle("EHH Decay Plot for K540E in DRC-East")
drc2.2_K540E_ehhdecay <- crossehhplotdf$ehhplot[[11]] + ggtitle("EHH Decay Plot for K540E in DRC-West")



jpeg(file = "results/final_figures/K540E_DRC_EvW_EHH_Bifur.jpg", width = 11, height = 8, units = "in", res = 800)
cowplot::plot_grid(drc1.2_K540E_ehhdecay, drc2.2_K540E_ehhdecay,
                   bifur_dhps_drc1.2_spiderplot, bifur_dhps_drc2.2_spiderplot,
                   labels = c("A", "B", "C", "D"), ncol = 2)
graphics.off()


#--------------------------------------------------------------------------------------
# Show All EHH
#--------------------------------------------------------------------------------------
plts_ehh <- crossehhplotdf %>%
  dplyr::filter(! mutation %in% "I356T") %>%
  dplyr::arrange(mutation) %>%
  .$ehhplot

jpeg(filename = "results/final_figures/All_EHHdecay_plots.jpg", width = 8, height = 11, units = "in", res=800)
do.call("grid.arrange", c(plts_ehh, ncol=2))
graphics.off()



#------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make Haplotype Plots
#------------------------------------------------------------------------------------------------------------------------------------------------------------
haplotypes_sub <- readRDS("~/Documents/GitHub/mip_bigbarcode_ehh/data/derived/05-haplotype_plots_allsmpls.rds")

#.....................
# K76T
#.....................
plts_hap_k76T <- haplotypes_sub %>%
  dplyr::filter(mut_name == "K76T") %>%
  .$happlot

jpeg(filename = "results/final_figures/K76T_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
do.call("grid.arrange", c(plts_hap_k76T, nrow=2))
graphics.off()



#.....................
# G437A
#.....................
plts_hap_G437A <- haplotypes_sub %>%
  dplyr::filter(mut_name == "G437A") %>%
  .$happlot

jpeg(filename = "results/final_figures/G437A_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
do.call("grid.arrange", c(plts_hap_G437A, nrow=2))
graphics.off()




#.....................
# K540E
#.....................
plts_hap_K540E <- haplotypes_sub %>%
  dplyr::filter(mut_name == "K540E") %>%
  .$happlot

jpeg(filename = "results/final_figures/K540E_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
do.call("grid.arrange", c(plts_hap_K540E, nrow=2))
graphics.off()




#.....................
# A581G
#.....................
plts_hap_A581G <- haplotypes_sub %>%
  dplyr::filter(mut_name == "A581G") %>%
  .$happlot

jpeg(filename = "results/final_figures/A581G_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
do.call("grid.arrange", c(plts_hap_A581G, nrow=2))
graphics.off()






















