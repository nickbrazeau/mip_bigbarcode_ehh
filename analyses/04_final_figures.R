library(tidyverse)
library(gridExtra)
library(rplasmodium)
source("~/Documents/GitHub/mip_bigbarcode_ehh/R/01-EHH_tools.R")
source("~/Documents/GitHub/mip_bigbarcode_ehh/R/02-spiderplot.R")


#------------------------------------------------------------------------------------------------------------------------------------------------------------
# EHH and Bifurcation plots
#------------------------------------------------------------------------------------------------------------------------------------------------------------
load("~/Documents/GitHub/mip_bigbarcode_ehh/data/derived/01-drugres_obj_rehh_subpopulationlevel.rda")

#--------------------------------------------------------------------------------------
# Remake haplohh to have bifurcation plots that can use pmap with bp instead of cM
#--------------------------------------------------------------------------------------
bifurmap <- drugregions_sub %>%
  dplyr::select(-c("ehh", "pmap", "thap", "haplohh", "scanhh")) %>%
  dplyr::filter(mut_name %in% c("K76T", "I356T", "G437A", "K540E", "A581G"))

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



#--------------------------------------------------------------------------------------------------------------------------------------
# Make Spider Plots
#--------------------------------------------------------------------------------------------------------------------------------------
#..................................................
#                     CRT
#..................................................
#.................................
# loci crt K76T DRC 1.2
#.................................
bifur_crt_drc1.2 <- bifurmap %>%
  dplyr::filter(region == "DRC-East" & mut_name == "K76T")



bifur_crt_K76T_drc1.2_spiderplot <- spiderplot(hh = bifur_crt_drc1.2$haplohh[[1]],
           focal = 403625,
           pmapobj = bifur_crt_drc1.2$pmap[[1]],
           nucleotides = T,
           left = 7, # marker here is 13
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
  dplyr::filter(region == "DRC-West" & mut_name == "K76T")


cols <- RColorBrewer::brewer.pal(11, "RdBu")
cols <- cols[ c(7:11, 4:1) ]

bifur_crt_K76T_drc2.2_spiderplot <- spiderplot(hh = bifur_crt_drc2.2$haplohh[[1]],
           focal = 403625,
           pmapobj = bifur_crt_drc2.2$pmap[[1]],
           nucleotides = T,
           left = 7, # marker here is 13
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
  dplyr::filter(region == "DRC-East" & mut_name == "K540E")

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
  dplyr::filter(region == "DRC-West" & mut_name == "K540E")

cols <- RColorBrewer::brewer.pal(11, "RdBu")
cols <- cols[ c(7:10, 4:1) ]

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
# Bring it home for DRC K76T crt
#--------------------------------------------------------------------------------------
drc1.2_K76T_ehhdecay <- crossehhplotdf %>%
  dplyr::filter(mutation == "K76T" & region == "DRC-East") %>%
  .$ehhplot
drc1.2_K76T_ehhdecay <- drc1.2_K76T_ehhdecay[[1]] + ggtitle("EHH Decay Plot for K76T in DRC-East")


drc2.2_K76T_ehhdecay <- crossehhplotdf %>%
  dplyr::filter(mutation == "K76T" & region == "DRC-West") %>%
  .$ehhplot
drc2.2_K76T_ehhdecay <- drc2.2_K76T_ehhdecay[[1]] + ggtitle("EHH Decay Plot for K76T in DRC-West")


jpeg(file = "results/final_figures/K76T_DRC_EvW_EHH_Bifur.jpg", width = 11, height = 8, units = "in", res = 800)
cowplot::plot_grid(drc1.2_K76T_ehhdecay, drc2.2_K76T_ehhdecay,
                   bifur_crt_K76T_drc1.2_spiderplot, bifur_crt_K76T_drc2.2_spiderplot,
                   labels = c("A", "B", "C", "D"), ncol = 2)
graphics.off()

#--------------------------------------------------------------------------------------
# Bring it home for DRC dhps
#--------------------------------------------------------------------------------------
drc1.2_K540E_ehhdecay <-  crossehhplotdf %>%
  dplyr::filter(mutation == "K540E" & region == "DRC-East") %>%
  .$ehhplot
drc1.2_K540E_ehhdecay <- drc1.2_K540E_ehhdecay[[1]] + ggtitle("EHH Decay Plot for K540E in DRC-East")

drc2.2_K540E_ehhdecay <-  crossehhplotdf %>%
  dplyr::filter(mutation == "K540E" & region == "DRC-West") %>%
  .$ehhplot
drc2.2_K540E_ehhdecay <- drc2.2_K540E_ehhdecay[[1]] + ggtitle("EHH Decay Plot for K540E in DRC-West")

jpeg(file = "results/final_figures/K540E_DRC_EvW_EHH_Bifur.jpg", width = 11, height = 8, units = "in", res = 800)
cowplot::plot_grid(drc1.2_K540E_ehhdecay, drc2.2_K540E_ehhdecay,
                   bifur_dhps_drc1.2_spiderplot, bifur_dhps_drc2.2_spiderplot,
                   labels = c("A", "B", "C", "D"), ncol = 2)
graphics.off()


#--------------------------------------------------------------------------------------
# EHH & Haplotype Plots
#--------------------------------------------------------------------------------------
haplotypes_sub <- readRDS("~/Documents/GitHub/mip_bigbarcode_ehh/data/derived/03-haplotype_plots_allsmpls.rds")


#...................................................
# crt
#...................................................
#.....................
# K76T
#.....................
plts_ehh_k76T <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "N75E") %>% #74,75,76 all same basically, can drop to 1
  .$ehhplot

plts_hap_k76T <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "K76T") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/K76T_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_k76T[[1]],     plts_ehh_k76T[[2]],
                   plts_hap_k76T[[1]],     plts_hap_k76T[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()

#.....................
# I356T
#.....................
plts_ehh_I356T <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "I356T") %>% #74,75,76 all same basically, can drop to 1
  .$ehhplot

plts_hap_I356T <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "I356T") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/I356T_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_I356T[[1]],     plts_ehh_I356T[[2]],
                   plts_hap_I356T[[1]],     plts_hap_I356T[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()

#...................................................
# dhps
#...................................................
#.....................
# G437A
#.....................
plts_ehh_G437A <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "G437A") %>% #74,75,76 all same basically, can drop to 1
  .$ehhplot

plts_hap_G437A <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "G437A") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/G437A_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_G437A[[1]],     plts_ehh_G437A[[2]],
                   plts_hap_G437A[[1]],     plts_hap_G437A[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()


#.....................
# K540E
#.....................
plts_ehh_K540E <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "K540E") %>% #74,75,76 all same basically, can drop to 1
  .$ehhplot

plts_hap_K540E <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "K540E") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/K540E_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_K540E[[1]],     plts_ehh_K540E[[2]],
                   plts_hap_K540E[[1]],     plts_hap_K540E[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()



#.....................
# A581G
#.....................
plts_ehh_A581G <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "A581G") %>% #74,75,76 all same basically, can drop to 1
  .$ehhplot

plts_hap_A581G <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "A581G") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/A581G_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_A581G[[1]],     plts_ehh_A581G[[2]],
                   plts_hap_A581G[[1]],     plts_hap_A581G[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()




#...................................................
# mdr1
#...................................................
#.....................
# N86Y
#.....................
plts_ehh_N86Y <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "N86Y") %>% #74,75,76 all same basically, can drop to 1
  .$ehhplot

plts_hap_N86Y <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "N86Y") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/N86Y_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_N86Y[[1]],     plts_ehh_N86Y[[2]],
                   plts_hap_N86Y[[1]],     plts_hap_N86Y[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()


#.....................
# Y184F
#.....................
plts_ehh_Y184F <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "Y184F") %>% #74,75,76 all same basically, can drop to 1
  .$ehhplot

plts_hap_Y184F <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "Y184F") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/Y184F_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_Y184F[[1]],     plts_ehh_Y184F[[2]],
                   plts_hap_Y184F[[1]],     plts_hap_Y184F[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()


#.....................
# D1246Y
#.....................
plts_ehh_D1246Y <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "D1246Y") %>% #74,75,76 all same basically, can drop to 1
  .$ehhplot

plts_hap_D1246Y <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "D1246Y") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/D1246Y_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_D1246Y[[1]],     plts_ehh_D1246Y[[2]],
                   plts_hap_D1246Y[[1]],     plts_hap_D1246Y[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()




#...................................................
# mdr2
#...................................................
#.....................
# F423Y
#.....................
plts_ehh_F423Y <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "F423Y") %>%
  .$ehhplot

plts_hap_F423Y <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "F423Y") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/F423Y_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_F423Y[[1]],     plts_ehh_F423Y[[2]],
                   plts_hap_F423Y[[1]],     plts_hap_F423Y[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()

#.....................
# I492V
#.....................
plts_ehh_I492V <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "I492V") %>%
  .$ehhplot

plts_hap_I492V <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "I492V") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/I492V_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_I492V[[1]],     plts_ehh_I492V[[2]],
                   plts_hap_I492V[[1]],     plts_hap_I492V[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()



#...................................................
# dhfr
#...................................................
#.....................
# N51I
#.....................
plts_ehh_N51I <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "N51I") %>%
  .$ehhplot

plts_hap_N51I <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "N51I") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/N51I_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_N51I[[1]],     plts_ehh_N51I[[2]],
                   plts_hap_N51I[[1]],     plts_hap_N51I[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()


#.....................
# C59R
#.....................
plts_ehh_C59R <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "C59R") %>%
  .$ehhplot

plts_hap_C59R <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "C59R") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/C59R_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_C59R[[1]],     plts_ehh_C59R[[2]],
                   plts_hap_C59R[[1]],     plts_hap_C59R[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()

#.....................
# S108N
#.....................
plts_ehh_S108N <- crossehhplotdf %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mutation == "S108N") %>%
  .$ehhplot

plts_hap_S108N <- haplotypes_sub %>%
  dplyr::filter(region %in% c("DRC-East", "DRC-West")) %>%
  dplyr::filter(mut_name == "S108N") %>%
  .$happlot

jpeg(filename = "results/final_figures/haplotypeplots/S108N_all_monohaplotypes_plots.jpg", width = 11, height = 8, units = "in", res=800)
cowplot::plot_grid(plts_ehh_S108N[[1]],     plts_ehh_S108N[[2]],
                   plts_hap_S108N[[1]],     plts_hap_S108N[[2]],
                   labels = LETTERS[1:4], nrow = 2)
graphics.off()





#------------------------------------------------------------------------------------------------------------------------------------------------------------
# Standardized XP-EHH plot
#------------------------------------------------------------------------------------------------------------------------------------------------------------
putdrugres_xpehh <- readRDS("data/derived/02-putdrugres_xpehh.rds")

jpeg(filename = "results/final_figures/XP-EHHStats.jpg", width = 11, height = 8, units = "in", res=800)
putdrugres_xpehh %>%
  dplyr::filter(!is.na(pprimexpehh_scale)) %>%
  ggplot() +
  geom_pointrange(aes(x = mut_name, y = crudexpehh, ymin = crudexpehh, ymax = crudexpehh,
                      color = factor(statsig), group = factor(statsig)), size=1.5, alpha = 0.8) +
  facet_grid(name ~ ., scales = "free_y") +
  scale_color_manual("Stat. Sig.", values = c("#313695", "#a50026")) +
  ggtitle("Crude XP-EHH for the Derived Allele \n Comparing the Log Ratio of DRC-East/DRC-West") +
  xlab("Mutation Name") + ylab("XP-EHH Estimate") +
  labs(caption = "DRC-East was classified as Pop1 while DRC-West was classified as Pop2") +
  coord_flip() +
  theme_bw() +
  theme(
    plot.title = element_text(family = "Arial", face = "bold", hjust = 0.5, vjust =0.5, size = 14),
    axis.title = element_text(family = "Arial", face = "bold", hjust = 0.5, vjust =0.5, size = 11)
  )

graphics.off()















