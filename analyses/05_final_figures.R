library(tidyverse)
library(gridExtra)
source("~/Documents/GitHub/mip_bigbarcode_ehh/R/02-spiderplot.R")
load("~/Documents/GitHub/mip_bigbarcode_ehh/data/derived/03-analyze_crt_dhps_drc.rda")


#------------------------------------------------------------------------------------------------------------------------------------------------------------
# EHH and Haplotype Maps
#------------------------------------------------------------------------------------------------------------------------------------------------------------
#.................................
# crt loci 76T
#.................................
plts_ehh <- crossehhplotdf %>%
  dplyr::filter(mutation == "K76T") %>%
  .$ehhplot


plts_hap <- crossehhplotdf %>%
  dplyr::filter(mutation == "K76T") %>%
  .$happlot



#................
# bring it home
#...............

jpeg(filename = "results/final_figures/K76T_haplotype_EHHdecay.jpg", width = 8, height = 11, units = "in", res=800)
plts <- c(plts_ehh, plts_hap)
do.call("grid.arrange", c(plts,ncol=2))
graphics.off()





#.................................
# dhps loci K540E
#................................
plts_ehh <- crossehhplotdf %>%
  dplyr::filter(mutation == "K540E") %>%
  .$ehhplot


plts_hap <- crossehhplotdf %>%
  dplyr::filter(mutation == "K540E") %>%
  .$happlot

#................
# bring it home
#...............

jpeg(filename = "results/final_figures/K540E_haplotype_EHHdecay.jpg", width = 8, height = 11, units = "in", res=800)
plts <- c(plts_ehh, plts_hap)
do.call("grid.arrange", c(plts,ncol=2))
graphics.off()



#.................................
# dhps loci A581G
#................................
plts_ehh <- crossehhplotdf %>%
  dplyr::filter(mutation == "A581G") %>%
  .$ehhplot


plts_hap <- crossehhplotdf %>%
  dplyr::filter(mutation == "A581G") %>%
  .$happlot

#................
# bring it home
#...............

jpeg(filename = "results/final_figures/A581G_haplotype_EHHdecay.jpg", width = 8, height = 11, units = "in", res=800)
plts <- c(plts_ehh, plts_hap)
do.call("grid.arrange", c(plts,ncol=2))
graphics.off()



#.................................
# bring it home
#.................................
plts <- c(plts_ehh, plts_hap)
do.call("grid.arrange", c(plts, ncol=3, nrow=4))






#...........................................................................
# make bifurcation plot for crt loci 76T
#...........................................................................
bifurmap <- crtdhps_sub %>%
  dplyr::filter(mut_name == "K76T") %>%
  dplyr::filter(! region %in% c("Uganda", "Zambia", "Ghana")) %>% # lack variation, can't make a bifurcation plot
  dplyr::select(c("haplohh", "pmap", "cM_Pos")) %>%
  dplyr::rename(hh = haplohh,
                pmapobj = pmap,
                focal = cM_Pos) %>%
  dplyr::mutate(nucleotides = T,
                nucleotidetrim_left = 15.5,
                nucleotidetrim_right = 16.0,
                left = 6, # marker here is 13
                right = 6, # max here would be 66-13
                max.haps = 2,
                palette = "RdBu",
                reverse = FALSE,
                relabel = NULL)

crtK76T_spdrplots <- pmap(bifurmap, spiderplot)

# check change fro DRC1.2
cols <- RColorBrewer::brewer.pal(11, "RdBu")
cols <- cols[ c(8:11, 4:1) ]

crtK76T_spdrplots[[1]] <- crtK76T_spdrplots[[1]] +
  scale_colour_manual(values = cols) +
  ggtitle("Bifurcation Plot for K76T in DRC1.2") +
  theme(
    plot.title = element_text(family = "Arial", hjust = 0.5, vjust = 0.5, size = 12, face = "bold")
  )

crtK76T_spdrplots[[2]] <- crtK76T_spdrplots[[2]] +
  ggtitle("Bifurcation Plot for K76T in DRC2.2") +
  theme(
    plot.title = element_text(family = "Arial", hjust = 0.5, vjust = 0.5, size = 12, face = "bold")
  )




crtK76T_spdrplots[[1]]
crtK76T_spdrplots[[2]]




#...........................................................................
# make bifurcation plot for dhps loci K540E
#...........................................................................
bifurmap <- crtdhps_sub %>%
  dplyr::filter(mut_name == "K540E") %>%
  dplyr::filter(! region %in% c("Uganda", "Ghana")) %>% # lack variation, can't make a bifurcation plot
  dplyr::select(c("haplohh", "pmap", "cM_Pos")) %>%
  dplyr::rename(hh = haplohh,
                pmapobj = pmap,
                focal = cM_Pos) %>%
  dplyr::mutate(nucleotides = T,
                nucleotidetrim_left = 30,
                nucleotidetrim_right = 33,
                left = 6, # marker here is 13
                right = 6, # max here would be 66-13
                max.haps = 2,
                palette = "RdBu",
                reverse = FALSE,
                relabel = NULL)

dhpsK540E_spdrplots <- pmap(bifurmap, spiderplot)



dhpsK540E_spdrplots[[1]] <- dhpsK540E_spdrplots[[1]] +
#  scale_colour_manual(values = cols) +
  ggtitle("Bifurcation Plot for K540E in DRC1.2") +
  theme(
    plot.title = element_text(family = "Arial", hjust = 0.5, vjust = 0.5, size = 12, face = "bold")
  )
dhpsK540E_spdrplots[[1]]

# check change for drc 2.2
cols <- RColorBrewer::brewer.pal(11, "RdBu")
cols <- cols[ c(10:11, 4:1) ]
dhpsK540E_spdrplots[[2]] <- dhpsK540E_spdrplots[[2]] +
  scale_colour_manual(values = cols) +
  ggtitle("Bifurcation Plot for K540E in DRC2.2") +
  theme(
    plot.title = element_text(family = "Arial", hjust = 0.5, vjust = 0.5, size = 12, face = "bold")
  )
dhpsK540E_spdrplots[[2]]

# check change for Zambia
cols <- RColorBrewer::brewer.pal(11, "RdBu")
cols <- cols[ c(10:11, 3:1) ]
dhpsK540E_spdrplots[[3]] <- dhpsK540E_spdrplots[[3]] +
  scale_colour_manual(values = cols) +
  ggtitle("Bifurcation Plot for K540E in Zambia") +
  theme(
    plot.title = element_text(family = "Arial", hjust = 0.5, vjust = 0.5, size = 12, face = "bold")
  )


dhpsK540E_spdrplots[[3]]




