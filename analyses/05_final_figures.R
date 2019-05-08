


#.................................
# make bifurcation plots
#.................................
bifurmap <- crtdhps_sub %>%
  dplyr::select(c("haplohh", "pmap", "cM_Pos")) %>%
  dplyr::rename(hh = haplohh,
                pmapobj = pmap,
                focal = cM_Pos) %>%
  dplyr::mutate(nucleotides = T,
                left = 8,
                right = 8,
                max.haps = 2,
                palette = "RdBu",
                reverse = FALSE,
                relabel = NULL)


crossehhplotdf$spdrplot <- pmap(bifurmap, spiderplot)
