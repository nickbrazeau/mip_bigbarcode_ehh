#-------------------------------------------------------------------------------------------------------
# Purpose of this script is to analyze the DR panel and BB panel in DRC samples
# loci for differences in selection
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
drc <- mc %>%
  magrittr::set_colnames(tolower(colnames(.))) %>%
  dplyr::filter(country == "DRC") %>%
  dplyr::rename(region = country) %>%
  group_by(region) %>%
  nest()
drc$smpls <- purrr::map(drc$data, "id")

drglong <- do.call("rbind", lapply(1:nrow(drugres), function(x){
  return(drc)
})) %>%
  dplyr::arrange(region) # make a copy of the region data for each drug res site (not memory efficient but fits within long tidy framework)

drugregions <- as_tibble(cbind.data.frame(drglong, drugres))


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
  filter(nvar > 20)

#.................................
# make map and haplotype files
#.................................
drugregions_sub$pmap <- purrr::map(drugregions_sub$vcfRobj, getpmap.polarize)
drugregions_sub$thap <- purrr::map(drugregions_sub$vcfRobj, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/ehhdrcpop/"
if(!dir.exists(outdir)){
  dir.create(outdir, recursive = T)
}

for(i in 1:nrow(drugregions_sub)){

  file = paste0(outdir, drugregions_sub$name[i], "-", drugregions_sub$region[i])

  write.table(x = drugregions_sub$thap[[i]], file = paste0(file, ".thap"),
              sep = " ", quote = F, row.names =F, col.names = F)

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
                                                                                               min_perc_geno.hap = 0,
                                                                                               min_perc_geno.snp = 1, # if entire SNP is missing information, it is due to the fact that we weren't able to polarize it (otherwise would have been dropped upstream); need to drop those as that information is useless
                                                                                               min_maf = 0)) # !! TWEAKING NEEDED


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



#..............................................
# Subset to Regions of Interest Identified by Bob's Final Drug Res Sites
#.............................................

# do a full join off of Bob's final Drug Res sites



#--------------------------------------------------------------------
# Manually look at DHPS
#--------------------------------------------------------------------
dhps <- drugregions_sub %>%
  dplyr::filter(name == "dhps")

dhps$ehh <- map(dhps$haplohh, ~rehh::calc_ehh(
  haplohh = .x,
  mrk = which(dhps$haplohh[[1]]@snp.name == "Pf3D7_08_v3_549993"),
  limhaplo = 2,
  limehh = 0.05,
  plotehh = F)
)



dhps$pos <- map(dhps$haplohh, "position")
dhps$ehh <- map(dhps$ehh, "ehh")
dhps$ehh <- map(dhps$ehh, function(x){ return(as.data.frame(t(x))) } )

dhps$ehhdf <- map2(dhps$pos, dhps$ehh, ~data.frame(pos = .x,
                                                   ancestral = .y[,1],
                                                   derived = .y[,2]))
dhps$ehhdf <- map(dhps$ehhdf, function(x){
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
    ggtitle(paste0("This is marker: ", marker, " for region ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

dhps$marker <- "Pf3D7_08_v3_549993" # TEMP
dhps$ehhplot <-  pmap(dhps[,c("region", "marker", "ehhdf")], plotehh)


#.................................
# make spider plots
#.................................
focal = dhps$haplohh[[1]]@position[ which(dhps$haplohh[[1]]@snp.name == "Pf3D7_08_v3_549993") ] # 11 is where

dhps$bifurcplot <- map2(dhps$haplohh, dhps$pmap, ~spiderplot(hh = .x,
                                        pmapobj = .y,
                                        focal = focal,
                                        nucleotides = T,
                                        nucleotidetrim = 5,
                                        left = 10,
                                        right = 10,
                                        max.haps = 2,
                                        palette = "RdBu",
                                        reverse = FALSE,
                                        relabel = NULL))



#.................................
# make haplotype plots
#.................................
dhps$hapdf <- map(dhps$thap, function(x){
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


plothap <- function(hapdf, region, name){
  plotObj <- ggplot() +
    geom_tile(data=hapdf, aes(x=loci, y=smpl, fill = factor(allele))) +
    scale_fill_manual("Allele", values = c("#d73027", "#fee090", "#66bd63", "#4575b4")) +
    ggtitle(paste0("Gene ", name, " for region ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
          axis.text.y = element_blank())
}
dhps$happlot <-  pmap(dhps[,c("region", "name", "hapdf")], plothap)





