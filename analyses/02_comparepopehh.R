
#-------------------------------------------------------------------------------------------------------
# Purpose of this script is to perform EHH controlling for population structure
#-------------------------------------------------------------------------------------------------------

#.................................
# imports
#.................................
source("analyses/00_metadata.R")
source("R/01-EHH_tools.R")
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
library(rehh)

mipbb_dr_panel_vcfR <- vcfR::read.vcfR(file = "data/polarized_mipbi_drugres_bigbarcode_panels.vcf.bgz")

#.................................
# datawrangle
#.................................
#.....................
# set up drug res sites
#.....................
drugreslist <- split(drugres, factor(1:nrow(drugres)))

#.....................
# remove sites that lack polarization
#.....................

nopolar <- getpmap.polarize(mipbb_dr_panel_vcfR) %>%
  dplyr::filter(is.na(anc)) %>%
  dplyr::mutate(
    CHROM = stringr::str_split_fixed(snpname, "_(?!.*_)", n=2)[,1],
    POS = stringr::str_split_fixed(snpname, "_(?!.*_)", n=2)[,2],
    POS = as.numeric(POS)
  ) %>%
  # needs to be in bed format
  dplyr::rename(chr = CHROM) %>%
  dplyr::mutate(
    start = POS,
    end = POS,
    gene_symbol = ".",
    score = ".",
    strand = ".",
    geneid = snpname,
    seqname = chr
  ) %>%
  dplyr::select(c("chr", "start", "end", "gene_symbol", "score", "geneid", "seqname"))

mipbb_dr_panel_vcfR_ancpolar <- vcfRmanip::vcffilter_ChromPos(vcfRobject = mipbb_dr_panel_vcfR,
                                                              chromposbed = nopolar)


#.....................
# Set up Samples by Regions
#.....................
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
# subset and expand vcf
#.................................

subsetlocismpls <- function(smpls, seqname, start, end, geneid){

  chromposbed = tibble(seqname = seqname, start = start, end = end, geneid = geneid)
  ret <- vcfRmanip::vcfR2SubsetChromPos(vcfRobject = mipbb_dr_panel_vcfR_ancpolar, # global
                                        chromposbed = chromposbed) # subset loci
  ret <- ret[, colnames(ret@gt) %in% c("FORMAT", smpls)] # subset samples
  return(ret)
}

drugregions$vcfRobj <- purrr::pmap(drugregions[,c("smpls", "seqname", "start", "end", "geneid")],
                                   subsetlocismpls)


#.................................
# remove bad sites
#.................................
drugregions$nvar <- unlist(purrr::map(drugregions$vcfRobj, function(x){ return(nrow(x@gt)) }))
drugregions$nsmpls <- unlist(purrr::map(drugregions$vcfRobj, function(x){ return(ncol(x@gt)-1) }))
drugregions_sub <- drugregions %>%
  filter(nvar > 20) %>%
  filter(nsmpls > 5)


#.................................
# make map and haplotype files
#.................................
drugregions_sub$pmap <- purrr::map(drugregions_sub$vcfRobj, getpmap.polarize)
drugregions_sub$thap <- purrr::map(drugregions_sub$vcfRobj, vcfRmanip::vcfR2thap)

outdir <- "~/Desktop/temp_ehhwork/ehhcrosspop/"
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
                                                                                         min_perc_geno.snp = 0,
                                                                                         min_maf = 0)) # Assuming Bob's filters upstream good enough


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
                                     limehh = 0.05,
                                     discard_integration_at_border = F)

# looks like NAs are being produced from
# https://github.com/cran/rehh/blob/master/src/hh_utils.c
# line 412 - 426
#   if (discard_integration_at_border && ((y_axis[0] > threshold) || (y_axis[n - 1] > threshold))) {  // If the EHH or EHHS is larger than the minimum value at either end of the chromosome, ...
#   return (UNDEFND);                                                           // ... then do not compute the integral, and quit
# Note, can turn this off by setting `discard_integration_at_border = F`
# Although this setting makes absolute sense in whole genome setting, I have borders at every loci
# by design with MIPs (e.g. really multiplexed targeted amplicon sequencing). This means that I will have hard
# cutoffs for every vcfRobject. As such, I am going to ignore boundaries under the assumption that because the MIP density is all constant
# this gave all regions a "fair chance" to have that same border.
# It will certainly inflate samples/regions integration with



# Find the marker with the largest integrated EHHS statistic by Sabeti (unnormalized)
drugregions_sub$marker <- purrr::map(drugregions_sub$scanhh, function(x){
  maxmarker <- which(x$iES_Sabeti_et_al_2007 == max(x$iES_Sabeti_et_al_2007, na.rm = T))
  maxmarker <- ifelse(purrr::is_empty(maxmarker), NA, maxmarker)
  return(maxmarker)
})


# Throw away loci that don't have a clear marker
drugregions_sub <- drugregions_sub %>%
  dplyr::filter(!is.na(marker))




#.................................
# call ehh based on ihs value
#.................................
drugregions_sub$ehh <- map2(drugregions_sub$haplohh, drugregions_sub$marker, ~rehh::calc_ehh(
                            haplohh = .x,
                            mrk = .y,
                            limhaplo = 2,
                            limehh = 0.05,
                            discard_integration_at_border = F,
                            plotehh = F)
)


#.................................
# make plots for ehh continous by region
#.................................
crossehhplotdf <- tibble(geneid = drugregions_sub$geneid,
                    marker = unlist( drugregions_sub$marker ),
                    region = drugregions_sub$region,
                    name = drugregions_sub$name )



crossehhplotdf$pos <- map(drugregions_sub$haplohh, "position")
crossehhplotdf$ehh <- map(drugregions_sub$ehh, "ehh")
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
plotehh <- function(ehhdf, region, marker, geneid){
  plotObj <- ggplot() +
    geom_line(data=ehhdf, aes(x=pos, y=ehh, colour = factor(allele), group = factor(allele))) +
    scale_color_manual("Allele Status", values = c("#4575b4", "#d73027")) +
    ggtitle(paste0("This is marker: ", marker, " on gene ", geneid, " for region ", region)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10))
}

crossehhplotdf$ehhplot <-  pmap(crossehhplotdf[,c("region", "marker", "geneid", "ehhdf")], plotehh)

#.................................
# make haplotype plots
#.................................
crossehhplotdf$hapdf <- map(drugregions_sub$thap, function(x){
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
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_blank()
    )
}
crossehhplotdf$happlot <-  pmap(crossehhplotdf[,c("region", "name", "hapdf")], plothap)


#.................................
# write out
#.................................
save(drugregions_sub, crossehhplotdf, file = "data/derived/02-drugres_obj_rehh_subpopulationlevel.rda")

