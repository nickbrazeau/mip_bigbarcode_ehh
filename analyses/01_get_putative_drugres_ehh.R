
#-------------------------------------------------------------------------------------------------------
# Purpose of this script is to perform EHH controlling for population structure
# as well as to subset to samples with no missingness. The rehh package stops
# "considering" a haplotype as soon as it hits its first missing base pair. This
# (logical) effect can affect the monotonicity of the EHH decay. Would have needed to impute
# (which is complicated given MOI affecting our AF/LD calculations and limited monoclonal samples) or drop them, so have dropped missing samples
#-------------------------------------------------------------------------------------------------------

#.................................
# imports
#.................................
source("analyses/00_metadata.R")
source("R/01-EHH_tools.R")
source("R/04-subset_smpls.R")
library(tidyverse)
devtools::install_github("IDEELResearch/vcfRmanip")
library(vcfRmanip)
library(rehh)

mipDRpanel <- readRDS("data/raw_snps_filtered/dr_monoclonal.rds") # just for checks later
mipbb_dr_panel_vcfR <- vcfR::read.vcfR(file = "data/derived/vcfs/polarized_mipbi_drugres_bigbarcode_panels.vcf.bgz")
mipbb_dr_panel_vcfR@gt[mipbb_dr_panel_vcfR@gt == "./.:.:."] <- NA # import issue. See "issue#1" in vcfdo https://github.com/IDEELResearch/vcfdo/issues/1 but resolved because issue lies with variant calling...either way, should be NA
#-------------------------------------------------------------------------------------------------------------------------------
# datawrangle
#-------------------------------------------------------------------------------------------------------------------------------
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

regions <- regions %>%
  dplyr::mutate(region = ifelse(region == "DRC1.2", "DRC-East",
                                ifelse(region == "DRC2.2", "DRC-West", as.character( region) )),
                region = factor(region))


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

# find the number of samples in each vcf
# TZ had no monocolonals, so drop that
drugregions$nsmpls <- unlist(purrr::map(drugregions$vcfRobj, function(x){ return(ncol(x@gt)-1) }))
drugregions <- drugregions %>%
  dplyr::filter(nsmpls > 0)


#.................................
# Drop Haplotypes with any missing site
# and recalculate n smpls and nvar
#.................................
drugregions_sub <- drugregions
drugregions_sub$vcfRobj_new <- purrr::map(drugregions$vcfRobj, remove_smpls_by_smpl_missingness, misscutoff = 0)
# recalculate nvar and nsmpl with new not missing vcfR
drugregions_sub$nvar <- unlist(purrr::map(drugregions_sub$vcfRobj_new, function(x){ return(nrow(x@gt)) }))
drugregions_sub$nsmpls <- unlist(purrr::map(drugregions_sub$vcfRobj_new, function(x){ return(ncol(x@gt)-1) }))


#.................................
# make map and haplotype files
#.................................
drugregions_sub$pmap <- purrr::map(drugregions_sub$vcfRobj_new, getpmap.polarize)
drugregions_sub$thap <- purrr::map(drugregions_sub$vcfRobj_new, vcfRmanip::vcfR2thap)

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





#---------------------------------------------------------------------------------------------------------------------------------
# find putiative markers markers
#---------------------------------------------------------------------------------------------------------------------------------
# This is the supplementary Table that Bob V put together
# Details the prevalence of putative drug resistance loci

chromliftover <- data.frame(CHROM = rplasmodium::chromnames()) %>%
  dplyr::filter(! CHROM %in% c("Pf_M76611", "Pf3D7_API_v3")) %>%
  dplyr::mutate(
    chrom = paste0("chr", 1:14))


putative_drugsites  <- read_csv("data/Supplemental Table 2_ Drug resistance MIP panel prevalences - Sheet1.csv") %>%
  dplyr::mutate(
    DRC_perc = stringr::str_split_fixed(DRC, "\\(", n=2)[,2],
    prev = as.numeric( stringr::str_extract(DRC_perc,  "\\d+\\.*\\d*") ) #https://stackoverflow.com/questions/19252663/extracting-decimal-numbers-from-a-string
  ) %>%
  dplyr::filter(prev > 5) %>% # subset to putative sites with reasonable support in the DRC
  dplyr::left_join(x=., y=chromliftover, by = "chrom") %>%
  dplyr::mutate(snpname = paste0(CHROM, "_", pos)) %>%
  dplyr::select(c("gene", "mut_name", "snpname")) %>%
  dplyr::mutate(gene = ifelse(gene == "dhfr-ts", "dhfr", gene),
                gene = ifelse(gene == "k13", "kelch", gene),
                gene = ifelse(gene == "atp6", "pfatp6", gene)) %>%
  dplyr::rename(name = gene)

if(!all(putative_drugsites$name %in% drugregions_sub$name)){
  stop("There is a conflict in the names in the putative drug resistance file (supp table 2) and the drug regions vcf map that you have made")
}

drugregions_sub <- drugregions_sub %>%
  dplyr::left_join(., putative_drugsites, by = "name") %>%
  dplyr::filter(!is.na(mut_name))


findmarker <- function(haplohh, snpname){
  marker <- which(haplohh@snp.name == snpname)
  if(purrr::is_empty(marker)){
    return(NA)
  } else{
    return(marker)
  }
}

findcMPosfromMarker <- function(haplohh, marker){
  cMPos <- haplohh@position[[marker]]
  return(cMPos)
}

#......................
# Get Markers
#......................
drugregions_sub$marker <- unlist( purrr::pmap(drugregions_sub[,c("haplohh", "snpname")], findmarker) )
# sanity check
check <- drugregions_sub %>%
  dplyr::filter(is.na(marker))
# multiallelic
mipDRpanel$loci[mipDRpanel$loci$CHROM %in% paste0("chr", check$chrom_fct[1]) & mipDRpanel$loci$POS == 404407, ] # manual check here since pos is in cM now
mipDRpanel$loci[mipDRpanel$loci$CHROM %in% paste0("chr", check$chrom_fct[2]) & mipDRpanel$loci$POS == 267882, ] # manual check here since pos is in cM now

#......................
# Subset marker (missing ones are non-biallelics)
#......................
drugregions_sub <- drugregions_sub %>%
  dplyr::filter(!is.na(marker))
#......................
# Get cMPos from marker (for spiderplots)
#......................
drugregions_sub$cMPos <- purrr::pmap(drugregions_sub[,c("haplohh", "marker")], findcMPosfromMarker)


#---------------------------------------------------------------------------------------------------------------------------------
# call ehh based on putative markers
#---------------------------------------------------------------------------------------------------------------------------------
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
                    marker = drugregions_sub$marker,
                    mutation = drugregions_sub$mut_name,
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
plotehh <- function( region, mutation, name, ehhdf){
  plotObj <- ehhdf %>%
    dplyr::mutate(allele = factor(allele),
                  allele = forcats::fct_drop(allele))

  chck <- plotObj %>%
    dplyr::group_by(allele) %>%
    dplyr::summarise( checknotflat = all(ehh == 0)) %>%
    dplyr::filter(checknotflat) %>%
    dplyr::select(allele) %>%
    unlist(.)

  if(!purrr::is_empty(chck)){
    plotObj <- plotObj %>%
      dplyr::filter(! allele %in% chck)
  }

  plotObj %>%
    ggplot() +
    geom_line(aes(x=pos, y=ehh, colour = allele, group = allele), alpha = 0.8) +
    scale_color_manual("Allele Status", values = c("#4575b4", "#d73027")) +
    ggtitle(paste0("EHH Decay Plot for ", mutation, " in ", region)) +
    labs(x = "position (cM)", y = "EHH Statistic") +
    scale_color_manual("Allele", labels = c("Ancestral", "Derived"), values = c("#2166ac", "#b2182b")) +
    theme_bw() +
    theme(plot.title = element_text(family = "Arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 12),
          axis.title = element_text(family = "Arial", face = "bold", size = 10),
          axis.text.x =  element_text(family = "Arial", face = "bold", size = 8),
          legend.title = element_text(family = "Arial", face = "bold", hjust = 0.5, vjust = 0.5, size = 10))
}

crossehhplotdf$ehhplot <-  pmap(crossehhplotdf[,c("region", "mutation", "name", "ehhdf")], plotehh)


#.................................
# make haplotype plots
#.................................
crossehhplotdf$hapdf <- pmap(drugregions_sub[, c("thap", "pmap", "haplohh", "marker")],
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


plothap <- function(hapdf, region, mutation){
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
    geom_text(data = asterickmarks, aes(x=loci, y=height, label = ast), size = 8, color = "#ff7f00") +
    facet_grid(allelepol ~ ., scale = "free_y", space = "free_y") +
    scale_fill_manual("Allele", values = c("#d73027", "#fee090", "#66bd63", "#4575b4", "#f0f0f0")) +
    ggtitle(paste0("Haplotype Plot for ", mutation, " in ", region)) +
    ylab("Samples") + xlab("Loci") +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 10),
          axis.title = element_text(face = "bold", hjust = 0.5, size = 8.5),
          axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 3),
          axis.text.y = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 3)
    )
}
crossehhplotdf$happlot <-  pmap(crossehhplotdf[,c("hapdf", "region", "mutation")], plothap)


#.................................
# write out
#.................................
save(drugregions_sub, crossehhplotdf, file = "data/derived/01-drugres_obj_rehh_subpopulationlevel.rda")






