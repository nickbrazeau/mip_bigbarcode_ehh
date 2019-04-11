#.................................
# Dependencies
#.................................
devtools::install_github("nickbrazeau/rplasmodium")
library(rplasmodium)
library(vcfR)
library(tidyverse)
library(stringr)
library(GenomicRanges)

#--------------------------------------------------------------
# Metadata files
#--------------------------------------------------------------

# monoclonal samples
mc <- readRDS("data/monoclonals.rds")
# drug res sites
drugres <- rplasmodium::pf_3d7_PutDrugRxSites
drugres$start <- drugres$start + 1 # take 1-based, vcfs are 1-based

#.................................
# make map file
#.................................
maplm <- function(data){
  lm(cM ~ pos,
     data = data)
}
mp <- rplasmodium::recomb_map_pf3d7
mapmodels <- mp %>%
  group_by(chr) %>%
  nest() %>%
  dplyr::mutate(model = purrr::map(data, maplm),
                chr_fct = rplasmodium::factor_chrom(chr),
                CHROM = chr
                ) %>%
  dplyr::select(-c("chr_fct", "chr"))


#.................................
# alter drug res
#.................................
drugres <- drugres %>%
  dplyr::filter(chr != "Pf_M76611")  # no recombo in mtdna, no EHH

## Find Overlap
drugres_10kb <- GRanges(seqnames = drugres$chr,
                        ranges = IRanges(start = drugres$start - 1e4,
                                         end = drugres$end + 1e4 ),
                        strand = factor(rep("+", nrow(drugres))),
                        mcols=drugres$geneid) # make genomic ranges object

drugres_50kb <- GRanges(seqnames = drugres$chr,
                        ranges = IRanges(start = drugres$start - 5e4,
                                         end = drugres$end + 5e4 ),
                        strand = factor(rep("+", nrow(drugres))),
                        mcols=drugres$geneid) # make genomic ranges object
drugres_100kb <- GRanges(seqnames = drugres$chr,
                        ranges = IRanges(start = drugres$start - 1e5,
                                         end = drugres$end + 1e5 ),
                        strand = factor(rep("+", nrow(drugres))),
                        mcols=drugres$geneid) # make genomic ranges object

# read in MIP design
ideelbar <- read_tsv(file = "~/Google_Drive/PhD_Work/IDEEL_Colab/Project_DRC_BIG_BARCODE/IDEELBarcode_MIP_marker_SNPs/bedLocations/allPf3d7_intersectedWithGenes.bed", col_names = F)
colnames(ideelbar) <- c("CHROM", "START", "END", "EXTRACT", "TARGETLENGTH", "SENSE", "INFO")
MIP_ranges <- GRanges(seqnames = ideelbar$CHROM,
                      ranges = IRanges(start = ideelbar$START,
                                       end = ideelbar$END),
                      strand = factor(rep("+", nrow(ideelbar))),
                      mcols=ideelbar$INFO)


overlap10KB <- GenomicRanges::countOverlaps(drugres_10kb, MIP_ranges) # find these overlaps
overlap50KB <- GenomicRanges::countOverlaps(drugres_50kb, MIP_ranges) # find these overlaps
overlap100KB <- GenomicRanges::countOverlaps(drugres_100kb, MIP_ranges) # find these overlaps
overlapdf <- cbind(drugres[, c("geneid", "name")],
                   overlap10KB, overlap50KB, overlap100KB)


## Make final drug res table
## after talking to group 100 kb was selected
drugres <- drugres %>%
  dplyr::mutate(chrom_fct = rplasmodium::factor_chrom(chr),
                seqname = paste0("chr", as.character(chrom_fct)),
                start = start - 1e4, # already made 1-based
                # TODO fix this
                end = end  + 1e5)
# remember for the vcfR2SubsetChromPos function to work the chromosome name in the VCF must be in the `seqname` field. To protect against multiple chrom names floating around
drugres$seqname <- drugres$chr # updated for our work

# make sure we are still on the genomic map
drugres$start[ min(drugres$start) < 0 ] <- 0
drugends <- aggregate(drugres$end, list(factor(drugres$chr)), max)
colnames(drugends)[1] <- "CHROM"
chromends <- tibble(chr = names( rplasmodium::chromsizes_3d7()),
                    chromend =  rplasmodium::chromsizes_3d7())
( drugends <- chromends %>%
    dplyr::rename(CHROM = chr) %>%
    dplyr::left_join(drugends, .) %>%
  dplyr::mutate(offend = x > chromend) )
if(any(drugends$offend)){
  warning("Have mapped beyone the end of the chromosome; automatic fix")

  drugres <- drugres %>%
    dplyr::left_join(x=., y=chromends, by =  "chr") %>%
    dplyr::mutate(end = ifelse(end > chromend, chromend, end)) %>%
    dplyr::select(-c("chromend", "x"))
}

# but false for crt and dhps so OK (on chrom 7 and 8)




