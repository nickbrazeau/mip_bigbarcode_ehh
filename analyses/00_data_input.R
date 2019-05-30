#.................................
# Dependencies
#.................................
source("R/00-recover_biallelics.R")
devtools::install_github("andrewparkermorgan/rplasmodium")
library(rplasmodium)
devtools::install_github("mrc-ide/MIPanalyzer")
library(MIPanalyzer)
library(vcfR)
library(tidyverse)
library(stringr)


#.................................
# read in filtered big-barcode panel
#.................................
mipbigpanel <- readRDS("data/raw_snps_filtered/biallelic_distances.rds")
#.................................
# read in drug res panel
#.................................
mipDRpanel <- readRDS("data/raw_snps_filtered/dr_monoclonal.rds")
#.................................
# recover biallelic snps in the subsetted dr panel
#.................................
mipDRpanel <- recoverbiallelics(mipDRpanel)$mipanalyzer_multiallelic


biallelic_sites <- !stringr::str_detect(mipDRpanel$loci$ALT, ",")  # this works because we only have SNPs for now
# long way to subset to a biallelic from multiallelic vcf but checked in raw vcf on command line, this works fine for now
mipDRpanel$loci <- mipDRpanel$loci[biallelic_sites, ]
mipDRpanel$coverage <- mipDRpanel$coverage[,biallelic_sites]
mipDRpanel$counts <- mipDRpanel$counts[1,,biallelic_sites] # this is the wsraf which is typically stored in the biallelic mip analyzer object
# note in this RDS, NA is both 0 and missing. Will set to 0 for referent allele freq if there is coverage at that site
mipDRpanel$counts[ is.na(mipDRpanel$counts) ] <- 0 # first make all NA sites to 0
mipDRpanel$counts[ is.na(mipDRpanel$coverage) ] <- NA # then if there is no coverage at that site, write over to NA since it should be coded as missing

# TODO -- ask Bob to fix for consistency with big-barcode
# raw example vs fixed
mipDRpanel_raw <- readRDS("data/raw_snps_filtered/dr_monoclonal.rds")
mipDRpanel_raw$loci[241,]
mipDRpanel_raw$counts[,,241]
mipDRpanel_raw$coverage[,241]


fxexample <- which(paste0(mipDRpanel_raw$loci[241,1],mipDRpanel_raw$loci[241,2]) == paste0(mipDRpanel$loci$CHROM, mipDRpanel$loci$POS))
mipDRpanel$counts[,fxexample]
mipDRpanel$coverage[,fxexample]



class(mipDRpanel) <- "mipanalyzer_biallelic"



# check to see that all samples are present in both
if(
  sum( mipDRpanel$samples$ID %in% mipbigpanel$samples$ID ) !=
  nrow(mipDRpanel$samples)){
  stop("Incompatibility between big panel vcf and dr panel vcf")
}

# eyeball test for smpl agreement
mipDRpanel$samples$ID[ mipDRpanel$samples$ID %in% mipbigpanel$samples$ID ]
mipbigpanel$samples$ID[ mipbigpanel$samples$ID %in% mipDRpanel$samples$ID ]

# SUBSET big panel vcf
mipbigpanel_sub <- MIPanalyzer::filter_samples(mipbigpanel, c( mipbigpanel$samples$ID %in% mipDRpanel$samples$ID ),
                                               description = "Drop BB Samples to DR filtered Samples")


#.................................
# error handle the overlapping sites
#.................................

if(any(
  paste0(mipbigpanel$loci[,1], mipbigpanel$loci[,2]) %in% paste0(mipDRpanel$loci[,1], mipDRpanel$loci[,2])
)){
  warning(paste("Overlapping loci in the big mip panel and DR panel -- there are", sum( paste0(mipbigpanel$loci[,1], mipbigpanel$loci[,2]) %in% paste0(mipDRpanel$loci[,1], mipDRpanel$loci[,2]) ),
             "overlapping sites between the DR panel and Big Panel VCF"))
}

mipBB_sub_repeats_loci <- which( paste0(mipbigpanel_sub$loci[,1], mipbigpanel_sub$loci[,2]) %in% paste0(mipDRpanel$loci[,1], mipDRpanel$loci[,2]) )
mipDR_sub_repeats_loci <- which( paste0(mipDRpanel$loci[,1], mipDRpanel$loci[,2]) %in%  paste0(mipbigpanel_sub$loci[,1], mipbigpanel_sub$loci[,2]) )

# confirm they are arranged in the same way
mipbigpanel_sub$loci[mipBB_sub_repeats_loci, c("CHROM", "POS")]
mipDRpanel$loci[mipDR_sub_repeats_loci, c("CHROM", "POS")]

all_equal(mipbigpanel_sub$loci[mipBB_sub_repeats_loci, c("CHROM", "POS")] ,  mipDRpanel$loci[mipDR_sub_repeats_loci, c("CHROM", "POS")])

#...........
# liftover count
#............

# merge the big barcode panel reads into DR reads
tempcounts <- mipbigpanel_sub$counts[, mipBB_sub_repeats_loci]
missing <- is.na(tempcounts)  +  is.na(mipDRpanel$counts[, mipDR_sub_repeats_loci])

# temporarily move all missing to 0 for these loci
tempcounts[ is.na(tempcounts) ] <- 0
mipDRpanel$counts[, mipDR_sub_repeats_loci][ is.na(mipDRpanel$counts[, mipDR_sub_repeats_loci]) ] <- 0

# add them up
mipDRpanel$counts[, mipDR_sub_repeats_loci] <- mipDRpanel$counts[, mipDR_sub_repeats_loci] + tempcounts
# now move counts back to missing
mipDRpanel$counts[ , mipDR_sub_repeats_loci][missing == 2] <- NA


#...........
# liftover count
#............

tempcov <- mipbigpanel_sub$coverage[, mipBB_sub_repeats_loci]
missing <- is.na(tempcov)  +  is.na(mipDRpanel$coverage[, mipDR_sub_repeats_loci])

# temporarily move all missing to 0 for these loci
tempcov[is.na(tempcov)] <- 0
mipDRpanel$coverage[, mipDR_sub_repeats_loci][ is.na(mipDRpanel$coverage[, mipDR_sub_repeats_loci]) ] <- 0

# add them up
mipDRpanel$coverage[, mipDR_sub_repeats_loci] <- mipDRpanel$coverage[, mipDR_sub_repeats_loci] + tempcov
# now move counts back to missing
mipDRpanel$coverage[ , mipDR_sub_repeats_loci][missing == 2] <- NA



## sanity check
if(sum(is.na(mipDRpanel$counts[ , mipDR_sub_repeats_loci])) != sum(is.na(mipDRpanel$coverage[ , mipDR_sub_repeats_loci]))){
  stop("Issue with additive merge of big-barcode and drug-res counts/coverage")
}


# drop these sites in big barcode now
mipbigpanel_sub$counts <- mipbigpanel_sub$counts[, -c(mipBB_sub_repeats_loci)]
mipbigpanel_sub$coverage <- mipbigpanel_sub$coverage[, -c(mipBB_sub_repeats_loci)]
mipbigpanel_sub$loci <- mipbigpanel_sub$loci[-c(mipBB_sub_repeats_loci), ]



#.................................
# Make the vcfs
#.................................
mipbivcfR <- MIPanalyzerbi2vcfR(input = mipbigpanel_sub, cutoff = 0.5) # given that these were identified as monoclonal, going to force to major haplotype
mipDRvcfR <- MIPanalyzerbi2vcfR(input = mipDRpanel, cutoff = 0.5)

#.................................
# remove sites that are all ref in the BB-DR combined panel
#.................................
mipdrgt <- vcfR::extract.gt(mipDRvcfR, element = "GT")
mipdr_ref_loci <- apply(mipdrgt, 1, function(x){ all(x == "0/0", na.rm = T) })

mipDRvcfR@fix <- mipDRvcfR@fix[!mipdr_ref_loci, ]
mipDRvcfR@gt <- mipDRvcfR@gt[!mipdr_ref_loci, ]

mipbbgt <- vcfR::extract.gt(mipbivcfR, element = "GT")
mipbb_ref_loci <- apply(mipbbgt, 1, function(x){ all(x == "0/0", na.rm = T) })

mipbivcfR@fix <- mipbivcfR@fix[!mipbb_ref_loci, ]
mipbivcfR@gt <- mipbivcfR@gt[!mipbb_ref_loci, ]

#.................................
# write out to be compatible with pf3d7
#.................................
liftover <- tibble(chrom = rplasmodium::chromnames(genome = "pf3d7")[1:14],
       chr = paste0("chr", seq(1:14)))
CHROMbb <- left_join(tibble(chr = vcfR::getCHROM(mipbivcfR)), liftover)
mipbivcfR@fix[,1] <- unlist(CHROMbb[,2])

CHROMdr <- left_join(tibble(chr = vcfR::getCHROM(mipDRvcfR)), liftover)
mipDRvcfR@fix[,1] <- unlist(CHROMdr[,2])

vcfR::write.vcf(mipbivcfR, file = "data/derived/vcfs/mipbi_bigbarcodepanel.vcf.gz")
vcfR::write.vcf(mipDRvcfR, file = "data/derived/vcfs/mipbi_drugrespanel.vcf.gz")

system("bash ~/Documents/Github/mip_bigbarcode_ehh/analyses/vcf_manip_vcfdo.sh")





