#.................................
# Dependencies
#.................................
devtools::install_github("nickbrazeau/rplasmodium")
library(rplasmodium)
devtools::install_github("OJWatson/MIPanalyzer")
library(MIPanalyzer)
library(vcfR)
library(tidyverse)
library(stringr)


#.................................
# read in filtered big-barcode panel
#.................................
mipbigpanel <- readRDS("data/biallelic_distances.rds")
#.................................
# read in drug res panel
#.................................
mipDRpanel <- readRDS("data/dr_processed.rds")
biallelic_sites <- !stringr::str_detect(mipDRpanel$loci$ALT, ",")  # this works because we only have SNPs for now
# long way to subset to a biallelic from multiallelic vcf but checked in raw vcf on command line, this works fine for now
mipDRpanel$loci <- mipDRpanel$loci[biallelic_sites, ]
mipDRpanel$coverage <- mipDRpanel$coverage[,biallelic_sites]
mipDRpanel$counts <- mipDRpanel$counts[2,,biallelic_sites] # this is the wsnraf which is typically stored in the biallelic mip analyzer object
class(mipDRpanel) <- "mipanalyzer_biallelic"



# know from conversations w/ Bob that there the mipDRpanel is the limiting sample set
# will subset the big panel sample list to the DR panel list so that we can concatenate
if(
  sum( mipDRpanel$samples$ID %in% mipbigpanel$samples$ID ) !=
  nrow(mipDRpanel$samples)){
  stop("Incompatibility between big panel vcf and dr panel vcf")
}
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

# merge the big barcode panel reads into DR reads
tempcounts <- mipbigpanel_sub$counts[, mipBB_sub_repeats_loci]
tempcounts[is.na(tempcounts)] <- 0
mipDRpanel$counts[, mipDR_sub_repeats_loci] <- mipDRpanel$counts[, mipDR_sub_repeats_loci] + tempcounts

tempcov <- mipbigpanel_sub$coverage[, mipBB_sub_repeats_loci]
tempcov[is.na(tempcov)] <- 0
mipDRpanel$coverage[, mipDR_sub_repeats_loci] <- mipDRpanel$coverage[, mipDR_sub_repeats_loci] + tempcov

# drop these sites in big barcode now
mipbigpanel_sub$counts <- mipbigpanel_sub$counts[, -c(mipBB_sub_repeats_loci)]
mipbigpanel_sub$coverage <- mipbigpanel_sub$coverage[, -c(mipBB_sub_repeats_loci)]
mipbigpanel_sub$loci <- mipbigpanel_sub$loci[-c(mipBB_sub_repeats_loci), ]



#.................................
# Make the vcfs
#.................................

mipbivcfR <- MIPanalyzer::MIPanalyzerbi2vcfR(input = mipbigpanel_sub, cutoff = 0.1)
mipDRvcfR <- MIPanalyzer::MIPanalyzerbi2vcfR(input = mipDRpanel, cutoff = 0.1)

#.................................
# remove sites that are all ref in the 112 sample for BB and DR
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

vcfR::write.vcf(mipbivcfR, file = "data/mipbi_bigbarcodepanel.vcf.gz")
vcfR::write.vcf(mipDRvcfR, file = "data/mipbi_drugrespanel.vcf.gz")

system("bash ~/Documents/MountPoints/mountedMeshnick/Projects/mip_bigbarcode_ehh/analyses/vcf_manip_vcfdo.sh")





