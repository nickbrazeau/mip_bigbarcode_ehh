
#.................................
# Functions
#.................................


#' @title Assign the GT calls from the NONREFERENT WSAF
#' @return A GT matrix
#' not exported

assignGTfromWSNRAF<- function(wsnraf, cutoff = 0.1){

  GT <- matrix(NA, dim(wsnraf)[1], dim(wsnraf)[2])

  GT <- ifelse(wsnraf > 1-cutoff, "1/1",
               ifelse(wsnraf < 0+cutoff, "0/0",
                      ifelse(!is.na(wsnraf), "0/1", NA)))

  return(GT)

}


#' @title Convert a MIPanalyzer Biallelic Object to a vcfR object
#' @param input is of class \code{mipanalyzer_biallelic}
#' @param cutoff is the within-sample non-referent allele frequency cutoff to transform your biallelic site to a genotype matrix.
#'
#' @export


MIPanalyzerbi2vcfR <- function(input = NULL, cutoff = cutoff){

  if(!inherits(input, c("mipanalyzer_biallelic"))){
    stop("This function only works on objects of class mipanalyzer_biallelic or mipanalyzer_multiallelic, not class ", class(mipobj))
  }

  # setup for gt
  wsnraf <- input$counts/input$coverage
  GT <- assignGTfromWSNRAF(wsnraf, cutoff = cutoff)
  ADref <- input$counts
  ADalt <- input$coverage - input$counts
  DP <- input$coverage

  gt <- t( ifelse(is.na(GT), NA, paste0(GT, ":", ADref, ",", ADalt, ":", DP)) ) # loci as rows, smpls as columns

  # append format column and sample names
  gt <- cbind(FORMAT = "GT:AD:DP", gt)
  colnames(gt)[2:ncol(gt)] <- input$samples$ID

  # getFix
  fix <- as.matrix(input$loci[,c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")]) # must be in this order and only these
  # https://stackoverflow.com/questions/15618527/why-does-as-matrix-add-extra-spaces-when-converting-numeric-to-character
  # not sure why trim is failing here. annoying
  fix[,2] <- stringr::str_replace_all(fix[,2], "\\s", "")
  # get meta
  # meta <- append(input$vcfmeta, "##MIPanalyzer=This vcf was filtered and modified by the MIPanalyzer R package")
  meta <- c("##fileformat=VCFv4.2")
  # write out new vcfRobj
  newvcfR <- new("vcfR", meta = meta, fix = fix, gt = gt)
  return(newvcfR)

}

