#' @param
#' @details The purpose of this function is to recover sites that were previously mutliallelic in a larger
#' dataset but have no become biallelic in a "subset". This will look at a \code{mipanalyzer_multiallelic}
#' object and try to recover sites that can now be considered biallelic before making a \code{mipanalyzer_biallelic} object

recoverbiallelics <- function(mpnlyzr_multi){

  if(!inherits(mpnlyzr_multi, "mipanalyzer_multiallelic")){
    stop(paste("This is not a mipanalyzer_multiallelic object, it is a", class(mipanalyzer_multiallelic), "object"))
  }

  mlti_sites <- data.frame(
    mlti_sites = which(stringr::str_count(mpnlyzr_multi$loci$ALT, ",") > 0),
    biallelic = NA, # truly biallelic logical
    altbase = NA # which alt base is the biallelic alt allele
  )

  # going to be slower but for clarity, putting in for-if loop
  nloci <- nrow(mlti_sites)
  nsmpls <- nrow(mpnlyzr_multi$samples)


  for(i in 1:nloci){
    # 4 bases but only 3 ALTs as REF (1) can never be all NA
    alt1 <- !all( is.na(mpnlyzr_multi$counts[2, , mlti_sites$mlti_sites[i]]) )
    alt2 <- !all( is.na(mpnlyzr_multi$counts[3, , mlti_sites$mlti_sites[i]]) )
    alt3 <- !all( is.na(mpnlyzr_multi$counts[4, , mlti_sites$mlti_sites[i]]) )

    if(sum(c(alt1, alt2, alt3)) == 1){

      mlti_sites$biallelic[i] <-  TRUE
      # find biallelic base
      bsnum <- which(c(alt1, alt2, alt3))

      # find alt site and overwrite it with the new biallelic base
      mlti_sites$altbase[i] <- unlist( stringr::str_split(mpnlyzr_multi$loci$ALT[mlti_sites$mlti_sites[i]], ",") )[ bsnum ]

      mpnlyzr_multi$loci$ALT[ mlti_sites$mlti_sites[i] ] <- mlti_sites$altbase[i]

    } else if(sum(c(alt1, alt2, alt3)) == 3){
      # all REF
      # Note, this line doesn't need to be explicit but
      # just noting this scenario as this site is no longer variable (and therefore not a biallelic)
      # user will need to drop other non-variable sites explicitly
      mlti_sites$biallelic[i] <-  FALSE
      warning(paste("Loci",
                    paste0(
                      mpnlyzr_multi$loci$CHROM[mlti_sites$mlti_sites[i]], "_",  mpnlyzr_multi$loci$POS[mlti_sites$mlti_sites[i]]),
                    "in this set has no variable sites, and should be dropped. \n"))

    } else{
      mlti_sites$biallelic[i] <-  FALSE
    }

  } # end for loop

  return(
    list(
      mlti_sites_ret = mlti_sites,
      mipanalyzer_multiallelic = mpnlyzr_multi
    )
  )


}

