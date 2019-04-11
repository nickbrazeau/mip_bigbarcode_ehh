#-------------------------------------------------------------------------------------------------------
# EHH Code Manipulations (not included in vcfRmanip)
#-------------------------------------------------------------------------------------------------------

getpmap.polarize <- function(x, chrommapdf = mapmodels){
  chrompos <- vcfR::getFIX(x)[, c("CHROM", "POS")] %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(pos = as.numeric(POS))

  predmp <- inner_join(chrommapdf, chrompos) %>%
    mutate(cMpred = map2_dbl(model, pos, ~predict(.x, newdata = data.frame(pos = .y))))

  fixtidy <- vcfR::extract_info_tidy(x, info_fields = "AA")
  refalt <- tibble::as_tibble(vcfR::getFIX(x)[,c("REF", "ALT")])
  alleles <- dplyr::bind_cols(refalt, fixtidy)

  ancderiv <- tibble::tibble(
    anc = alleles$AA,
    der = ifelse(alleles$REF == alleles$AA,
                 alleles$ALT, alleles$REF)

  )
  # APM specific masks
  ancderiv$anc[fixtidy$AA %in% c("N", "X")] <- NA
  ancderiv$der[fixtidy$AA %in% c("N", "X")] <- NA

  ret <- tibble::tibble(snpname = paste0(predmp$CHROM, "_", predmp$pos),
                        chrom = predmp$CHROM,
                        pos = predmp$cMpred,
                        anc = ancderiv$anc,
                        der = ancderiv$der
  )

  return(ret)
}



getpmap.polarize.nucbp <- function(x){
  chrompos <- vcfR::getFIX(x)[, c("CHROM", "POS")] %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(pos = as.numeric(POS))


  fixtidy <- vcfR::extract_info_tidy(x, info_fields = "AA")
  refalt <- tibble::as_tibble(vcfR::getFIX(x)[,c("REF", "ALT")])
  alleles <- dplyr::bind_cols(refalt, fixtidy)

  ancderiv <- tibble::tibble(
    anc = alleles$AA,
    der = ifelse(alleles$REF == alleles$AA,
                 alleles$ALT, alleles$REF)

  )
  # APM specific masks
  ancderiv$anc[fixtidy$AA %in% c("N", "X")] <- NA
  ancderiv$der[fixtidy$AA %in% c("N", "X")] <- NA

  ret <- tibble::tibble(snpname = paste0(chrompos$CHROM, "_", chrompos$pos),
                        chrom = chrompos$CHROM,
                        pos = chrompos$pos,
                        anc = ancderiv$anc,
                        der = ancderiv$der
  )

  return(ret)
}
