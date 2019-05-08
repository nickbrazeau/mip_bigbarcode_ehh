#-------------------------------------------------------------------------------------------------------
# EHH Code Manipulations (not included in vcfRmanip)
#-------------------------------------------------------------------------------------------------------

getpmap.polarize <- function(x, chrommapdf = mapmodels){
  chrompos <- vcfR::getFIX(x)[, c("CHROM", "POS")] %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(pos = as.numeric(POS),
                  chromosome = rplasmodium::factor_chrom(CHROM),
                  start = pos,
                  end = pos)

  chrommapdf <- chrommapdf %>%
    dplyr::rename(chromosome = chr)



  predmp <- fuzzyjoin::genome_inner_join(x = chrompos,
                                         y = chrommapdf,
                                         by = c("chromosome", "start", "end"))

  # sanity check
  if(nrow(predmp) != nrow(chrompos)){
    stop("Issue with the Fuzzy Genomic Join")
  }

  # local interpolation
  predmp <- predmp %>%
    mutate(cMpred = cM_slope*(pos-start.y) + cM_intercept)

  # sanity check
  ismonotonic <- predmp %>%
    group_by(CHROM) %>%
    dplyr::summarise( ismonotonic = MonoInc::monotonic(cMpred, direction = "inc")) %>%
    dplyr::select(ismonotonic) %>%
    base::unlist(.)
  if(any(! ismonotonic )){
    stop("Your cM predictions are not monotonically increasing along the genome")
  }


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
