#' @param misscutoff numeric; This is the proportion of missing that you will allow. Any more missing and the sample will be excluded

remove_smpls_by_smpl_missingness <- function(vcfRobject, misscutoff){

  miss <- vcfRmanip::calc_loci_missingness_by_smpl(vcfRobject)
  if(misscutoff > 0){
    miss <- miss %>%
      dplyr::filter(missprop < misscutoff) %>%
      dplyr::select(c("sample")) %>%
      unlist(.)
  } else if(misscutoff == 0){
    miss <- miss %>%
      dplyr::filter(missprop <= misscutoff) %>%
      dplyr::select(c("sample")) %>%
      unlist(.)
  }

  vcfRobject_new <- vcfRmanip::select_samples(vcfRobject = vcfRobject,
                                              smplvctr = miss)

  return(vcfRobject_new)

}


