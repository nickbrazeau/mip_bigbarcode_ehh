

#' ... for other argument passed to rehh::calc_ehhs

unnormalized_xpehh_derivedallele_logratio <- function(haplohh.pop1, haplohh.pop2, marker, discard_integration_at_border = F, limhaplo = 2, limehh = 0.05, ...){

  ehh_d.pop1 <- rehh::calc_ehh(haplohh.pop1, mrk = marker, plotehh = F,
                              discard_integration_at_border = discard_integration_at_border,
                              limhaplo = limhaplo,
                              limehh = limehh, ...)

  ehh_d.pop2 <- rehh::calc_ehh(haplohh.pop2, mrk = marker, plotehh = F,
                              discard_integration_at_border = discard_integration_at_border,
                              limhaplo = limhaplo,
                              limehh = limehh, ...)

  base::message("The logratio of Pop2/Pop1 is considered using the IHH Derived Allele statistic")

  ret <- log( ehh_d.pop2$ihh["Derived allele"] / ehh_d.pop1$ihh["Derived allele"] )

  return(ret)
}


# basic function for taking a pairwise distance matrix to a square distance matrix
pairwisedistdf_to_squaredist <- function(pairwisedistdf){
  pairwisedistdf <- as.matrix(pairwisedistdf)
  ret <- tapply(pairwisedistdf[,3], list(pairwisedistdf[,1],pairwisedistdf[,2]), "[[", 1)
  ret <- stats::as.dist(ret)
  return(ret)
}
