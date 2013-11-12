ChangeHaploSuff <- function(x, y, haplosuff.mut, loci.s) {   # Mutations that extend the SLR(sex limited region)   # x = population    y = result of function RadDose
  newhap.val <- runif(1, min = 0, max = haplosuff.mut)
  foobish <- x[[y[[2]]]]  #this pulls out the genome we are mutating
  which.locus <- sample(1:loci.s, 1)
  which.hap <- sample(c(13, 14), 1)
  foobish[which.hap, which.locus] <- newhap.val
  foobish -> x[[y[[2]]]]
  return(x)
}
