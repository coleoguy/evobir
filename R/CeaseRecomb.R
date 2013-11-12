CeaseRecomb <- function(x, recom.mu, loci.s, evaluation) {   # Mutations that extend the SLR(sex limited region)   # x = population    y = result of function RadDose
  y <- sample(which(evaluation[,1]=="Male"), 1)
  foobish <- x[[y]]                         # this pulls out the genome we are mutating
  inv.size <- sample(1:recom.mu, 1)
  chrom.pick <- sample(11:12, 1)
  inv.start <- which.max(foobish[chrom.pick, ])
  if({inv.start + inv.size} > loci.s){
    foobish[chrom.pick, ] <- 0
  }else{    
    foobish[chrom.pick, 1:{inv.start + inv.size}] <- 0
    x[[y]] <- foobish
  }
  return(x)
}
