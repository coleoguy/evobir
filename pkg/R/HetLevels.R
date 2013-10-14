# This script will slide across all of the assemblies
# and record a 0,1,2,3 based on the level of inferred 
# heterozygosity in the genome.  My thoughts are that 
# the sex chromosome should jump right out at you.  The
# one confounding factor that is really obvious is that
# sequence depth on the X chromosome is going to be lower
# if you sequence males or males and females.
# 1 for A,T,C,G
# 2 for R,Y,S,W,K,M
# 3 for B,D,H,V
# 4 for N
HetLevels <- function(fasta){
  raw.assemb <- read.fasta(fasta, as.string=F)
  het.count <- matrix(,length(raw.assemb),3)
  colnames(het.count) <- c("length","polymorphic.sites","Ns")
    for(i in 1:length(raw.assemb)){
      site.count <- sum(nchar(raw.assemb[[i]]))
      foo <- sum(unlist(lapply(raw.assemb[[i]], '%in%', c("a", "c", "g", "t"))))
      foo.poly <- sum(unlist(lapply(raw.assemb[[i]], '%in%', c("r", "y", "s", "w", "k", "m","b", "d", "h", "v"))))
      foo.n <- sum(unlist(lapply(raw.assemb[[i]], '%in%', "n")))
      het.count[i,1] <- foo + foo.poly + foo.n
      het.count[i,2] <- foo.poly 
      het.count[i,3] <- foo.n
      cat("\nProcessing Sequence:",i)
    }
  cat("\n")
  return(het.count)
}
