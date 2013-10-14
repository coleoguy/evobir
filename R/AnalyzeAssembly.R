#this code is just for basic stats about a sequenced genome
AnalyzeAssembly <- function(genome, max_N=25, plot=F){
  test<-unlist(lapply(genome, nchar))                     #record length of all seqs
  scaf.gsz <- sum(test)                                   #genome size
  big.scaf <- sum(test>1000000)                           #number of scaffolds greater than 1MB
  #calc the gc content
  gc <- vector()
  atcg <- vector()
  for(i in 1:length(genome)){
    gc[i] <- str_count( genome[[i]], "c") + str_count( genome[[i]], "g")
    atcg[i] <- nchar(genome[[i]])
  }
  gc.cont <- sum(gc)/sum(atcg)*100
  
  
  test<-(sort(test, decreasing = T))
  foo<-0
  i<-1
  while(foo<(scaf.gsz*.5)){
    foo<-foo+test[i]
    i<-i+1
  }
  n50.scaf <- as.numeric(test[i])                    #N50 size
  
  #accumulation plot
  coverage<-vector("numeric",length(genome))
  for(i in 1:length(genome)){
    coverage[i]<-sum(test[1:i])
  }
  if(plot==T){
    plot(1:length(genome),coverage,pch=19,cex=.3,xlab="scaffolds",col="red",main="Cummulative Coverage" )
  }

  #create an assembly of contigs
  genome.cont <- list()
  z<-1
  gap<-paste(rep("n",max_N),sep="",collapse="")
  
  for(i in 1:length(genome)){
    foo2 <- strsplit(genome[[i]],gap)[[1]]
    foo2 <- foo2[nchar(foo2)>1]
    for(j in 1: length(foo2)){
      genome.cont[z] <- gsub("^n*?[acgt]","",foo2[j])
      z<-z + 1
    }
  }
  test.cont <- unlist(lapply(genome.cont, nchar))                     #record length of all seqs
  test.cont <- (sort(test.cont, decreasing = T))
  cont.gsz <-  sum(test.cont) 
  foo<-0
  i<-1
  while(foo<(cont.gsz*.5)){
    foo<-foo+test.cont[i]
    i<-i+1
  }
  n50.cont <- as.numeric(test.cont[i])                    #N50 size
  results <- data.frame()
  results[1,1] <- "SCAFFOLDS"
  results[1,2] <- ""
  results[2,1] <- "Number of Scaffolds:"
  results[2,2] <- length(genome)
  results[3,1] <- "Assembly Size Based on Scaffolds: "
  results[3,2] <- paste(round(scaf.gsz/1000000,3)," MB", sep="")
  results[4,1] <- "Number of Scaffolds over 1MB:"
  results[4,2] <- big.scaf
  results[5,1] <- "N50 Scaffold Size:"
  results[5,2] <- n50.scaf
  results[6,1] <- "CONTIGS"
  results[6,2] <- ""
  results[7,1] <- "Number of Contigs:"
  results[7,2] <- length(genome.cont)
  results[8,1] <- "Assembly Size Based on Contigs: "
  results[8,2] <- paste(round(cont.gsz/1000000,3)," MB", sep="")
  results[9,1] <- "N50 Contig Size:"
  results[9,2] <- n50.cont
  results[10,1] <- "Minimum Contig Size:"
  results[10,2] <- min(test.cont)
  results[11,1] <- "Percent GC:"
  results[11,2] <- round(gc.cont,2)
  colnames(results) <- c('Statistic', 'Value')
return(results)
}
  