CalcD <- function(alignment = "alignment.fasta", 
                  sig.test="B",                                                    # options are "N", "B", "J"
                  block.size = 100,                                                # size of blocks to drop in jacknife
                  replicate=1000){
  # this function is used regardless of approach
  d.calc <- function(alignment){
    abba <- 0                                                                         #  set up my variables
    baba <- 0                                                                         #  set up my variables
    for(i in 1:ncol(alignment)){                                               #  run through all sites
      if(length(unique(alignment[, i])) == 2){                                 #  unique(c(p1,p2,p3,o))==2 aka biallelic
        if(alignment[1, i] != alignment[2, i]){                         #  p1 != p2   aka different resolutions in p1 and p2
          if(alignment[4, i] != alignment[3, i]){                       #  o != p3    durand says "less likely pattern due to seq. errors
            if(alignment[3, i] == alignment[1, i]) {baba <- baba + 1}   #  add to the count of baba sites
            if(alignment[2, i] == alignment[3, i]) {abba <- abba + 1}   #  add to the count of abba sites
          } 
        }
      }
    }
    d <- (abba - baba) / (abba + baba)   #what its all about
    results <- list()
    results[[1]] <- d
    results[[2]] <- abba
    results[[3]] <- baba
    return(results)
  }
  
  #### Test of empirical data
  alignment <- read.alignment(alignment, format = "fasta")                         #  read in the alignment
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
  }
  results <- d.calc(alignment.matrix)
  d <- results[[1]]
  abba <- results[[2]]
  baba <- results[[3]]
  
## THIS SECTION WILL CALCULATE THE P-VAL BASED ON BOOTSTRAPPING
## SITES ARE SAMPLED WITH REPLACEMENT TO MAKE A NEW DATASET OF
## OF EQUAL SIZE TO THE ORIGINAL DATASET THIS ALLOWS US TO CALCULATE
## THE STANDARD DEVIATION AND THUS A Z SCORE.
  if(sig.test=="B"){
    sim.d<-vector()
    foo <- ncol(alignment.matrix)
    sim.matrix<-matrix(,4,foo)
    cat("\nperforming bootstrap")
    for(k in 1:replicate){
      cat(".")
      sim.matrix[1:4,1:foo] <- alignment.matrix[1:4, sample(1:foo, replace=T)]
      results <- d.calc(sim.matrix)
      sim.d[k] <- results[[1]]
    }
    z <- abs(d-0/sd(sim.d))
    new.pval <- 2 * (1 - pnorm(z))
    ## NOW WE MAKE THE OUTPUTS  
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic / Z-score = ", d, " / ", z)
    cat("\n\nResults from ", replicate, "bootstraps")
    cat("\nSD D statistic =", sd(sim.d))
    cat("\nP-value (that D=0) = ",new.pval) #after Eaton and Ree 2013 
  }

## THIS SECTION WILL CALCULATE THE P-VAL BASED ON JACKKNIFING
## THE DATA IS REANALYZED WHILE DROPPING A PORTION OF THE DATA
## THE SIZE OF THE DROPPED PORTION IS DETERMINED BY THE BLOCK SIZE
## ARGUMENT THIS PROCEDURE ALLOWS US TO CALCULATE
## THE STANDARD DEVIATION AND THUS A Z SCORE. THIS APPROACH IS PARTICULARLY 
## IMPORTANT WHEN WE WILL BE USING DATA WHERE THE SNPs MAY BE IN LINKAGE WITH
## ONE ANOTHER
if(sig.test=="J"){
  
  #first lets test whether we are limited in the number of reps by block size
  max.rep <- ncol(alignment.matrix) - block.size
  if(block.size >= (ncol(alignment.matrix)/2)){
    stop(call. = F, paste("\nWith a block size of", block.size, 
              "and an alignment of", ncol(alignment.matrix), 
              "sites \nsome sites would never be included in the \nanalysis",
              "\n\nThe maximum block size is 1/2 the alignment length"))
  }
  if(max.rep > replicate){
    drop.pos <- seq.int(from=1, to=(max.rep-1), length.out=replicate)
    replicate2 <- replicate
  }
  sim.d<-vector()
  foo <- ncol(alignment.matrix)
  sim.matrix<-matrix(,4,foo-block.size)
  cat("\nperforming jackknife")
  for(k in 1:replicate2){  
    cat(".")
    sim.matrix[1:4,1:(foo-block.size-1)] <-alignment.matrix[1:4, -drop.pos[k]:-(drop.pos[k]+block.size)]
    results <- d.calc(sim.matrix)
    sim.d[k] <- results[[1]]
  }
  sd.sim.d <- round(sqrt(var(sim.d)),5)
  mn.sim.d <- round(mean(sim.d),5)
  new.pval <- 2*(pnorm(-abs(d/sd.sim.d)))
  ## NOW WE MAKE THE OUTPUTS  
  cat("\nSites in alignment =", ncol(alignment.matrix))
  cat("\nNumber of sites with ABBA pattern =", abba)
  cat("\nNumber of sites with BABA pattern =", baba)
  cat("\nD raw statistic", d)
  cat("\nZ-score = ", d/sd.sim.d)
  cat("\n\nResults from", replicate2, "jackknifes with block size of", block.size)
  cat("\nSD D statistic =", sd.sim.d)
  cat("\nP-value (that D=0) = ",new.pval) #after Eaton and Ree 2013 
}
  if(sig.test=="N"){
    cat("\nSites in alignment =", ncol(alignment.matrix))
    cat("\nNumber of sites with ABBA pattern =", abba)
    cat("\nNumber of sites with BABA pattern =", baba)
    cat("\n\nD raw statistic = ", d)
  }
}
CalcPopD <- function(alignment = "alignment.fasta"){
  ##  Now we have eqn. 2 from page 2240
  ##  input is an alignment the can take multiple sequences from each 
  ##  population of interest.  IMPORTANT MAKE SURE SEQUENCES ARE IN ORDER
  ##  P1, P2, P3, OUTGROUP!  Again we find the biallelic sites but now 
  ##  those biallelic sites need not be fixed and we will calculate frequencies
  ##  of SNP for each population.  The way the function is set up we do need to 
  ##  feed in an alignment where each sequence from a population has the same name:
  ##  pop1
  ##  AACCACAAGCCAGCTCAGCTACAG
  ##  pop1
  ##  TACAACAAGCGAGCTCAGCTACAG
  ##  pop1
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop2
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop2
  ##  GGCCACAAGCCAGCTCAGCTACAG
  ##  pop3
  ##  TACCACAAGCCAGCTCAGCTACAG
  ##  OUTGROUP
  ##  TACCAGGAGCCAGCTCTTCTACCC
  Mode <- function(x) {                                                                      #  i need a little mode function which R is lacking ugh
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }  
  alignment<-read.alignment(alignment, format="fasta")                                       #  read in the alignment
  alignment.matrix<-matrix(,length(alignment$nam),nchar(alignment$seq[[1]])+1)               #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i,2:ncol(alignment.matrix)]<-unlist(strsplit(alignment$seq[[i]],""))    #  fill in the matrix
  }
  alignment.matrix[,1]<-alignment$nam                                                        #  get those names into our matrix row names dont work :(
  groups<-unique(alignment$nam)
  p1 <- p2 <- p3 <- p4 <- 0                                                                  #  lets just set up the variable names from the durand paper
  numerator <- denominator <- 0
  useful<-0                                                                                  #  plus some of my own
  segregating<-0                                                                             #  plus some of my own  
  seg.pos<-F                                                                                 #  plus some of my own  
  for(i in 2:ncol(alignment.matrix)){                                                        #  run through all sites
    seg.pos<-F                                                                               #  reset this switch  
    if(length(unique(alignment.matrix[,i]))==2){                                             #  unique(c(p1,p2,p3,o))==2 aka biallelic
      A <- Mode(alignment.matrix[alignment.matrix[, 1] == groups[4], i])                     #  lets treat the more common variant in the outgroup as "A"
      B <- unique(alignment.matrix[,i])[unique(alignment.matrix[, i]) != A]                  #  not purposely obfuscating... the other variant in variable "B"
      if(B %in% unique(alignment.matrix[alignment.matrix[, 1] == groups[3], i])){            #  makes sure that we have at least some indication of an ABBA/BABA pattern
        if(length(unique(alignment.matrix[alignment.matrix[, 1] %in% groups[1:2], i])) == 2){  #  makes sure that we've got some different resolutions in the ingroups
          useful <- useful + 1                                                                 #  lets just keep track of how many sites are even useful 
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[1], i])) == 2) {seg.pos<-T}#  next 5 lines are a lame way of counting sites that are segregating
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[2], i])) == 2) {seg.pos<-T}#  vs those that are fixed another words is population sampling
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[3], i])) == 2) {seg.pos<-T}#  really of any value within the data set that we are examining
          if(length(unique(alignment.matrix[alignment.matrix[, 1] == groups[4], i])) == 2) {seg.pos<-T}
          if(seg.pos == T){segregating <- segregating + 1} 
          #print(segregating)
          p1 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[1], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[1], i])  #  freq of A snp in first population
          p2 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[2], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[2], i])  #  freq of A snp in second population
          p3 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[3], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[3], i])  #  freq of A snp in third population
          p4 <- (sum(alignment.matrix[alignment.matrix[, 1] == groups[4], i] == A))/length(alignment.matrix[alignment.matrix[, 1] == groups[4], i])  #  freq of A snp in outgroup population
          #  Durands explanation of eqn 2 is lacking... as least to my feable mind!
          #  it appears to me that as written p hat is actually the frequency of SNP "B" so....
          #  snap...  vindicated my interpretation matches that found in the supplemental material of the 
          #  heliconius genome paper supplement... too cool
          p1 <- 1-p1  #convert these over from proportion A to proportion B
          p2 <- 1-p2  #convert these over from proportion A to proportion B
          p3 <- 1-p3  #convert these over from proportion A to proportion B
          p4 <- 1-p4  #convert these over from proportion A to proportion B
          numerator <- ((1 - p1) * p2 * p3 * (1 - p4)) - (p1 * (1 - p2) * p3 * (1 - p4)) + numerator              #  build up our numerator sum
          denominator <- ((1 - p1) * p2 * p3 * (1 - p4)) + (p1 * (1 - p2) * p3 * (1 - p4)) + denominator          #  build up our denominator sum
        }
      }
    }
  }
  d <- numerator / denominator    #what its all about
  
  user.result <- list()
  user.result$d.stat <- d
  user.result$pval <- "HELP"
  user.result$align.length <- ncol(alignment.matrix) - 1
  user.result$useful.sites <- useful
  user.result$seg.sites <- segregating
  print(paste("Sites in alignment =", ncol(alignment.matrix) - 1))
  print(paste("Number of sites with ABBA or BABA patterns =", useful))
  print(paste("Number of ABBA or BABA sites that are still segregating in at least one population =", segregating))  
  print(paste("D statistic =", d))
}
CalcPartD <- function(alignment = "alignment.fasta", boot=F, replicate = 1000, alpha =.05){
  alignment <- read.alignment(alignment, format = "fasta")                          #  read in the alignment
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
  }
  abbaa <- babaa <- 0    ## d1                                                                       #  set up my variables
  ababa <- baaba <- 0    ## d2                                                                       #  set up my variables
  abbba <- babba <- 0    ## d12 
  for(i in 1:ncol(alignment.matrix)){                                               #  run through all sites
    if(length(unique(alignment.matrix[, i])) == 2){                                 #  unique(c(p1,p2,p3.1,p3.2,O))==2 aka biallelic
      if(alignment.matrix[1, i] != alignment.matrix[2, i]){                         #  p1 != p2   aka different resolutions in p1 and p2
        if(alignment.matrix[5, i] != alignment.matrix[3, i] | alignment.matrix[5, i] != alignment.matrix[4, i] ){#  o != p3.1 or o is !=p3.2   
          ## D1
          if(alignment.matrix[4, i] == alignment.matrix[5, i]){
            if(alignment.matrix[1, i] == alignment.matrix[5, i]){abbaa <- abbaa+1}
            if(alignment.matrix[2, i] == alignment.matrix[5, i]){babaa <- babaa+1} 
          }
          ## D2        
          if(alignment.matrix[3, i] == alignment.matrix[5, i]){
            if(alignment.matrix[1, i] == alignment.matrix[5, i]){ababa <- ababa+1}
            if(alignment.matrix[2, i] == alignment.matrix[5, i]){baaba <- baaba+1} 
          }
          ##D12
          if(alignment.matrix[3, i] == alignment.matrix[4, i]){
            if(alignment.matrix[1, i] == alignment.matrix[5, i]){abbba <- abbba+1}
            if(alignment.matrix[2, i] == alignment.matrix[5, i]){babba <- babba+1}
          }
        } 
      }
    }
  }
  d1 <- (abbaa - babaa) / (abbaa + babaa)   
  d2 <- (ababa - baaba) / (ababa + baaba)   
  d12 <- (abbba - babba) / (abbba + babba)
  
  ## THIS SECTION WILL CALCULATE THE P-VAL BASED ON BOOTSTRAPPING
  ## SITES ARE SAMPLED WITH REPLACEMENT TO MAKE A NEW DATASET OF
  ## OF EQUAL SIZE TO THE ORIGINAL DATASET THIS ALLOWS US TO CALCULATE
  ## THE STANDARD DEVIATION AND THUS A Z SCORE.
  
  if(boot==T){
    sim.d1<-vector()
    sim.d2<-vector()
    sim.d12<-vector()
    foo <- ncol(alignment.matrix)
    sim.matrix<-matrix(,5,foo)
    for(k in 1:replicate){      
      for(j in 1:5){sim.matrix[j,1:foo] <-sample(alignment.matrix[j,1:foo],replace=T)}
      ##NOW JUST RERUN OUR WHOLE ALGORITHM     
      t.abbaa <- t.babaa <- 0     ## d1
      t.ababa <- t.baaba <- 0     ## d2                                                                    #  set up my variables
      t.abbba <- t.babba <- 0     ## d12
      for(i in 1:ncol(sim.matrix)){                                                              #  run through all sites
        if(length(unique(sim.matrix[, i])) == 2){                                                #  unique(c(p1,p2,p3.1,p3.2,O))==2 aka biallelic
          if(sim.matrix[1, i] != sim.matrix[2, i]){                                              #  p1 != p2   aka different resolutions in p1 and p2
            if(sim.matrix[5, i] != sim.matrix[3, i] | sim.matrix[5, i] != sim.matrix[4, i] ){    #  o != p3.1 or o is !=p3.2   
              ## D1
              if(sim.matrix[4, i] == sim.matrix[5, i]){
                if(sim.matrix[1, i] == sim.matrix[5, i]){t.abbaa <- t.abbaa+1}
                if(sim.matrix[2, i] == sim.matrix[5, i]){t.babaa <- t.babaa+1} 
              }
              ## D2        
              if(sim.matrix[3, i] == sim.matrix[5, i]){
                if(sim.matrix[1, i] == sim.matrix[5, i]){t.ababa <- t.ababa+1}
                if(sim.matrix[2, i] == sim.matrix[5, i]){t.baaba <- t.baaba+1} 
              }
              ##D12
              if(sim.matrix[3, i] == sim.matrix[4, i]){
                if(sim.matrix[1, i] == sim.matrix[5, i]){t.abbba <- t.abbba+1}
                if(sim.matrix[2, i] == sim.matrix[5, i]){t.babba <- t.babba+1}
              } 
            }
          }
        }
      }
      sim.d1[k] <- (t.abbaa - t.babaa) / (t.abbaa + t.babaa)   
      sim.d2[k] <- (t.ababa - t.baaba) / (t.ababa + t.baaba)   
      sim.d12[k] <- (t.abbba - t.babba) / (t.abbba + t.babba)
    }
    sd.sim.d1 <- round(sqrt(var(sim.d1)),5)
    mn.sim.d1 <- round(mean(sim.d1),5)
    new.pval.d1 <- 2*(pnorm(-abs(d1/sd.sim.d1)))
    
    sd.sim.d2 <- round(sqrt(var(sim.d2)),5)
    mn.sim.d2 <- round(mean(sim.d2),5)
    new.pval.d2 <- 2*(pnorm(-abs(d2/sd.sim.d2)))
    
    sd.sim.d12 <- round(sqrt(var(sim.d12)),5)
    mn.sim.d12 <- round(mean(sim.d12),5)
    new.pval.d12 <- 2*(pnorm(-abs(d12/sd.sim.d12)))
    
    if(is.nan(d1)) d1<- "Error Missing Data"
    if(is.nan(d2)) d2<- "Error Missing Data"
    if(is.nan(d12)) d12<- "Error Missing Data"
    
    ## NOW WE MAKE THE OUTPUTS  
    cat("Sites in alignment =", ncol(alignment.matrix))
    cat("\nD1 sites with ABBAA/BABAA pattern =", abbaa,"/",babaa)
    cat("\nD2 sites with ABABA/BAABA pattern =", ababa,"/",baaba)
    cat("\nD12 sites with ABABA/BAABA pattern =", abbba,"/",babba)
    cat("\n\nD1 raw statistic / Z-score =", d1,"/",d1/sd.sim.d1)
    cat("\nD2 raw statistic / Z-score =", d2,"/",d2/sd.sim.d2)
    cat("\nD12 raw statistic / Z-score =", d12,"/",d12/sd.sim.d12)
    if(!(d1=="Error Missing Data")){
      cat("\n\nD1 Bootstrap Statistics: ")
      cat("SD = ", sd.sim.d1)
      cat("  P-val = ", new.pval.d1)
    }
    if(!(d2=="Error Missing Data")){
      cat("\nD2 Bootstrap Statistics: ")
      cat("SD = ", sd.sim.d2)
      cat("  P-val = ", new.pval.d2)
    }
    if(!(d12=="Error Missing Data")){
      cat("\nD12 Bootstrap Statistics: ")
      cat("SD = ", sd.sim.d12)
      cat("  P-val = ", new.pval.d12)
    }
    cat("\n\nBonferroni adjustment: alpha selected:",alpha," number of tests:",3,"\nSo P-value of less than ",round(alpha/3,4)," should be considered significant", sep="")
  }
  if(boot==F){
    cat("\nD1 sites with ABBAA/BABAA pattern =", abbaa,"/",babaa)
    cat("\nD2 sites with ABABA/BAABA pattern =", ababa,"/",baaba)
    cat("\nD12 sites with ABABA/BAABA pattern =", abbba,"/",babba)
    cat("\n\nD1 raw statistic =", d1)
    cat("\nD2 raw statistic =", d2)
    cat("\nD12 raw statistic =", d12)
    
  } 
}