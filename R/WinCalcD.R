WinCalcD <- function(alignment = "alignment.fasta", win.size = 100, step.size=50,
                     boot = F, replicate = 1000){
  alignment <- read.alignment(alignment, format = "fasta")                         #  read in the alignment
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))    #  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
  }
  full.align <- alignment.matrix
  
  total <- ncol(full.align)
  spots <- seq(from = 1, to = (total - win.size), by = step.size)
  results.matrix <- as.data.frame(matrix(,1,6))
  colnames(results.matrix) <- c("range", "abba", "baba", "d", "Z", "pval")
  for(q in 1:win.num){
    alignment.matrix <- full.align[spots[q]:(spots[q] + win.size - 1)]
    abba <- 0                                                                         #  set up my variables
    baba <- 0                                                                         #  set up my variables
    for(i in 1:ncol(alignment.matrix)){                                               #  run through all sites
      if(length(unique(alignment.matrix[, i])) == 2){                                 #  unique(c(p1,p2,p3,o))==2 aka biallelic
        if(alignment.matrix[1, i] != alignment.matrix[2, i]){                         #  p1 != p2   aka different resolutions in p1 and p2
          if(alignment.matrix[4, i] != alignment.matrix[3, i]){                       #  o != p3    durand says "less likely pattern due to seq. errors
            if(alignment.matrix[3, i] == alignment.matrix[1, i]) {baba <- baba + 1}   #  add to the count of baba sites
            if(alignment.matrix[2, i] == alignment.matrix[3, i]) {abba <- abba + 1}   #  add to the count of abba sites
          } 
        }
      }
    }
    d <- (abba - baba) / (abba + baba)   #what its all about    
  ## THIS SECTION WILL CALCULATE THE P-VAL BASED ON BOOTSTRAPPING
  ## SITES ARE SAMPLED WITH REPLACEMENT TO MAKE A NEW DATASET OF
  ## OF EQUAL SIZE TO THE ORIGINAL DATASET THIS ALLOWS US TO CALCULATE
  ## THE STANDARD DEVIATION AND THUS A Z SCORE.
    if(boot==T){
      sim.d<-vector()
      foo <- ncol(alignment.matrix)
      sim.matrix<-matrix(,4,foo)
      for(k in 1:replicate){      
        for(j in 1:4){
          sim.matrix[j,1:foo] <-sample(alignment.matrix[j,1:foo],replace=T)
        }
        t.abba <- t.baba <- 0                                                                         #  set up my variables
        for(i in 1:ncol(sim.matrix)){                                               #  run through all sites
          if(length(unique(sim.matrix[, i])) == 2){                                 #  unique(c(p1,p2,p3,o))==2 aka biallelic
            if(sim.matrix[1, i] != sim.matrix[2, i]){                               #  p1 != p2   aka different resolutions in p1 and p2
              if(sim.matrix[4, i] != sim.matrix[3, i]){                             #  o != p3    durand says "less likely pattern due to seq. errors
                if(sim.matrix[3, i] == sim.matrix[1, i]) {t.baba <- t.baba + 1}     #  add to the count of baba sites
                if(sim.matrix[2, i] == sim.matrix[3, i]) {t.abba <- t.abba + 1}     #  add to the count of abba sites
              } 
            }
          }
        }
        sim.d[k] <- (t.abba - t.baba) / (t.abba + t.baba)   #what its all about   
      }
      sd.sim.d <- round(sqrt(var(sim.d)),5)
      mn.sim.d <- round(mean(sim.d),5)
      new.pval <- 2*(pnorm(-abs(d/sd.sim.d)))
      ## NOW WE MAKE THE OUTPUTS  
      cat("\nSites in alignment =", ncol(alignment.matrix))
      cat("\nNumber of sites with ABBA pattern =", abba)
      cat("\nNumber of sites with BABA pattern =", baba)
      cat("\n\nD raw statistic / Z-score = ", d, " / ", d/sd.sim.d)
      cat("\n\nResults from ", replicate, "bootstraps")
      cat("\nSD D statistic =", sd.sim.d)
      cat("\nP-value (that D=0) = ",new.pval) #after Eaton and Ree 2013 
      results.matrix[q, 1:6] <- c(paste(starting,":",ending,sep=""), 
                                  abba, baba, d, d/sd.sim.d, new.pval)
      
    }
    if(boot==F){
      cat("\nSites in alignment =", ncol(alignment.matrix))
      cat("\nNumber of sites with ABBA pattern =", abba)
      cat("\nNumber of sites with BABA pattern =", baba)
      cat("\n\nD raw statistic = ", d)
      results.matrix[q, 1:4] <- c(paste(starting,":",ending,sep=""), 
                                  abba, baba, d)
    }
  }
  return(results.matrix)
}