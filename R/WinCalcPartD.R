WinCalcPartD <- function(alignment = "alignment.fasta", boot=F, replicate = 1000, 
                      win.size = 100, step.size=50, alpha =.05){
  alignment <- read.alignment(alignment, format = "fasta")                          #  read in the alignment
  alignment.matrix <- matrix(, length(alignment$nam), nchar(alignment$seq[[1]]))#  make a matrix for the alignment
  for(i in 1:length(alignment$nam)){
    alignment.matrix[i, ] <- unlist(strsplit(alignment$seq[[i]], ""))               #  fill in the matrix
  }  
  full.align <- alignment.matrix
  total <- ncol(full.align)
  spots <- seq(from = 1, to = (total - win.size), by = step.size)
  results.matrix <- as.data.frame(matrix(,1,16))
  colnames(results.matrix) <- c("range", "ABBAA", "BABAA", "ABABA", 
                                "BAABA", "ABABA", "BAABA", "d1", "d2", 
                                "d12","Z-d1", "Z-d2", "Z-d12", "d1-pval",  
                                "d2-pval",  "d12-pval")
  for(q in 1:length(spots)){  
    alignment.matrix <- full.align[,spots[q]:(spots[q] + win.size - 1)]
    starting <- spots[q]
    ending <- spots[q] + win.size - 1 
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
      results.matrix[q,1:16] <- c(paste(starting,":",ending,sep=""),
                                  abbaa, babaa, ababa, baaba, abbba,
                                  babba, d1, d2, d12, d1/sd.sim.d1, 
                                  d2/sd.sim.d2, d12/sd.sim.d12, 
                                  new.pval.d1, new.pval.d2, 
                                  new.pval.d12)
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
}