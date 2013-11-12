Monitor <- function(population, report.style, report.freq, k, loci.s, generations, evaluation, counter, rep.cols){
  if(report.style == "population"){
    return(population)
  }
  if(report.style == "y-chrom"){
    ychrom.s <-  ychrom.d <- matrix(,length(population), loci.s)
    for(i in 1:length(population)){
      ychrom.s[i,] <- population[[i]][9, 1:loci.s]
      ychrom.d[i,] <- population[[i]][10, 1:loci.s]
    }
    data <- list(ychrom.s, ychrom.d)
    return(data)
  }
  if(report.style =="PAR-active"){
    if(counter == 1){
      plot(0,0, col = "white", main = "Proportion of sex chromosome loci in\n the population that are restricted to the Y",
           xlab = "sites", ylab = "proportion", xlim = c(1, loci.s), 
           ylim = c(0,1))
      Sys.sleep(.2)
    }
    males <- which(evaluation[,1]=="Male")
    foob <- matrix(, length(males), loci.s)
    for(i in 1:length(males)){
      foob[i, ] <- population[[males[i]]][12,1:loci.s]
    }
    baq <- colMeans(foob)
    lines(baq, col=rep.cols[counter])
    Sys.sleep(.2)
  }
  if(report.style =="PAR-save"){
    if(counter == 1) dir.create("results")
    if(counter <10) fix.name <- paste("000", counter, ".png", sep="")
    if((counter < 100) & (counter > 9)) fix.name <- paste("00", counter, ".png", sep="")
    if(counter > 99) fix.name <- paste("0", counter, ".png", sep="")
    png(file= paste("results/", fix.name, sep=""))
    plot(0,0, col = "white", main = "Proportion of sex chromosome loci in\n the population that are restricted to the Y",
           xlab = "sites", ylab = "proportion", xlim = c(1, loci.s), 
           ylim = c(0,1))
    males <- which(evaluation[,1]=="Male")
    foob <- matrix(, length(males), loci.s)
    for(i in 1:length(males)){
      foob[i, ] <- population[[males[i]]][12,1:loci.s]
    }
    baq <- colMeans(foob)
    lines(baq, col=rep.cols[counter], lwd=3)
    dev.off()
  }
  
  
  
  
}