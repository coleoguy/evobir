NextGen <- function(x, y, census, gametes, loci.a) {  # Assembles the next generation with selection   # x = gamete.pool    y = evaluation
  error.check <- 0
  while(error.check == 0){
    new.pop <- list()
    # SELECTION ACTS NOW BY EFFECTING THE PROBABILITY OF LEAVING OFFSPRING
    foo <- which(y[, 1] == "Male")
    foo2 <- which(y[, 1] == "Female")
    dads<-moms<-vector()
    for(i in 1:census){
      dads[i] <- sample(foo, 1, prob = y[foo, 2], replace = T) * 2 - sample(0:1, 1)   # select male gametes based on fathers fitness
      moms[i] <- sample(foo2, 1, prob = y[foo2, 2], replace = T) * 2 - sample(0:1, 1) # select female gametes based on mohters fitness
    }
    # NOW LETS ASSEMBLE OUR NEW GENOMES
    genome <- matrix(0,14,loci.a)
    rownames(genome) <- c("auto.1.s", "auto.1.d", "auto.2.s", "auto.2.d", "x.1.s",    "x.1.d",    "x.2.s",    "x.2.d", 
                          "y.1.s", "y.1.d", "xy.recom1", "xy.recom2", "haplosufficiency.x1","haplosufficiency.xy")
    for(i in 1:census){
      baz <- genome
      # and the rest of fertilization
      baz[c(3, 4, 5, 6, 11, 13), ] <- x[[moms[i]]][c(1, 2, 3, 4, 5, 6), ]
      baz[c(1, 2, 7, 8, 9, 10, 12, 14), ]  <- x[[dads[i]]][c(1, 2, 3, 4, 5, 6, 7, 8), ]
      new.pop[[i]] <- baz
    }
    qwert <- vector()
    for(k in 1:census){
      qwert[k] <- sex(new.pop[[k]])
    }
    if(length(unique(qwert))>1) error.check <- 1
  }
  return(new.pop)
}