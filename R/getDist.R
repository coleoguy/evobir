getDist <- function(map1, map2, method=NULL){
  condHist <- function(x, case.s){
    x <- cumsum(x)
    z <- x[1]
    for(i in 2:length(x)){
      if(names(z)[length(z)] == names(x)[i]){
        z[length(z)] <- x[i]
      }
      if(names(z)[length(z)] != names(x)[i]){
        z <- c(z, x[i])
      }
    }
    c.type <- vector()
    for(i in 1:(length(z)-1)){
      c.type <- c(c.type, paste(names(z)[i:(i+1)],collapse=""))
    }  
    if(case.s == "U"){
      c.type <- replace(c.type, c.type=="12", "F")       
      c.type <- replace(c.type, c.type=="21", "R")  
    }
    if(case.s == "L"){
      c.type <- replace(c.type, c.type=="12", "f")       
      c.type <- replace(c.type, c.type=="21", "r")  
    }
    z <- z[-length(z)]
    names(z) <- c.type
    return(z)
  }
  ###### END HELPER FUNCTIONS ######
  
  if(!identical(map1$edge, map2$edge)){
    stop("The phylogenies for the two mappings are not the same.  Only distances measured across the same phylogenies are possible.")
  }
  # getPaths is in evobiR
  pathways <- getPaths(map1, type="branch")
  results <- list()
  # I cant figure how to initialize this list
  # everything breaks this is a total hack but
  # it works
  results[[9]] <- 9999
  results <- results[-9]
  # to track the two traits capitol letters
  # refer to the first trait while lower case
  # refer to trait one "F" always indicates transition
  # from 1 to 2 "R" indicates transition from 
  # 2 to 1.
  names(results) <- c("Ff", "Fr", "Rf", "Rr",
                      "fF", "fR", "rF", "rR")
  char1 <- map1$maps
  char2 <- map2$maps
  states1 <- dimnames(map1$mapped.edge)[[2]]
  states2 <- dimnames(map2$mapped.edge)[[2]]
  for(i in 1:length(pathways)){
    branches <- pathways[[i]]
    hist1 <- unlist(char1[branches])
    hist2 <- unlist(char2[branches])
    # check to see if there are changes in both pathways
    if(length(unique(names(hist1))) > 1 &
       length(unique(names(hist2))) > 1){
      # combine paths to a single history
      c.hist <- sort(c(condHist(hist1, "U"), condHist(hist2, "L")))
      for(k in 1:(length(c.hist)-1)){
        # check to see if we are spanning a change in both traits
        if(sum(isupper(names(c.hist)[k:(k+1)])) ==
           sum(islower(names(c.hist)[k:(k+1)]))){
          c.type <- which(names(results) == paste(names(c.hist)[k:(k+1)], collapse=""))
          c.dist <- as.numeric(c.hist[k+1] - c.hist[k])
          results[[c.type]] <- c(results[[c.type]], c.dist)
        }
      }
    }
  }
  # changes that occur in more than one evol
  # pathway need to be cleared out but they
  # are always exactly identical so
  for(i in 1:length(results)){
    if(!is.null(results[[i]])) results[[i]] <- unique(results[[i]])
  }
  return(results)
}
