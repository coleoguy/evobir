getData <- function(tree, dat, id.col, tr.col, method){
  if(!method %in% c("R", "Me", "Mo")) stop("must specify a valid method")
  # first we create a vector to hold our new data
  # and give the indices names that match our
  # tree tip labels
  tip.matches <- tree$tip.label[tree$tip.label %in% dat[,id.col]]
  tip.vals <- vector(length=length(tip.matches))
  names(tip.vals) <- tip.matches
  
  # lets make this thing handle multiple traits too  
  v.or.m <- "v"
  if(length(tr.col > 1)){
    v.or.m <- "m"  
    tip.mat <- matrix(,length(tip.vals), length(tr.col))
    row.names(tip.mat) <- tip.matches
    colnames(tip.mat) <- colnames(dat)[tr.col]
  }
  
    if(method=="R"){
      for(i in 1:length(tip.vals)){
        x <- sample(which(tip.matches[i] == dat[,id.col]), 1)
        if(v.or.m == "v") tip.vals[i] <- as.numeric(dat[x, tr.col])
        if(v.or.m == "m") tip.mat[i, 1:length(tr.col)] <- as.numeric(dat[x, tr.col])
      }
    }
    if(method=="Me"){
      for(i in 1:length(tip.vals)){
        x <- which(tip.matches[i] == dat[,id.col])
        if(v.or.m == "v") tip.vals[i] <- mean(as.numeric(dat[x, tr.col]))
        if(v.or.m == "m"){
          for(j in 1:length(tr.col)){
            tip.mat[i, j] <- mean(as.numeric(dat[x, tr.col[j]]))
          }
        }
      }
    }
    if(method == "Mo"){
      for(i in 1:length(tip.vals)){
        x <- which(tip.matches[i] == dat[,id.col])
        if(v.or.m == "v") tip.vals[i] <- as.numeric(Mode(dat[x, tr.col]))
        if(v.or.m == "m"){
          for(j in 1:length(tr.col)){
            tip.mat[i, j] <- Mode(as.numeric(dat[x, tr.col[j]]))
          }
        }
      }
    }
  drop <- tree$tip.label[!tree$tip.label %in% tip.matches]
  new.tree <- drop.tip(phy=tree, tip=drop)
  if(v.or.m == "m") return(list(new.tree, tip.mat))
  if(v.or.m == "v") return(list(new.tree, tip.vals))
}

