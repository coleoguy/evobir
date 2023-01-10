ReOrderAlignment <- function(file, newfile, ref=NULL){
  dat <- read.fasta(file)
  # validate inputs
  if(!is.null(ref)){
    if(is.numeric(ref[1])){
      ref <- names(dat)[ref]
    }
    if(!sum(ref %in% names(dat)) == length(ref)){
      stop("One of your reference sequences was not found in your file")
    }
  }
  # find correct order
  startspots <- c()
  for(i in 1:length(dat)){
    if(dat[[i]][1] != "-"){
      startspots[i] <- 1
    }else{
      bar <- which(dat[[i]] == "-")
      vec <- bar[-1] - bar[-length(bar)]
      if(length(vec) == sum(vec)){
        startspots[i] <- length(vec)
      }else{
        startspots[i] <- min(which(vec != 1))
      }
    }
  }
  seqorder <- order(startspots)
  starts <- which(names(dat) %in% ref)
  seqorder <- c(starts, seqorder[!seqorder %in% starts])
  # make new seq list
  newlist <- dat[seqorder]
  write.fasta(newlist, names = names(newlist), file.out = newfile)
  print("Done!")
}