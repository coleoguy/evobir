cdcModel <- function(x, y, tree, model, initial.vals=NULL){
  # this is the constraint for a model where trait 2 changes at different speed right
  # after transition in character 1
  make343 <- function(tree, data){
    lik <- make.musse(tree, data, k=6, strict=F)
    lik.c <- constrain(lik, q13 ~ q24, q35 ~ q46, q12 ~ q56,
                       q14 ~ 0, q15 ~ 0, q16 ~ 0, q21 ~ 0, 
                       q23 ~ 0, q25 ~ 0, q26 ~ 0, q31 ~ 0, 
                       q32 ~ 0, q36 ~ 0, q41 ~ 0, q42 ~ 0, 
                       q43 ~ 0, q45 ~ 0, q51 ~ 0, q52 ~ 0, 
                       q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, 
                       q63 ~ 0, q64 ~ 0, q65 ~ 0, 
                       lambda6 ~ 1, lambda5 ~ 1, lambda4 ~ 1,
                       lambda3 ~ 1, lambda2 ~ 1, lambda1 ~ 1,
                       mu6 ~ .4, mu5 ~ .4, mu4 ~ .4, mu3 ~ .4,
                       mu2 ~ .4, mu1 ~ .4)
    return(lik.c)
  }
  # this is the constraint for a model where all transitions in trait 2 are the same
  make333 <- function(tree, data){
    lik <- make.musse(tree, data, k=6, strict=F)
    lik.c <- constrain(lik, q13 ~ q24, q35 ~ q46, q34~q56, q12 ~ q56,
                       q14 ~ 0, q15 ~ 0, q16 ~ 0, q21 ~ 0, 
                       q23 ~ 0, q25 ~ 0, q26 ~ 0, q31 ~ 0, 
                       q32 ~ 0, q36 ~ 0, q41 ~ 0, q42 ~ 0, 
                       q43 ~ 0, q45 ~ 0, q51 ~ 0, q52 ~ 0, 
                       q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, 
                       q63 ~ 0, q64 ~ 0, q65 ~ 0, 
                       lambda6 ~ 1, lambda5 ~ 1, lambda4 ~ 1,
                       lambda3 ~ 1, lambda2 ~ 1, lambda1 ~ 1,
                       mu6 ~ .4, mu5 ~ .4, mu4 ~ .4, mu3 ~ .4,
                       mu2 ~ .4, mu1 ~ .4)
    return(lik.c)
  }
  # this is the constraint for a model where all transitions sensu pagel dependent
  make344 <- function(tree, data){
    lik <- make.musse(tree, data, k=6, strict=F)
    lik.c <- constrain(lik, q13 ~ q24, q35 ~ q46, q34~q56,
                       q14 ~ 0, q15 ~ 0, q16 ~ 0, q21 ~ 0, 
                       q23 ~ 0, q25 ~ 0, q26 ~ 0, q31 ~ 0, 
                       q32 ~ 0, q36 ~ 0, q41 ~ 0, q42 ~ 0, 
                       q43 ~ 0, q45 ~ 0, q51 ~ 0, q52 ~ 0, 
                       q53 ~ 0, q54 ~ 0, q61 ~ 0, q62 ~ 0, 
                       q63 ~ 0, q64 ~ 0, q65 ~ 0, 
                       lambda6 ~ 1, lambda5 ~ 1, lambda4 ~ 1,
                       lambda3 ~ 1, lambda2 ~ 1, lambda1 ~ 1,
                       mu6 ~ .4, mu5 ~ .4, mu4 ~ .4, mu3 ~ .4,
                       mu2 ~ .4, mu1 ~ .4)
    return(lik.c)
  }
  # this function take a tree and two binary trait vectors with data in the same order as the trees
  # and returns a matrix recoded into a 6 state character with row names matching the tree
  # state 1 of character 1 is considered to have a hidden state and so is parsed into having 50%
  # prob of observing either state
  dataPrep <- function(data1, data2, tree){
    x <- apply(cbind(data1, data2), 1, paste, collapse="")
    convertBtoHex <- function(x){
      y <- vector()
      if(x=="00") y <- 1
      if(x=="01") y <- 2
      if(x=="10") y <- 3
      if(x=="11") y <- 4
      return(y)
    }
    x <- vapply(x, convertBtoHex, FUN.VALUE=1)
    tip.mat <- matrix(0, length(x), 6)
    row.names(tip.mat) <- tree$tip.label
    colnames(tip.mat) <- c("AB", "Ab", "a1B", "a1b", "a2B", "a2b")
    for(i in 1:nrow(tip.mat)){
      if(x[i]==1) tip.mat[i,1]<- 1
      if(x[i]==2) tip.mat[i,2]<- 1
      if(x[i]==3) tip.mat[i,c(3,5)]<- .5
      if(x[i]==4) tip.mat[i,c(4,6)]<- .5
    }
    return(tip.mat)
  }

  ## and now the actual code to perform the analysis
  tip.mat <- dataPrep(x, y, tree)
  if(model == 333) lik <- make333(tree, tip.mat)
  if(model == 344) lik <- make344(tree, tip.mat)
  if(model == 343) lik <- make343(tree, tip.mat)
  if(is.null(initial.vals)) initial.vals <- rep(.1, length(argnames(lik)))
  result <- find.mle(lik, x.init=initial.vals)
  return(result)
}
