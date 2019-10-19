##################################################################
#                                                                #
#  Nathan W Anderson Heath Blackmon & Richard Adams              #
#  Continuous value at nodes producing a                         #
#  derived state: August 10 2015                                 #
#                                                                #
##################################################################
##
## INPUT DATA
## trees: a phylo object

## data: a dataframe with three columns (Labels, cont trait, disc trait)
## - Labels should match taxa labels in the phylogeny
## - continuous trait should be numeric values
## - discrete trait must be coded as 1 and 2 if one is ancestral then it must be coded as 1
## There should be no missing data. If a taxa does not have available cont
## and discrete data, please prune it from the tree

## mc: the number of Monte Carlo simulations per simulated dataset
## used to calc p-value

## drop.state: should be NULL unless working under the assumption 
## that one state is ancestral and the other derived and back 
## transitions are not possible. Using this assumption will
## ignore continuous data from taxa in the derived state

## mat: transition matrix. Should contain the rate 
## matrix for evolution of the discrete trait. Acceptable matrices are
## c(0,0,1,0), c(0,1,1,0), c(0,2,1,0)

## pi: The probabilities the root of the tree are either of the
## discrete character states same values possible as make.simmap:
## "equal", "estimated", or vector length 2 with probabilities 
## for each state

## n.tails: either 1 or 2 depending on whether user has apriori hypothesis about a certain state


AncCond <- function(tree, data, mc = 1000, drop.state=NULL, mat=c(0,2,1,0), pi="equal", n.tails = 1, message = T) {
  ##### testing inputs #####
  if(class(tree) != 'phylo') {stop('tree must be class phylo')}
  if(!is.data.frame(data) & ncol(data) == 3){stop('data should be a dataframe with 3 columns\n(tip labels, cont data, discrete data)')}
  if(class(mc) != 'numeric' | round(mc) != mc | mc < 1){stop('mc should be a numeric positive integer integer')}
  if(!is.null(drop.state)) if(!drop.state %in% c(1,2)){stop('drop.state must be NULL, or numeric 1 or 2')}
  if(!sum(mat == c(0,0,1,0)) == 4 & !sum(mat == c(0,1,1,0)) == 4 & !sum(mat == c(0,2,1,0)) == 4){
    stop('mat must be a vector of the form c(0,0,1,0), c(0,1,1,0), or c(0,2,1,0)')
  }
  
  if((!pi %in% c('equal', 'estimated'))[1]){
    if(!is.numeric(pi)) stop('pi must be equal, estimated or a vector of length 2\nwith probabilities for the state of the discrete character at the root')
    if(length(pi) != 2 | sum(pi) != 1) stop('pi must be equal, estimated or a vector of length 2\nwith probabilities for the state of the discrete character at the root')
  }
  
  if(n.tails != 1 & n.tails != 2){stop('n.tails should be numeric 1 or 2')}
  #####
  #####

  ##### create named vector for disc trait for all taxa #####
  dt.vec <- data[, 3]
  names(dt.vec) <- data[, 1]
  
  ##### create named vector for cont trait taxa not in derived state #####
  if(!is.null(drop.state)){
    ct.data <- data[(data[, 3] != drop.state),]
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
  }else{
    ct.data <- data
    ct.vec <- as.numeric(ct.data[, 2])
    names(ct.vec) <- ct.data[, 1]
  }
  if(sum(is.na(ct.vec)) > 0 | sum(is.na(dt.vec)) > 0){
    stop('There exists missing trait data for some species in the phylogeny.\n
         Please remove such taxa from the tree.')
  }
  ## ASR for the continuous trait
  anc.states.cont.trait <- anc.ML(tree, ct.vec, model = "BM")
  
  ## ASR for discrete trait
  ## using stochastic mappings to nail down specific transition points
  anc.state.dt <- make.simmap(tree, dt.vec,
                              model = matrix(mat, 2),
                              nsim = 1,
                              pi = pi,
                              message = F)
  
  ## Parse simmap to get producing nodes
  # the mapped edge object has time spent in a state in
  # two columns so only branches with a change have an entry
  # in both columns
  ss_nodes <- anc.state.dt$mapped.edge[, 1] > 0 &
    anc.state.dt$mapped.edge[, 2] > 0
  
  # this returns the node pairs describing a branch with origins
  wanted_branches <- ss_nodes[ss_nodes == T]
  wanted_nodes <- names(wanted_branches)
  
  if(sum(mat) > 1){
    # for the general model we partition the producing nodes for 1->2 and 1<-2 transitions
    producing.nodes12 <- c()
    producing.nodes21 <-c()
    trans.maps <- anc.state.dt$maps[ss_nodes == T]
    # now we take the rootward node of each branch and get rid of duplicates
    wanted_nodes <- gsub(",.*", "", wanted_nodes)
    ##### Just realized we can do this with describe.simmap :( 
    ##### But i dont want to change it, it would require match function
    for(i in 1:length(wanted_nodes)){
      if(names(trans.maps[[i]])[1] == '1'){
        producing.nodes12 <- c(producing.nodes12, wanted_nodes[i])
      }else if(names(trans.maps[[i]])[1] == '2'){
        producing.nodes21 <- c(producing.nodes21, wanted_nodes[i])
      }
    }
    producing.nodes12 <- unique(producing.nodes12)
    producing.nodes21 <- unique(producing.nodes21)
    ## get the mean ancestral value for the cont trait
    ## at nodes producing the derived state marginalizing across trees
    anc.states <- anc.states.cont.trait
    orig.val12 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                        producing.nodes12])
    orig.val21 <- mean(anc.states$ace[names(anc.states$ace) %in%
                                        producing.nodes21])
    ## Produce the null distribution of nodes in ancestral cond
    null.orig.val12 <- vector(length = mc)
    null.orig.val21 <- vector(length = mc)
    number.of.trans12 <- length(producing.nodes12)
    number.of.trans21 <- length(producing.nodes21)
    anc.dt <- anc.state.dt
    anc.ct <- anc.states.cont.trait
    node.states <- describe.simmap(anc.dt)$states
    anc.cond.nodes12 <- anc.ct$ace[names(anc.ct$ace) %in%
                                     names(node.states)[node.states != '2']]
    anc.cond.nodes21 <- anc.ct$ace[names(anc.ct$ace) %in%
                                     names(node.states)[node.states != '1']]
    
    for (j in 1:mc){
      # set.seed(j)
      null.orig.val12[j] <- mean(sample(anc.cond.nodes12,
                                        length(producing.nodes12)))
      null.orig.val21[j] <- mean(sample(anc.cond.nodes21,
                                        length(producing.nodes21)))
    }
    ## how many more extreme
    
    bigger12 <- (sum(null.orig.val12 >= orig.val12) / mc)
    smaller12 <- (sum(null.orig.val12 <= orig.val12) / mc)
    if(!is.null(producing.nodes12)){ 
      if (bigger12 <= smaller12){pval12 <- bigger12}
      if (smaller12 < bigger12){pval12 <- smaller12}
      if (n.tails == 2){pval12 <- 2 * pval12}
    }else{
      pval12 <- NA
    }
    
    bigger21 <- (sum(null.orig.val21 >= orig.val21) / mc)
    smaller21 <- (sum(null.orig.val21 <= orig.val21) / mc)
    if(!is.null(producing.nodes21)){ 
      if (bigger21 <= smaller21){pval21 <- bigger21}
      if (smaller21 < bigger21){pval21 <- smaller21}
      if (n.tails == 2){pval21 <- 2 * pval21}
    }else{
      pval21 <- NA
    }
    
    ## print results to terminal
    if (message == T){
      cat(paste(
        "Mean value for the continuous trait at origin oftrait 2:",
        round(orig.val12, digits = 4),
        "\n"
      ))
      cat(paste(
        "Mean value for the continuous trait at origin of trait 1:",
        round(orig.val21, digits = 4),
        "\n"
      ))
      cat(paste("Number of producing nodes 1->2:", round(mean(number.of.trans12), 
                                                    digits = 4), "\n"))
      cat(paste("Number of producing nodes 2->1:", round(mean(number.of.trans21), 
                                                         digits = 4), "\n"))
      cat(paste("Mean of null dist 1->2:", round(mean(null.orig.val12), 
                                            digits = 4), "\n"))
      cat(paste("Mean of null dist 2->1:", round(mean(null.orig.val21), 
                                                 digits = 4), "\n"))
      cat(paste("SD of null dist 1->2:", round(sd(null.orig.val12), digits = 4), "\n"))
      cat(paste("SD of null dist 2->1:", round(sd(null.orig.val21), digits = 4), "\n"))
      
      cat(paste("pvalue 1->2:", round(pval12, digits = 4), "\n"))
      cat(paste("pvalue 2->1:", round(pval21, digits = 4), "\n\n"))
      if(is.null(producing.nodes12)){cat('No 1 -> 2 transitions occured NA and NaN values produced.')}
      if(is.null(producing.nodes21)){cat('No 2 -> 1 transitions occured NA and NaN values produced.')}
    }
    
    ## return results to user
    results <- list()
    results[[1]] <- orig.val12
    results[[2]] <- number.of.trans12
    results[[3]] <- null.orig.val12
    results[[4]] <- pval12
    results[[5]] <- orig.val21
    results[[6]] <- number.of.trans21
    results[[7]] <- null.orig.val21
    results[[8]] <- pval21
    names(results) <- c("OriginatingVals1->2", "NTrans1->2",
                        "NullDist1->2", "pval1->2","OriginatingVals2->1", "NTrans2->1",
                        "NullDist2->1", "pval2->1")
  }else{
    # now we take the rootward node of each branch and get rid of duplicates
    wanted_nodes <- gsub(",.*", "", wanted_nodes)
    producing.nodes <- unique(wanted_nodes)
    ## get the mean ancestral value for the cont trait
    ## at nodes producing the derived state marginalizing across trees
    anc.states <- anc.states.cont.trait
    orig.val <- mean(anc.states$ace[names(anc.states$ace) %in%
                                      producing.nodes])
    
    ## Produce the null distribution of nodes in ancestral cond
    null.orig.val <- vector(length = mc)
    number.of.trans <- length(producing.nodes)
    anc.dt <- anc.state.dt
    anc.ct <- anc.states.cont.trait
    node.states <- describe.simmap(anc.dt)$states
    anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in%
                                   names(node.states)[node.states != '2']]
    
    for (j in 1:mc){
      # set.seed(j)
      null.orig.val[j] <- mean(sample(anc.cond.nodes,
                                      length(producing.nodes)))
    }
    ## how many more extreme
    
    bigger <- (sum(null.orig.val >= orig.val) / mc)
    smaller <- (sum(null.orig.val <= orig.val) / mc)
    if (bigger <= smaller){pval <- bigger}
    if (smaller < bigger){pval <- smaller}
    if (n.tails == 2){
      pval <- 2 * pval
    }
    ## print results to terminal
    if(message == T){cat(paste(
      "Mean value for the continuous trait at origin of derived trait:",
      round(orig.val, digits = 4),
      "\n"
    ))
      cat(paste("Number of producing nodes:", round(mean(number.of.trans), 
                                                    digits = 4), "\n"))
      cat(paste("Mean of null dist:", round(mean(null.orig.val), 
                                            digits = 4), "\n"))
      cat(paste("SD of null dist:", round(sd(null.orig.val), digits = 4), "\n"))
      
      cat(paste("pvalue:", round(pval, digits = 4), "\n\n"))
    }
    
    ## return results to user
    results <- list()
    results[[1]] <- orig.val
    results[[2]] <- number.of.trans
    results[[3]] <- null.orig.val
    results[[4]] <- pval
    names(results) <- c("OriginatingVals", "NTrans",
                        "NullDist", "pval")
  }
  return(results)
}



























