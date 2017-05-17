################################################
#                                              #
#  Heath Blackmon & Richard Adams              #
#  Continuous value at nodes producing a       #
#  derived state: August 10 2015               #
#                                              #
################################################
##
## INPUT DATA
## trees: a phylo or multiPhylo object
## data: a dataframe with three collumns (tips, cont trait, disc trait)
## derived.state: a text string or numeric matching one of the
##                entries in column 3 of data
## iterations the number of Monte Carlo simulations per tree
##                used to calc p-value

AncCond <- function(trees, data, iter = 1000) {
  ## create named vector for disc trait for all taxa
  dt.vec <- data[, 3]
  names(dt.vec) <- data[, 1]
  
  ## create named vector for cont trait taxa not in derived state
  ct.data <- data[data[, 3] != 1,]
  ct.vec <- as.numeric(ct.data[, 2])
  names(ct.vec) <- ct.data[, 1]
  
  ## ASR for the continuous trait
  anc.states.cont.trait <- anc.ML(trees, ct.vec, model = "BM")
  
  ## ASR for discrete trait
  ## using stochastic mappings to nail down specific transition points
  anc.state.dt <-
    make.simmap(
      trees,
      dt.vec,
      model = matrix(c(0, 0, 1, 0), 2),
      nsim = 1,
      pi = c(1, 0),
      message = F
    )
  
  ## Parse simmap to get producing nodes
  # the mapped edge object has time spent in a state in
  # two columns so only branches with a change have an entry
  # in both columns
  ss_nodes <- anc.state.dt$mapped.edge[, 1] > 0 &
    anc.state.dt$mapped.edge[, 2] > 0
  
  # this returns the node pairs describing a branch with origins
  wanted_branches <- ss_nodes[ss_nodes == T]
  wanted_nodes <- names(wanted_branches)
  
  # now we take the rootward node of each branch and get rid of duplicates
  wanted_nodes <- gsub(",.*", "", wanted_nodes)
  producing.nodes <- unique(wanted_nodes)
  
  ## get the mean ancestral value for the cont trait
  ## at nodes producing the derived state marginalizing across trees
  anc.states <- anc.states.cont.trait
  orig.val <- mean(anc.states$ace[names(anc.states$ace) %in%
                                    producing.nodes])
  orig.val <- mean(orig.val)
  
  ## Produce the null distribution of nodes in ancestral cond
  null.orig.val <- vector(length = iter)
  number.of.trans <- length(producing.nodes)
  anc.dt <- anc.state.dt
  anc.ct <- anc.states.cont.trait
  node.states <- describe.simmap(anc.dt)$states
  anc.cond.nodes <- anc.ct$ace[names(anc.ct$ace) %in%
                                 names(node.states)[node.states != 1]]
  for (j in 1:iter) {
    null.orig.val[j] <- mean(sample(anc.cond.nodes,
                                    length(producing.nodes)))
  }
  
  ## how many more extreme
  bigger <- (sum(null.orig.val >= orig.val) / iter)
  smaller <- (sum(null.orig.val <= orig.val) / iter)
  if (bigger < smaller)
    pval <- bigger
  if (smaller < bigger)
    pval <- smaller
  
  ## print results to terminal
  cat(paste(
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
  
  ## return results to user
  results <- list()
  results[[1]] <- orig.val
  results[[2]] <- number.of.trans
  results[[3]] <- null.orig.val
  results[[4]] <- pval
  names(results) <- c("OriginatingVals", "NTrans",
                      "NullDist", "pval")
  return(results)
}
