ReplTable <- function(tree, data.taxonomy, tree.taxonomy, levels=c("Genus", "Family"), verbose=TRUE){
  match.higher.level <- function(tree.repl, tree, leftover.tips, level, data.taxonomy, tree.taxonomy, verbose){    
    counter <- 0
    for (j in 1:nrow(tree.repl)){
      if (tree.repl[j, "in.tree"] == "0"){
        # extract the clade we want to match to
        level.data <- data.taxonomy[data.taxonomy[,'Binomial'] == tree.repl[j, "in.data"], level]
        level.tree <- vector()
        for(k in 1:length(leftover.tips)){ ## this pulls out the possibly usable ones
          level.tree[k] <- tree.taxonomy[tree.taxonomy[,'Binomial'] == leftover.tips[k], level]
        }
        used.tree <- vector()
        used.names <- tree.repl[!tree.repl[, "in.tree"] %in% c("-", "0"), 'in.tree']
        for(k in 1:length(used.names)){
          used.tree[k] <- tree.taxonomy[tree.taxonomy[,'Binomial'] == used.names[k], level]
        }
        level.tree <- level.tree[!level.tree %in% used.tree]
        #now we need to find if level.data is in level.tree 
        if (length(which(level.data == level.tree)) >= 1) {
          tip <- tree.taxonomy[tree.taxonomy[, level] == level.data, 'Binomial']
          if (length(tip) > 1){
            tip <- sample(tip, 1)
          }
          leftover.tips <- leftover.tips[!leftover.tips == tip]
          tree.repl[j, "in.tree"] <- tip
          tree.repl[j, "match.level"] <- level
          counter <- counter + 1  
        }
      }
    }
    if (verbose == TRUE) cat("\n", (paste("There are", counter, "matches at the level of", level, sep=" ")))
    return(list(tree.repl=tree.repl, leftover.tips=leftover.tips))
  }
  ## create empty file -- this will be our eventual output to the user
  tree.repl <- cbind(data.taxonomy[, "Binomial"], rep("0", nrow(data.taxonomy)), rep("0", nrow(data.taxonomy)))
  colnames(tree.repl) <- c("in.data", "in.tree", "match.level")
  ## create a leftover.tips object to ensure that there are no replicates
  leftover.tips <- tree$tip.label  
  ## first lets get the perfect matches at a species level
  counter <- 0
  for (i in 1:nrow(tree.repl)){
    if (tree.repl[i, "in.data"] %in% tree$tip.label){
      temp <- which(tree$tip.label == tree.repl[i, "in.data"])
      if(length(temp)>1)
        temp <- temp[1]
      tip <- tree$tip.label[temp]
      leftover.tips <- leftover.tips[-temp] # make sure that these are not sampled
      tree.repl[i, "in.tree"] <- tip
      tree.repl[i, "match.level"] <- "Species"
      counter <- counter + 1
    }
  }
  if (verbose == TRUE) cat("\n", paste("There are", counter, "matches at the level of Species", sep=" "))
  for (i in 1:length(levels)){	
    level <- levels[i]
    matched <- tree.taxonomy[tree.taxonomy[,'Binomial'] %in% tree.repl[tree.repl[,2] != 0, 2], level] #reports either the genera or family that have been used
    for(j in 1:nrow(tree.repl)){ # goes through each line of replacemnt table
      if(tree.repl[j,'in.tree'] == 0){ # if the replacement table indicates no match
        if(data.taxonomy[data.taxonomy[, 'Binomial'] == tree.repl[j,'in.data'], level] %in% matched){
          tree.repl[j,'in.tree'] <- '-'
          if(i == 1)tree.repl[j,'match.level'] <- 'Genus'
          if(i == 2)tree.repl[j,'match.level'] <- 'Family'
        }
      }
    }
    ref <- match.higher.level(tree.repl, tree, leftover.tips, level, data.taxonomy, tree.taxonomy, verbose=T)
    tree.repl <- ref$tree.repl
    leftover.tips <- ref$leftover.tips
  } 
  repl.table <- tree.repl
  return(repl.table)
}
## Matching data and tree
MaxTips <- function(tree, data, repl.table){
  `%ni%` <- Negate(`%in%`)
  for(i in 1:length(tree$tip.label)){
    if(tree$tip.label[i] %in%  repl.table[, "in.tree"]){      # is the tip going to be used
      tree$tip.label[i] <- repl.table[which(repl.table[, "in.tree"] == tree$tip.label[i]), "in.data"]
    }
  }
  to.drop <- tree$tip.label[tree$tip.label %ni% repl.table[,"in.data"]]
  pruned.tree <- drop.tip(tree, to.drop)
  pruned.data <- data[data[, 'Binomial'] %in% pruned.tree$tip.label, ]
  count.matches <- repl.table[!repl.table[,'in.tree'] %in% c(0, '-'),]
  perfect.match <- sum(count.matches[,'match.level'] == 'Species')
  genus.level <- sum(count.matches[,'match.level'] == 'Genus')
  family.level <- sum(count.matches[,'match.level'] == 'Family')
  cat(paste("\nThe pruned tree has ", length(pruned.tree$tip.label), " tips.\n", sep = ""))
  cat(paste(perfect.match, "of these were species level matches\n"))
  cat(paste(genus.level, "of these were genus level matches\n"))
  cat(paste(family.level, "of these were family level matches\n"))
  return(list(pruned.tree=pruned.tree, pruned.data=pruned.data)) 
}