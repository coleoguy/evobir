GetTaxonomy <- function(tree, database = "ncbi"){
  tree.names <- as.data.frame(tree$tip.label)                                   # lets get the names on the tree
  tree.names <- cbind(0, as.data.frame(gsub(pattern = "'", replacement = "", tree.names[,1])))  # some names in tree are quoted lets get rid of that
  for(i in 1:nrow(tree.names)){                                                                  # this loop will pull the genus names from the tree
    tree.names[i, 1] <- strsplit(as.character(tree.names[i, 2]), split = "_", fixed = T)[[1]][1]
  }
  colnames(tree.names) <- c("Genus", "species")  
  ## lets get the families for the species in the tree from NCBI
  famnames <- sapply(tree.names[,2], tax_name, get = "family", db = database)
  foo <- unlist(famnames)
  tree.names <- cbind(foo, tree.names)
  colnames(tree.names)[1] <- "Family"
  write.csv(tree.names, file = 'tree.taxonomy.csv', row.names=F)
}