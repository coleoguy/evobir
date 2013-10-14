HighLevelTree <- function(taxa.table, tree, cur.tips, new.tips){
  ## some functions
  `%ni%` <- Negate(`%in%`)
  ## get one species per family
  new<-subset(taxa.table,!duplicated(taxa.table[,new.tips]))
  ## drop all but one rep per clade
  tree.pruned <- drop.tip(tree, tree$tip.label[tree$tip.label %ni% new[,cur.tips]])
  ## lame but swaps species and family names
  for(i in 1:length(tree.pruned$tip.label)){
    tree.pruned$tip.label[i]<-as.character(new[new[,cur.tips] == tree.pruned$tip.label[i],new.tips])
    }
  return(tree.pruned)
}
