ShowTree <- function(tree, tip.vals, pch = 16, cols = NULL, tip.cex=1, ...){
  if(is.null(cols)) cols <- viridis(length(unique(tip.vals)))
  plot(tree, show.tip.label = F, no.margin=T)
  # make a color vector for the tips
  tiplabels(pch=pch,col=cols[factor(tip.vals)], frame="none", cex=tip.cex)
}