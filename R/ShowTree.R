ShowTree <- function(tree, 
                     tip.vals, 
                     pch = 16, 
                     cols = NULL, 
                     show.tip.label = F, 
                     tip.cex = 1, 
                     label.cex = 1, 
                     type = "phylogram"){
  if(is.null(cols)) cols <- viridis(length(unique(tip.vals)))
  .pardefault <- par(no.readonly = T)
  plot(tree, 
       show.tip.label = show.tip.label, 
       no.margin = T, 
       type = type, 
       cex = label.cex)
  
  tiplabels(pch = pch,
            col = cols[factor(tip.vals)], 
            frame = "none", 
            cex = tip.cex)
  par(.pardefault)
}