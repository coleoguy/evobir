ResSel <- function(data, traits, percent = 10, identifier=1, model = "linear"){
  ## now lets fit a linear model
  lm.norm <- lm(data[,traits[2]] ~ data[,traits[1]])
  foo<-rank(lm.norm$residuals)
  ## lets plot the data just for grins
  plot(data[,traits[1]:traits[2]], pch=19, cex=.5, main=paste("Top and Bottom ",percent,"%",sep=""), col="gray")
  points(data[foo>(nrow(data)-(nrow(data)*.01*percent)),traits[1]:traits[2]],col="blue",pch=19,cex=1.1)
  points(data[foo<(nrow(data)*.01*percent),traits[1]:traits[2]],col="red",pch=19,cex=1.1)
  abline(lm.norm, col='orange', lwd=5)
  high <- data[foo>(nrow(data)-(nrow(data)*.01*percent)),1]
  low <- data[foo<(nrow(data)*.01*percent),1]
  results <- list()
  results[[1]] <- high
  results[[2]] <- low
  names(results) <- c('high line', 'low line')
  return(results)
}