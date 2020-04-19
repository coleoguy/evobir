getNe <- function(males, females, locus="A"){
  if(locus == "A"){
    ne <- (4 * males * females) / (males + females)
  }
  if(locus == "X"){
    ne <- (9 * males * females) / (4*males + 2*females)
  }
  return(ne)
}
