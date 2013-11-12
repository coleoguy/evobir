DesignLife <- function(loci.a) {    # this will create a genome file that our population will start with
  genome <- matrix(1, 14, loci.a)
  # naming convention for the rows is auto.1.s is autosome copy 1 fitness in sires
  # x1.d is the fitness of x chromosome copy 1 in a dam
  rownames(genome) <- c("auto.1.s", "auto.1.d", "auto.2.s", "auto.2.d", 
                        "x.1.s",    "x.1.d",    "x.2.s",    "x.2.d", 
                        "y.1.s", "y.1.d",                                                         # we have to have the y.1.d because recombination can turn thes into X alleles
                        "xy.recom1", "xy.recom2",                                                               # thinking about how recombination ceases this should be the union of the parents
                        "haplosufficiency.x1","haplosufficiency.xy")        # haplo-sufficiency should be inherited with an X or Y
  genome[13:14, ] <- 0                                                                            # we begin with no haplosufficiency
  return(genome)
}
