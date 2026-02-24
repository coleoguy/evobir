## ----echo=F, tidy=T-----------------------------------------------------------
library(evobiR)


## -----------------------------------------------------------------------------
data(hym.tree)
names <- c("Pepsis_elegans", "Plagiolepis_alluaudi", "Pheidele_lucreti", "Meliturgula_scriptifronsi", "Andrena_afimbriat")
FuzzyMatch(tree = hym.tree, data = names, max.dist=3)


## ----eval=F-------------------------------------------------------------------
# SampleTrees(trees = system.file("trees.nex", package = "evobiR"),
#             burnin = .1, final.number = 20, format = 'new', prefix = 'sample')


## -----------------------------------------------------------------------------
CalcD(alignment = system.file("1.fasta", package = "evobiR"), sig.test = "N")

CalcPopD(alignment = system.file("3.fasta", package = "evobiR"), sig.test = "N")



## ----fig.height=5, fig.width=5------------------------------------------------
data(sunspot.month)
foo <- as.vector(sunspot.month)
plot(foo, ylab="Number of sunspots", xlab="record number")


## ----fig.height=5, fig.width=5------------------------------------------------
# first lets use the sliding window function to get the number sun spots

# average over 11 years and do this moving in 1 year steps through time
sunspots <- SlidingWindow(FUN = "mean",
                          data = foo,
                          window = 132,
                          step = 12,
                          strict = F)

# build a year vector matching the data length dynamically
start.yr <- start(sunspot.month)[1]
end.yr <- start.yr + length(foo) / 12
year.vec <- seq(start.yr, end.yr, length.out = length(foo))

# we repeat this on the years so we can plot against a sensible x axis
years <- round(SlidingWindow("mean",
                             data = year.vec,
                             window = 132,
                             step = 12,
                             strict = F))
{plot(x = years, y = sunspots, type = "l", lwd = 3)
abline(v = 1810, col = "red", lwd = 3)
text(y = 90, x = 1810, "Dalton Min.", col = "red", pos = 4)}


## -----------------------------------------------------------------------------
data <- read.csv(file = system.file("horn.beetle.csv", package = "evobiR"))


## -----------------------------------------------------------------------------
data[1:10,]


## ----fig.height=5, fig.width=5------------------------------------------------
ResSel(data = data, traits = c(2,3), percent = 15, identifier = 1, model = "linear")

