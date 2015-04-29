# evobiR test

library(devtools)
install_github('coleoguy/evobir')
library(evobiR)

CalcD(alignment = system.file("1.fasta", package = "evobiR"), boot = TRUE, replicate=10)

CalcPartD(alignment = system.file("2.fasta", package = "evobiR"), boot = TRUE, 10)

CalcPopD(alignment = system.file("3.fasta", package = "evobiR"))

data(hym.tree)
names <- c("Pepsis_elegans", "Plagiolepis_alluaudi", "Pheidele_lucreti",
           "Meliturgula_scriptifronsi", "Andrena_afimbriat")
FuzzyMatch(tree = hym.tree, data = names, max.dist=3)
AICc(-32, 3, 100)
data(data.mite)
data(trees.mite)
AncCond(trees[[1]], data, derived.state = "haplodiploidy", iterations=10) 

Even(c(1,2,3,4,5,6,2,5))
Mode(c(1,2,3,4,5,6,2,5))
data(trees)
data(mcmc2)
data(mcmc3)
# 1 tree 100 q-mats 3 states
PPSDiscrete(trees[[1]], MCMC=mcmc3[,2:10], states=c(.5,.2,.3), N=2)
data(finch)
phy <- finch$phy
data <- finch$data[,"wingL"]
# we create a second dataset with the order of taxa randomized
data <- sample(data, size=13)
names(data) == phy$tip.label ## order is different
data <- ReorderData(phy, data)
names(data) == phy$tip.label ## order is the same
data <- read.csv(file = system.file("horn.beetle.csv", package = "evobiR"))
ResSel(data = data, traits = c(2,3), percent = 15, identifier = 1, model = "linear")
SampleTrees(trees = system.file("trees.nex", package = "evobiR"), 
            burnin = .1, final.number = 20, format = 'new', prefix = 'sample')
data <- c(1,2,1,2,10,2,1,2,1,2,3,4,5,6,2,5)
SlidingWindow("mean", data, 3, 1)
SuperMatrix(missing = "N", prefix = "DATASET2", save = T)
WinCalcD(alignment = system.file("1.fasta", package = "evobiR"), win.size=100, step.size=50, boot = TRUE, replicate=10)

WinCalcPartD(alignment = system.file("2.fasta", package = "evobiR"), 
             win.size=1000, step.size=500, boot = F, replicate=100) 



