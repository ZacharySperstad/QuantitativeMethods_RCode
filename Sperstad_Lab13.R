#################################################
### Lab 13: Phylogenetic Comparative Methods  ###
#################################################
#Zachary Sperstad

require(geomorph)
require(RColorBrewer)
require(phytools)
tree<-read.tree('SpringerRAxML.tre')
tips<-tree$tip.label

master2018<-read.csv("master2018.csv")
nope<-subset(tips, )