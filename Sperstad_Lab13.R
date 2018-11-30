#################################################
### Lab 13: Phylogenetic Comparative Methods  ###
#################################################
#Zachary Sperstad

###################### Exercise 1 ###################### 

require(geomorph)
require(RColorBrewer)
require(phytools)
tree<-read.tree('SpringerRAxML.tre')
tips<-tree$tip.label
master<-read.csv("master2018.csv")

nope<-subset(tips,tree$tip.label%in%master$PhyloSpp==F)

sub.tree<-drop.tip(tree,nope)
plot(sub.tree)
ladderized.tree<-ladderize(sub.tree)
plot(ladderized.tree,edge.width=2)
nodelabels(bg='white',frame='circle',cex=0.5)

#Q1- The root node is node #17.This node represents the most recent
#    common ancestor of the taxa included in this tree.

#Q2- The node that represents the common ancestor of Hominidae is node
#    #19.

sub.tips<-ladderized.tree$tip.label

###################### Exercise 2 ######################

quant<-master[,4:198]
p<-195/3
k<-3
array<-arrayspecs(quant,p,k)
gpa<-gpagen(array)
coords<-gpa$coords
aligned<-aperm(coords)
dim(aligned)<-c(547,195)
ID<-master[1:3]
ID.coord.matrix<-cbind(ID,aligned)


aggregated.data<-aggregate(ID.coord.matrix,list(ID.coord.matrix$PhyloSpp),
                           FUN=mean)
row.names(aggregated.data)<-aggregated.data$Group.1
aggregated.data$Fam<-c('Cercopethicidae','Hominidae','Hominidae','Hominidae',
                 'Hylobatidae','Hylobatidae','Hylobatidae','Hylobatidae',
                 'Hominidae','Hominidae','Cercopethicidae','Cercopethicidae',
                 'Cercopethicidae','Hominidae','Hominidae','Hylobatidae')
sig<-physignal(aggregated.data[5:199],ladderized.tree)
sig

#Q3-  I obtain a K value of 0.7496.
#Q4-  These data have less phylogenetic signal than is expected even under a
#     purely Brownian Motion process of evolution.

###################### Exercise 3 ######################

csize<-gpa$Csize
mean.csize<-aggregate(csize,list(ID.coord.matrix$PhyloSpp),FUN=mean)
mean.csize.extract<-mean.csize$x
names(mean.csize.extract)<-mean.csize$Group.1
mean.csize.extract
asr<-ace(mean.csize.extract,ladderized.tree,type='continuous')
contMap(ladderized.tree,mean.csize.extract,lwd=5)
nodelabels(bg='white',frame='circle',cex=0.5)
asr$ace
asr$CI95

#Q5- The estimated ancestral centroid size of the MCRA
#    of the catarrhine primates is ~37.5. The 95% CI for
#    this estimation is ~29.4-45.5.

mean.csize.extract

#Q6- It looks like Pan paniscus (~38.5) has a centroid size 
#    closest to the estimatd centroid size of the MRCA.

