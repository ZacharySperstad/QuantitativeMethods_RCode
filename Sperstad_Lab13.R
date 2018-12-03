#################################################
### Lab 13: Phylogenetic Comparative Methods  ###
#################################################
#Zachary Sperstad

###################### Exercise 1 ###################### 

require(geomorph) #Allows geomorph functions to be performed.
require(RColorBrewer) #Allows RColorBrewer functions to be performed.
require(phytools) #Allows phytools functions to be performed.
tree<-read.tree('SpringerRAxML.tre') #Reads in tree.
tips<-tree$tip.label #Extracts tip labels as an object.
master<-read.csv("master2018.csv") #Reads in master2018.csv file.

nope<-subset(tips,tree$tip.label%in%master$PhyloSpp==F) #Removes taxa in tree that are not
                                                        #present in the CSV file.
sub.tree<-drop.tip(tree,nope) #Creates a subtree without taxa in the "nope" object.
plot(sub.tree) #Plots the subtree.
ladderized.tree<-ladderize(sub.tree) #Ladderizes the subtree.
plot(ladderized.tree,edge.width=2) #Plots the labberized subtree with thicker
                                   #branches.
nodelabels(bg='white',frame='circle',cex=0.5) #Adds white, circular node
                                              #labels to the tree.

#Q1- The root node is node #17.This node represents the most recent
#    common ancestor of the taxa included in this tree.

#Q2- The node that represents the common ancestor of Hominidae is node
#    #19.

sub.tips<-ladderized.tree$tip.label #Creates an object with node labels
                                    #from our subtree.

###################### Exercise 2 ######################

quant<-master[,4:198] #Extracts the quantitative data from the  original
                      #dataset.
p<-195/3 #Denotes how many 3D landmarks we will have.
k<-3 #Denotes the number of dimensions in which our landmarks are.
array<-arrayspecs(quant,p,k) #Creates an array out of the 3D landmarks of
                             #each individual.
gpa<-gpagen(array) #Performs a GPA on the array.
coords<-gpa$coords #Extracts the GPA landmarks from the GPA object.
aligned<-aperm(coords) #Allows us to adjust the dimensions of the matrix
                       #containing our landmarks.
dim(aligned)<-c(547,195) #Changes our landmark matrix dimensions to 547x195.
ID<-master[1:3] #Creates an object containing each individual's ID, their binomial
                #the name of the family to which they belong.
ID.coord.matrix<-cbind(ID,aligned)

aggregated.data<-aggregate(ID.coord.matrix,list(ID.coord.matrix$PhyloSpp),
                           FUN=mean) #Aggregates coordinates by species identity.
row.names(aggregated.data)<-aggregated.data$Group.1 #Assigns Group.1 names to aggregated.data. object.
aggregated.data$Fam<-c('Cercopethicidae','Hominidae','Hominidae','Hominidae',
                       'Hylobatidae','Hylobatidae','Hylobatidae','Hylobatidae',
                       'Hominidae','Hominidae','Cercopethicidae','Cercopethicidae',
                       'Cercopethicidae','Hominidae','Hominidae','Hylobatidae') #Applies family name to each row.
sig<-physignal(aggregated.data[5:199],ladderized.tree)
sig

#Q3-  I obtain a K value of 0.7496.
#Q4-  These data have less phylogenetic signal than is expected even under a
#     purely Brownian Motion process of evolution.

###################### Exercise 3 ######################

csize<-gpa$Csize #Extracts centroid size from GPA object
mean.csize<-aggregate(csize,list(ID.coord.matrix$PhyloSpp),FUN=mean) #Finds mean size value 
                                                                     #for each species.
mean.csize.extract<-mean.csize$x #Removes size from mean.csize object.
names(mean.csize.extract)<-mean.csize$Group.1 #Assigns names of Group.1 to
                                              #mean.csize extract object.
mean.csize.extract #Shows contents of mean.csize.extract.
asr<-ace(mean.csize.extract,ladderized.tree,type='continuous')
contMap(ladderized.tree,mean.csize.extract,lwd=5) #Plots a tree with the
                                                  #reconstruction of "size"
                                                  #through the phylogeny...RIGHTEOUS!
nodelabels(bg='white',frame='circle',cex=0.5) #Adds white, circular node
                                              #labels to my pretty tree.
asr$ace
asr$CI95

#Q5- The estimated ancestral centroid size of the MCRA
#    of the catarrhine primates is ~37.5. The 95% CI for
#    this estimation is ~29.4-45.5.

mean.csize.extract

#Q6- It looks like Pan paniscus (~38.5) has a centroid size 
#    closest to the estimatd centroid size of the MRCA.

mean.quant<-aggregated.data[5:199] #Creats a subset of the data with only
                                   #quantitative data.
mean.array<-arrayspecs(mean.quant,p,k) #Creates an array from the (averaged by species) quantitative data.
#dim(mean.array)<-c(16,195) #Changes dimension of arrays to 16x195.
#rownames(mean.array)<-dimnames(aggregated.data$Group1)[[3]] #Assigns names in "Group1" as 
dimnames(mean.array)[[3]]<-rownames(aggregated.data) #row names of the mean.array object.

fam.col<-c("green1","blue1",
           "purple1")[as.factor(aggregated.data[sub.tips,4])] #Assigns specific colors to each family.


pms<-plotGMPhyloMorphoSpace(ladderized.tree,mean.array,ancState=T,
                            plot.param=list(t.bg=fam.col))


#Q7-

ase<-pms$
  
plot3d()

###################### Exercise 4 ######################

ppca<-phyl.pca(ladderized.tree,mean.quant,method='BM',mode='cov')
eigvec<-ppca$Evec
eigval<-diag(ppca$Eval)/sum(diag(ppca$Eval))
eigval
#        PC1         PC2         PC3         PC4         PC5         PC6         PC7         PC8         PC9        PC10 
#0.475988799 0.232262041 0.093148258 0.052392977 0.035449733 0.021582167 0.019456914 0.016108058 0.011786924 0.011407753 
#       PC11        PC12        PC13        PC14        PC15 
#0.009830644 0.008085034 0.005991433 0.004533553 0.001975714
sum(eigval[1:8]) #0.9463889

##Q8- I would keep the first 8 because these 8 explain ~95% of the variation explained by all pPC axes.

ppc.scores<-ppca$S


