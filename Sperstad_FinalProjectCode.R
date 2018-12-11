############################################################
###                     Final Project                    ###
############################################################
### Zachary Sperstad

cranial<-read.csv("neurocranium_averaged.csv") #Reads in landmark data.
behavior<-c("sd","sdr","nsd") #Creates a vector of behavior categories.
behavior.known<-subset(cranial,
                       cranial$Sand.Diving%in%behavior) #Subsets the data to include
                                            #only the labrids with known
                                            #behavior.
summary(behavior.known$Sand.Diving) #Shows counts of taxa with each behavior.
#nsd  sd sdr   u 
# 61  32   2   0 
basic.coords<-behavior.known[6:44] #Subsets the dataset to include only quantitative
#data.
require(geomorph) #Allows geomorph functions to be performed.
k<-3 #Denotes the dimentionality of the data.
p<-39/3 #Denotes how many 3D landmarks will be in the analysis.
array<-arrayspecs(basic.coords,p,k) #Creates an array of the landmarks.
gpa<-gpagen(array) #Performs a GPA on the array.
shape<-gpa$coords #Extracts the coordinates from the GPA object.
plotAllSpecimens(shape) #Plots 3D landmarks of all individuals.
csize<-gpa$Csize #Extracts centroid size from GPA object.
hist(csize,breaks=95,col='blue',xlab='Centroid Size') #Shows a histogram of
                                                      #centroid sizes.
require(Morpho) #Allows Morpho functions to be performed.
behavior.diffs<-permudist(shape,groups=behavior.known$Sand.Diving,rounds=1000) 
#Divides the taxa by behavioral designation and performs a permutation 
#test to find a significant differenc in shape.

behavior.diffs #Shows result of permutation test.
#$`p.value`
#      nsd    sd
#sd  0.001      
#sdr 0.699 0.525

#$p.adjust.method
#[1] "none"

#$dist
#           nsd         sd
#sd  0.09988606           
#sdr 0.07712307 0.08928416

gdf<-geomorph.data.frame(shape=shape,behavior=behavior.known$Sand.Diving,
                         size=csize) #Makes new group variables to run the procD.lm command.
behavior.lm<-procD.lm(shape~behavior,data=gdf) #Performs an ANOVA on the data based on behavioral group identity.
behavior.lm #Shows results of this test.
#Call:
#procD.lm(f1 = shape ~ behavior, data = gdf)

#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#          Df      SS       MS     Rsq      F    Z Pr(>F)   
#behavior   2 0.21801 0.109007 0.11644 6.0621 4.33  0.001 **
#Residuals 92 1.65432 0.017982 0.88356                      
#Total     94 1.87233                                       
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

allo<-procD.allometry(shape~size,data=gdf) #Performes an allometry test.
allo #Shows results of allometry test.
#Call:
#procD.allometry(f1 = shape ~ size, data = gdf) 

#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#          Df      SS      MS     Rsq      F      Z Pr(>F)   
#log(size)  1 0.37161 0.37161 0.19848 23.029 5.6921  0.001 **
#Residuals 93 1.50072 0.01614 0.80152                        
#Total     94 1.87233                                        
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

plot(allo,method='RegScore',pch=19,cex=1.5,cex.lab=1.5) #Plots results of
                                                        #the allometry test.
ts<-plotTangentSpace(shape) #Plots GPA coordinates onto tangent space.
eigval<-(ts$sdev)^2/sum((ts$sdev)^2) #Creates a vectos containing the eigenvalues for each PC.
eigvec<-ts$rotation #Extracts a matrix containing the eigenvectors for each of the PCs
scores<-ts$pc.scores #Extracts the matrix of PC scores for each specimen on the new PC axes.
shapes<-ts$pc.shapes #Extracts landmark coordinates for the maximum and minimum
# specimen on each of the new PCs.
eigval #Shows eigenvalue values.
eigval<-round(eigval, digits=3) #Changes decimal place.
eigval #Shows eigenvalue values again.
# [1] 0.381 0.168 0.120 0.089 0.041 0.036 0.024 0.018 0.015 0.013 0.012 0.011 0.010
#[14] 0.009 0.008 0.007 0.006 0.005 0.005 0.004 0.004 0.003 0.002 0.002 0.002 0.001
#[27] 0.001 0.001 0.001 0.001 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000
barplot(eigval,ylab='Proportion of Variance',ylim=c(0,0.4), 
        xlab="Principal Components") #Makes a scree plot.
require(RColorBrewer) #Allows RColorBrewerFunctions to be performed.
shape.col<-brewer.pal(n=3, name='Set1') #Creates a vector with 3 colors.
behavior.known$Sand.Diving<-as.factor(behavior.known$Sand.Diving) #Makes behavior variable a factor.
behavior.col<-shape.col[behavior.known$Sand.Diving] #Assigns specific colors to behavior groups.
require(rgl) #Allows rgl functions to be performed.
scores123<-scores[,1:3] #Creates a vector with the scores from the first three PCs.
plot3d(scores123,xlim=c(-0.1,0.2),ylim=c(-0.175,0.125),zlim=c(-0.075,0.225),
       col=behavior.col, size=6, pch=19, xlab="PC1 (38.1%)", ylab="PC2 (16.8%)",
       zlab="PC3 (12.0%)") #Makes a 3D plot of the 
plot(scores123[,1],scores123[,2],col=behavior.col,pch=16,ylim=c(-.25,.25),
     xlab='PC1 (38.1%)',ylab='PC2 (16.8%)') #Plots PC1 and PC2
legend(x=0.27,y=0.24,legend=unique(behavior.known$Sand.Diving),pch=16,
       col=unique(shape.col)) #Adds a legend
plot(scores123[,1],scores123[,3],col=behavior.col,pch=16,
     ylim=c(-0.25,0.25),xlab='PC1 (38.1%)',ylab='PC3 (12.0%)') #Plots PC1 and PC3.
legend(x=0.27,y=0.24,legend=unique(behavior.known$Sand.Diving),pch=16,
       col=unique(shape.col)) #Adds a legend
plot(scores123[,2],scores123[,3],col=behavior.col,pch=16,
     ylim=c(-0.17,0.17),xlab='PC2 (38.1%)',ylab='PC3 (12.0%)') #Plots PC2 and PC3 
legend(x=0.08,y=0.155,legend=unique(behavior.known$Sand.Diving),pch=16,cex=0.75,
       col=unique(shape.col),ncol=3) #Adds legend.
ID<-behavior.known$Sand.Diving #Creates an object with behavior 
                               #identifications of each individual.
ID.coord.matrix<-cbind(ID,basic.coords) #Combines the ID and coordinate matrix.
aggregated.data<-aggregate(ID.coord.matrix,list(ID.coord.matrix$ID),
                           FUN=mean)#Finds mean coordinates for each
                                    #behavioral group.
mean.quant<-aggregated.data[3:41] #Extracts only quantitative data.
mean.array<-arrayspecs(mean.quant,p,k) #Creates an array of the mean 
                                       #coordinates of each behavioral 
                                       #group
plotRefToTarget(mean.array[,,1],mean.array[,,2],method="vector") #Shows non-sand-diver to sand-diver 
                                                                 #transformation.
plotRefToTarget(mean.array[,,1],mean.array[,,3],method="vector") #Shows non-sand-diver to sand-diver relative 
                                                                 #transformation.
plotRefToTarget(mean.array[,,2],mean.array[,,3],method="vector") #Shows sand-diver to sand-diver relative
                                                                 #transformation.