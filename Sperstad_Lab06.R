#######################################################################
###     Lab 06: Procrustean Analyses, Allometry, & Tangent Space    ###
#######################################################################
# Zachary Sperstad

############################ Exercise 1 ###############################

#Q1-  We will have 188 degrees of freedom. We will have this many because
#     we substract 7 from the total number of variables (i.e. 195-7=188, which
#     is equivalent to subtracting 7 from the number of 3D landmarks (65),
#     which is multiplied by the number of dimensions the point represents.

master<-read.csv("master.csv")
master$species<-substr(master$cat,1,2)
master$sex<-substr(master$cat,4,4)
cat<-master[,1]
gs<-master[,197]
sex<-master[,198]

primate<-master[2:196]
p<-195/3
k<-3

require(geomorph)

array<-arrayspecs(primate, p, k)
dim(array)
# [1]  65   3 547; There are 65 landmarks, these are 3D landmarks, and we have
# 547 individuals.

gpa<-gpagen(array)
coords<-gpa$coords

aligned<-aperm(coords)
dim(aligned)<-c(547,195)
write.csv(aligned,'aligned.csv')

#Q2-  D; I chose D based on the function description on www.rdocumentation.org.
#     As I understand it, GPA rotates, aligns, and scales the points without
#     using PCA. It looks like this function as projects the points into
#     tangent space by default; however, this option can be bypassed by
#     'proj=FALSE'.

############################ Exercise 2 ###############################

#install.packages('Morpho')

require(Morpho)
shapediff<-permudist(coords,groups=master$species, rounds=1000)
shapediff

#Q3- What this test is doing is testing for a significant difference
#    in cranial shape by looking at the values of each species and seeing
#    if the differ from other species (ANOVA). However, because this is a
#    permutation test, the initial difference between species is calculated
#    and then all of the data are pooled. Thereafter, new datasets are created
#    for each species by randomly drawing values out of the pool of data
#    without replacement. The difference between species is calculated again.
#    This process, in our case, is repeated 1,000 times. The significance
#    of differences is calculated by seeing how many times greater differences
#    were obtained by randomly drawing values from a pool of values. If only
#    one iteration had a greater difference than the original difference,
#    we would say the significant (p) is 1/1000 = 0.001.
#Q4- It looks like the greatest difference is between Ho (humans) and
#    Ph (Olive baboon), which was 0.38812743. Wow.

gdf<-geomorph.data.frame(coords=coords, species=master$species, sex=master$sex)

species.lm<-procD.lm(coords~species, data=gdf)
species.lm

#Call:
#procD.lm(f1 = coords ~ species, data = gdf)

#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#           Df      SS      MS     Rsq      F      Z  Pr(>F)   
#species    11  8.6312 0.78465 0.76632 159.49 19.125  0.001 **
#Residuals 535  2.6320 0.00492 0.23368                        
#Total     546 11.2632                                        
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Q5- Species ID explains ~76.6% of the variation in shape.

sex.lm<-procD.lm(coords~sex, data=gdf)
sex.lm

#Call:
#procD.lm(f1 = coords ~ sex, data = gdf)

#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#           Df      SS       MS     Rsq      F      Z  Pr(>F)   
#sex         2  0.2892 0.144588 0.02567 7.1675 4.4721  0.001 **
#Residuals 544 10.9740 0.020173 0.97433                        
#Total     546 11.2632                                         
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Q6- Sex explains ~2.6% of the variation in shape.

############################ Exercise 3 ###############################

gdf<-geomorph.data.frame(coords=coords, species=master$species, 
                         sex=master$sex, size=gpa$Csize)
allo<-procD.allometry(coords~size,data=gdf)
allo

#Call:
#procD.allometry(f1 = coords ~ size, data = gdf) 

#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#           Df      SS     MS     Rsq      F      Z  Pr(>F)   
#log(size)   1  3.5277 3.5277 0.31321 248.55 9.0685  0.001 **
#Residuals 545  7.7354 0.0142 0.68679                        
#Total     546 11.2632                                       
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Q7- Yes; Although much of the variation is still unexplained, the centroid
#    size values obtain from our GPA of all species explains (Rsq) ~31.3%
#    of the variation in the data. The ANOVA table indicates that this
#    is significant (0.001).

#Ask Tyler about the question above.


plot(allo,method = 'RegScore',pch=19,cex=1.5,cex.lab=1.5)

############################ Exercise 4 ###############################

ts<-plotTangentSpace(coords)

eigval<-(ts$sdev)^2
eigvec<-ts$rotation
scores<-ts$pc.scores
shapes<-ts$pc.shapes

dim(eigvec)
#[1] 195 188

#Q8- The dimensions of the eigenvectpr matrix are 195 by 188. The matrix
#    has these dimensions because ...

eigval
# [1] 8.964346e-03 4.718319e-03 1.291676e-03

#Ask Tyler about this too.

require(RColorBrewer)
require(rgl)
