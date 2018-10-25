####################################################
###            Midterm Take-Home Exam            ###
####################################################
# Zachary Sperstad

master<-read.csv("master.csv")
master$spp<-substr(master$cat,1,2)
master$sex<-substr(master$cat,4,4)
hominine<-c("Gg","Ho","Pp","Pt")
hdata<-subset(master, master$spp%in%hominine)

quant<-hdata[,2:196]
sub<-quant[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
              133:138, 160:165, 169:171, 175:177, 184:186, 190:192)]

require(geomorph)
p<-72/3
k<-3
array<-arrayspecs(sub, p, k)
gpa<-gpagen(array)
coords<-gpa$coords #Creates an object with the GPA coordinates from the primate data.
plotAllSpecimens(coords)
plot(gpa)

consensus<-gpa$consensus


require(Morpho) #Allows Morpho functions to be performed.
shapediff<-permudist(coords,groups=hdata$spp, rounds=1000, p.adjust.method="bonferroni" )
shapediff
#$`p.value`
#      Gg    Ho    Pp
#Ho 0.006            
#Pp 0.006 0.006      
#Pt 0.006 0.006 0.006

#$p.adjust.method
#[1] "bonferroni"

#$dist
#           Gg         Ho         Pp
#Ho 0.23244152                      
#Pp 0.11624028 0.17312055           
#Pt 0.10650330 0.17574679 0.05857514


gdf<-geomorph.data.frame(coords=coords, species=hdata$spp, sex=hdata$sex) #Creates a new object with the 
#variables we are interested in examining for significant differences.

species.lm<-procD.lm(coords~species, data=gdf) #Tests for a significant difference between species.
species.lm #Shows results.
#Call:
#procD.lm(f1 = coords ~ species, data = gdf)

#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions


#           Df     SS      MS     Rsq      F      Z Pr(>F)   
#species     3 2.2236 0.74120 0.56436 122.64 11.866  0.001 **
#Residuals 284 1.7164 0.00604 0.43564                        
#Total     287 3.9400                                        
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

sex.lm<-procD.lm(coords~sex, data=gdf) #Tests for a significant difference between species.
sex.lm
#Call:
#procD.lm(f1 = coords ~ sex, data = gdf)

#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions


#           Df     SS       MS     Rsq      F      Z Pr(>F)  
#sex         2 0.0506 0.025289 0.01284 1.8531 1.7581  0.055 .
#Residuals 285 3.8895 0.013647 0.98716                       
#Total     287 3.9400                                        
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

sexandspp.lm<-procD.lm(coords~species*sex, data=gdf) #Tests for a significant difference between species.
sexandspp.lm
#Call:
#procD.lm(f1 = coords ~ species * sex, data = gdf)

#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions


#             Df     SS      MS     Rsq        F       Z Pr(>F)   
#species       3 2.2236 0.74120 0.56436 125.0045 11.9019  0.001 **
#sex           2 0.0274 0.01372 0.00696   2.3132  7.1460  0.001 **
#species:sex   3 0.0347 0.01157 0.00881   1.9505  7.9058  0.001 **
#Residuals   279 1.6543 0.00593 0.41987                           
#Total       287 3.9400                                           
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1



ts<-plotTangentSpace(coords) #Creates a bivariate plot using the scores extracted from the GPA 
#earlier.

eigval<-(ts$sdev)^2/sum((ts$sdev)^2) #Creates a vectos containing the eigenvalues for each PC.
eigvec<-ts$rotation #Extracts a matrix containing the eigenvectors for each of the PCs
scores<-ts$pc.scores #Extracts the matrix of PC scores for each specimen on the new PC axes.
shapes<-ts$pc.shapes #Extracts landmark coordinates for the maximum and minimum
# specimen on each of the new PCs.

dim(eigvec) #Shows dimensions of the eigenvector matrix
#[1] 72 65

eigval<-round(eigval, digits=3)
head(eigval)
#[1] 0.474 0.113 0.058 0.045 0.034 0.024; the first three PCs explain 
barplot(eigval, ylim=c(0,1))


require(RColorBrewer) #Allows RColorBrewerFunctions to be performed.
sppcol<-brewer.pal(n=4, name='Set1') #Creates a vector with 4 colors.
hdata$spp<-as.factor(hdata$spp) #Makes species variable a factor.
hom.col<-sppcol[hdata$spp]
pch<-c(16,25)
hdata$sex<-as.factor(hdata$sex)
sex.pch<-pch[hdata$sex]

require(rgl) #Allows rgl functions to be performed.

scores123<-scores[,1:3] #Creates a vector with the scores from the first three PCs.

plot(scores123[,1], scores[,2],xlab="PC1",ylab="PC2", pch=sex.pch,
     col=hom.col,xlim=c(-0.2,0.2),ylim=c(-0.2,0.2))
plot(scores123[,1], scores[,3],xlab="PC1",ylab="PC3", pch=sex.pch,
     col=hom.col,xlim=c(-0.15,0.2),ylim=c(-0.175,0.175))
plot(scores123[,2], scores[,3],xlab="PC2",ylab="PC3", pch=sex.pch,
     col=hom.col,xlim=c(-0.1,0.1),ylim=c(-0.1,0.1))
