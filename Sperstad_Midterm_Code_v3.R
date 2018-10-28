####################################################
###            Midterm Take-Home Exam            ###
####################################################
# Zachary Sperstad

master<-read.csv("master.csv") #Reads in master.csv file.
master$spp<-substr(master$cat,1,2) #Creates a variable for species ID.
master$sex<-substr(master$cat,4,4) #Creates a variable for sex.
hominine<-c("Gg","Ho","Pp","Pt") #Creates a vector of the taxonomic abbreviations for the taxa of interest.
hdata<-subset(master, master$spp%in%hominine) #Subsets the original dataset to include only taxa of interest.
ho.data<-subset(hdata, hdata$spp=='Ho') #Creates a dataset with only humans.
gg.data<-subset(hdata, hdata$spp=='Gg') #Creates a dataset with only gorillas.
pp.data<-subset(hdata, hdata$spp=='Pp') #Creates a dataset with only bonobos.
pp.data<-subset(pp.data, pp.data$sex!='u') #Removes individuals of unknown (i.e. 'u') sex.
pt.data<-subset(hdata, hdata$spp=='Pt') #Creates a dataset with only humans.

ho.sub<-ho.data[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
                    133:138, 160:165, 169:171, 175:177, 184:186, 190:192)] #Makes a subset of desired data.
gg.sub<-gg.data[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
                    133:138, 160:165, 169:171, 175:177, 184:186, 190:192)] #Makes a subset of desired data.
pp.sub<-pp.data[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
                    133:138, 160:165, 169:171, 175:177, 184:186, 190:192)] #Makes a subset of desired data.
pt.sub<-pt.data[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
                    133:138, 160:165, 169:171, 175:177, 184:186, 190:192)] #Makes a subset of desired data.

require(geomorph) #Allows geomorph fuctions to be performed.
p<-72/3 #Creates a vector stating how many 3D landmarks we intend to have.
k<-3 #Denotes the number of dimensions we desire.
ho.array<-arrayspecs(ho.sub, p, k) #Creates an array of landmarks for each human individual.
gg.array<-arrayspecs(gg.sub, p, k) #Creates an array of landmarks for each gorilla individual.
pp.array<-arrayspecs(pp.sub, p, k) #Creates an array of landmarks for each bonobo individual.
pt.array<-arrayspecs(pt.sub, p, k) #Creates an array of landmarks for each chimp individual.

ho.gpa<-gpagen(ho.array) #Runs gpa on human array.
gg.gpa<-gpagen(gg.array) #Runs gpa on gorilla array.
pp.gpa<-gpagen(pp.array) #Runs gpa on bonobo array.
pt.gpa<-gpagen(pt.array) #Runs gpa on chimp array.

ho.shape<-ho.gpa$coords #Creates an object with the GPA coordinates from the human data.
gg.shape<-gg.gpa$coords #Creates an object with the GPA coordinates from the gorilla data.
pp.shape<-pp.gpa$coords #Creates an object with the GPA coordinates from the bonobo data.
pt.shape<-pt.gpa$coords #Creates an object with the GPA coordinates from the chimp data.

ho.sex.diff<-permudist(ho.shape, groups=ho.data$sex, rounds=1000) 
ho.sex.diff #Shows results.
#$`p.value`
#  f
#m 0.683

#$p.adjust.method
#[1] "none"

#$dist
#  f
#m 0.02124118

gg.sex.diff<-permudist(gg.shape, groups=gg.data$sex, rounds=1000)
gg.sex.diff #Shows results.
#$`p.value`
# f
#m 0.001

#$p.adjust.method
#[1] "none"

#$dist
#  f
#m 0.05517037

pp.sex.diff<-permudist(pp.shape, groups=pp.data$sex, rounds=1000)
pp.sex.diff #Shows results.
#$`p.value`
#  f
#m 0.162

#$p.adjust.method
#[1] "none"

#$dist
#  f
#m 0.02897106

pt.sex.diff<-permudist(pt.shape, groups=pt.data$sex, rounds=1000)
pt.sex.diff #Shows results.
#$`p.value`
#  f
#m 0.001

#$p.adjust.method
#[1] "none"

#$dist
#  f
#m 0.03364047

ho.csize<-ho.gpa$Csize #Extracts human centroid sizes.
gg.csize<-gg.gpa$Csize #Extracts gorilla centroid sizes.
pp.csize<-pp.gpa$Csize #Extracts bonobo centroid sizes.
pt.csize<-pt.gpa$Csize #Extracts chimp centroid sizes.
  
ho.gdf<-geomorph.data.frame(ho.shape=ho.shape,ho.sex=ho.data$sex,ho.size=ho.csize) #Creates a human data frame for the test to see if shape is explained by size.
gg.gdf<-geomorph.data.frame(gg.shape=gg.shape,gg.sex=gg.data$sex,gg.size=gg.csize) #Creates a gorilla data frame for the test to see if shape is explained by size.
pp.gdf<-geomorph.data.frame(pp.shape=pp.shape,pp.sex=pp.data$sex,pp.size=pp.csize) #Creates a bonobo data frame for the test to see if shape is explained by size.
pt.gdf<-geomorph.data.frame(pt.shape=pt.shape,pt.sex=pt.data$sex,pt.size=pt.csize) #Creates a chimp data frame for the test to see if shape is explained by size.

ho.allo<-procD.allometry(ho.shape~ho.size,data=ho.gdf) 
ho.allo #Shows results
#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#Df      SS        MS     Rsq      F      Z Pr(>F)  
#log(size)  1 0.01777 0.0177672 0.03986 2.2002 2.0812  0.033 *
#Residuals 53 0.42799 0.0080752 0.96014                       
#Total     54 0.44575                                         
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
plot(ho.allo,method="RegScore",pch=19,cex=1.5,cex.lab=1.5) #Plots results.


gg.allo<-procD.allometry(gg.shape~gg.size,data=gg.gdf)
gg.allo #Shows results.
#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#Df      SS       MS    Rsq      F      Z Pr(>F)   
#log(size)  1 0.07014 0.070135 0.1178 9.2132 5.9661  0.001 **
#Residuals 69 0.52526 0.007613 0.8822                        
#Total     70 0.59540                                        
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
plot(gg.allo,method="RegScore",pch=19,cex=1.5,cex.lab=1.5) #Plots results.


pp.allo<-procD.allometry(pp.shape~pp.size,data=pp.gdf)
pp.allo #Shows results.
#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#Df      SS        MS     Rsq      F       Z Pr(>F)
#log(size)  1 0.00713 0.0071300 0.02341 1.0306 0.23825  0.387
#Residuals 43 0.29749 0.0069184 0.97659                      
#Total     44 0.30462  
plot(pp.allo,method="RegScore",pch=19,cex=1.5,cex.lab=1.5) #Plots results.


pt.allo<-procD.allometry(pt.shape~pt.size,data=pt.gdf)
pt.allo #Shows results.
#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#Df      SS       MS     Rsq      F      Z Pr(>F)   
#log(size)   1 0.04088 0.040879 0.04552 5.3414 3.9825  0.002 **
#Residuals 112 0.85717 0.007653 0.95448                        
#Total     113 0.89804                                         
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
plot(pt.allo,method="RegScore",pch=19,cex=1.5,cex.lab=1.5) #Plots results.


spp.sex<-paste(hdata$spp,hdata$sex) #Creates groups based on sex and species ID.
quant<-hdata[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
             133:138, 160:165, 169:171, 175:177, 184:186, 190:192)] #Makes new data matrix with desired data.
hominine.array<-arrayspecs(quant, p, k) #Creates an array with all hominine taxa.
hominine.gpa<-gpagen(hominine.array) #Runs GPA on  array with all hominine taxa.
hominine.shape<-hominine.gpa$coords #Extracts coords for all hominine taxa from GPA.
group.shape<-coords.subset(hominine.shape,group=spp.sex) #Applied group ID (sex and species) to coordinate data.
group.mean<-lapply(group.shape,mshape) #Finds mean landmark coordinates for each defined group.

hof<-as.matrix(group.mean$`Ho f`) #Creates a subset of data with just female humans.
hom<-as.matrix(group.mean$`Ho m`) #Creates a subset of data with just male humans.
ggf<-as.matrix(group.mean$`Gg f`) #Creates a subset of data with just female gorillas.
ggm<-as.matrix(group.mean$`Gg m`) #Creates a subset of data with just male gorillas.
ppf<-as.matrix(group.mean$`Pp f`) #Creates a subset of data with just female bonobo.
ppm<-as.matrix(group.mean$`Pp m`) #Creates a subset of data with just male bonobo.
ptf<-as.matrix(group.mean$`Pt f`) #Creates a subset of data with just female chimp.
ptm<-as.matrix(group.mean$`Pt m`) #Creates a subset of data with just male chimp.

plotRefToTarget(hof,hom,method='vector') #Plots change in human female to male mean landmark coordinates.
plotRefToTarget(ggf,ggm,method='vector') #Plots change in gorilla female to male mean landmark coordinates.
plotRefToTarget(ppf,ppm,method='vector') #Plots change in bonobo female to male mean landmark coordinates.
plotRefToTarget(ptf,ptm,method='vector') #Plots change in chimp female to male mean landmark coordinates.
