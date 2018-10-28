####################################################
###            Midterm Take-Home Exam            ###
####################################################
# Zachary Sperstad

master<-read.csv("master.csv")
master$spp<-substr(master$cat,1,2)
master$sex<-substr(master$cat,4,4)
hominine<-c("Gg","Ho","Pp","Pt")
hdata<-subset(master, master$spp%in%hominine)
ho.data<-subset(hdata, hdata$spp=='Ho')
gg.data<-subset(hdata, hdata$spp=='Gg')
pp.data<-subset(hdata, hdata$spp=='Pp')
pp.data<-subset(pp.data, pp.data$sex!='u')
pt.data<-subset(hdata, hdata$spp=='Pt')

ho.sub<-ho.data[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
                    133:138, 160:165, 169:171, 175:177, 184:186, 190:192)]
gg.sub<-gg.data[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
                    133:138, 160:165, 169:171, 175:177, 184:186, 190:192)]
pp.sub<-pp.data[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
                    133:138, 160:165, 169:171, 175:177, 184:186, 190:192)]
pt.sub<-pt.data[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
                    133:138, 160:165, 169:171, 175:177, 184:186, 190:192)]

require(geomorph)
p<-72/3
k<-3
ho.array<-arrayspecs(ho.sub, p, k)
gg.array<-arrayspecs(gg.sub, p, k)
pp.array<-arrayspecs(pp.sub, p, k)
pt.array<-arrayspecs(pt.sub, p, k)

ho.gpa<-gpagen(ho.array)
gg.gpa<-gpagen(gg.array)
pp.gpa<-gpagen(pp.array)
pt.gpa<-gpagen(pt.array)

ho.shape<-ho.gpa$coords #Creates an object with the GPA coordinates from the primate data.
gg.shape<-gg.gpa$coords #Creates an object with the GPA coordinates from the primate data.
pp.shape<-pp.gpa$coords #Creates an object with the GPA coordinates from the primate data.
pt.shape<-pt.gpa$coords #Creates an object with the GPA coordinates from the primate data.

ho.sex.diff<-permudist(ho.shape, groups=ho.data$sex, rounds=1000)
ho.sex.diff
#$`p.value`
#  f
#m 0.683

#$p.adjust.method
#[1] "none"

#$dist
#  f
#m 0.02124118

gg.sex.diff<-permudist(gg.shape, groups=gg.data$sex, rounds=1000)
gg.sex.diff
#$`p.value`
# f
#m 0.001

#$p.adjust.method
#[1] "none"

#$dist
#  f
#m 0.05517037

pp.sex.diff<-permudist(pp.shape, groups=pp.data$sex, rounds=1000)
pp.sex.diff
#$`p.value`
#  f
#m 0.162

#$p.adjust.method
#[1] "none"

#$dist
#  f
#m 0.02897106

pt.sex.diff<-permudist(pt.shape, groups=pt.data$sex, rounds=1000)
pt.sex.diff
#$`p.value`
#  f
#m 0.001

#$p.adjust.method
#[1] "none"

#$dist
#  f
#m 0.03364047

ho.csize<-ho.gpa$Csize
gg.csize<-gg.gpa$Csize
pp.csize<-pp.gpa$Csize
pt.csize<-pt.gpa$Csize  
  
ho.gdf<-geomorph.data.frame(ho.shape=ho.shape,ho.sex=ho.data$sex,ho.size=ho.csize)
gg.gdf<-geomorph.data.frame(gg.shape=gg.shape,gg.sex=gg.data$sex,gg.size=gg.csize)
pp.gdf<-geomorph.data.frame(pp.shape=pp.shape,pp.sex=pp.data$sex,pp.size=pp.csize)
pt.gdf<-geomorph.data.frame(pt.shape=pt.shape,pt.sex=pt.data$sex,pt.size=pt.csize)

ho.allo<-procD.allometry(ho.shape~ho.size,data=ho.gdf)
ho.allo
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
plot(ho.allo,method="RegScore",pch=19,cex=1.5,cex.lab=1.5)


gg.allo<-procD.allometry(gg.shape~gg.size,data=gg.gdf)
gg.allo
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
plot(gg.allo,method="RegScore",pch=19,cex=1.5,cex.lab=1.5)


pp.allo<-procD.allometry(pp.shape~pp.size,data=pp.gdf)
pp.allo
#Type I (Sequential) Sums of Squares and Cross-products
#Randomized Residual Permutation Procedure Used
#1000 Permutations
#ANOVA effect sizes and P-values based on empirical F distributions

#Df      SS        MS     Rsq      F       Z Pr(>F)
#log(size)  1 0.00713 0.0071300 0.02341 1.0306 0.23825  0.387
#Residuals 43 0.29749 0.0069184 0.97659                      
#Total     44 0.30462  
plot(pp.allo,method="RegScore",pch=19,cex=1.5,cex.lab=1.5)


pt.allo<-procD.allometry(pt.shape~pt.size,data=pt.gdf)
pt.allo
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
plot(pt.allo,method="RegScore",pch=19,cex=1.5,cex.lab=1.5)


spp.sex<-paste(hdata$spp,hdata$sex)
quant<-hdata[,c(16:21, 25:27, 34:36, 40:45, 55:57, 61:69, 79:84, 109:114, 118:120, 127:129, 
             133:138, 160:165, 169:171, 175:177, 184:186, 190:192)]
hominine.array<-arrayspecs(quant, p, k)
hominine.gpa<-gpagen(hominine.array)
hominine.shape<-hominine.gpa$coords
group.shape<-coords.subset(hominine.shape,group=spp.sex)
group.mean<-lapply(group.shape,mshape)

hof<-as.matrix(group.mean$`Ho f`)
hom<-as.matrix(group.mean$`Ho m`)  
ggf<-as.matrix(group.mean$`Gg f`)  
ggm<-as.matrix(group.mean$`Gg m`)  
ppf<-as.matrix(group.mean$`Pp f`)  
ppm<-as.matrix(group.mean$`Pp m`)  
ptf<-as.matrix(group.mean$`Pt f`)  
ptm<-as.matrix(group.mean$`Pt m`)  

plotRefToTarget(hof,hom,method='vector')  
plotRefToTarget(ggf,ggm,method='vector')
plotRefToTarget(ppf,ppm,method='vector')
plotRefToTarget(ptf,ptm,method='vector')
