#####################################################
###      Lab 11: Elliptical Fourier Analyses      ###
#####################################################
# Zachary Sperstad

################### Exercise 1 ######################
orbits<-read.csv("Orb2D_upload.csv")
orbits$spp<-substr(orbits$cat,1,2)
hominines<-c("Gorilla","Homo","Pan","Pongo")
hominine_orbits<-subset(orbits,orbits$genus%in%hominines)
hominine_orbits<-droplevels(hominine_orbits)
quant<-hominine_orbits[5:58]
p<-54/2
k<-2
array<-arrayspecs(quant,p,k)
hom.out<-Out(array,fac=as.factor(hominine_orbits$genus))

hom.efa<-efourier(hom.out)
hom.efa
#99% 
#9 

#Q1- 9 harmonics are needed to ereach 99% of the harmonic
#      power.

#Q2- Because we are using 9 harmonics, we will have 36
#    coefficients, since each harmonic has four coefficents
#    (9x4=36).
#Q3-


require(RColorBrewer)
speciescol1<-brewer.pal(n=4, name='Set1')
hominine_orbits$genus<-as.factor(hominine_orbits$genus)
genus.col<-speciescol1[hominine_orbits$genus]

hom.pca<-PCA(hom.efa)
plot(hom.pca,col=genus.col)
#legend(legend=unique(hominine_orbits$genus), 
#       title = "Genera", col=unique(genus.col), pch=16, ncol=2, 
#       cex=0.70)  ###Work on this.

eigval<-hom.pca$eig
eigval
#[1] 5.572011e-01 1.432159e-01
hom.scores<-hom.pca$x
plot(hom.scores[,1],hom.scores[,2],ylim=c(-0.125,0.125),
     col=genus.col,pch=16, xlab="PC1 (55.7%)",ylab="PC2 (14.3%)",
     main="Hominine Orbital PCA")
legend("bottomright", legend=unique(hominine_orbits$genus), 
       title = "Genera", col=unique(genus.col), pch=16, ncol=2, 
       cex=0.70)
