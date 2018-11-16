#####################################################
###      Lab 11: Elliptical Fourier Analyses      ###
#####################################################
# Zachary Sperstad

################### Exercise 1 ######################
orbits<-read.csv("Orb2D_upload.csv")
hominids<-c("Gorilla","Homo","Pan","Pongo")
hominid_orbits<-subset(orbits,orbits$genus%in%hominids)
hominid_orbits<-droplevels(hominid_orbits)
quant<-hominid_orbits[5:58]
p<-54/2
k<-2

require(geomorph)
require(Rcpp)
require(rgl)
array<-arrayspecs(quant,p,k)
require(Momocs)
hom.out<-Out(array,fac=as.factor(hominid_orbits$genus))

hom.efa<-efourier(hom.out)
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
hominid_orbits$genus<-as.factor(hominid_orbits$genus)
genus.col<-speciescol1[hominid_orbits$genus]

hom.pca<-PCA(hom.efa)
plot(hom.pca,col=genus.col)
#legend(legend=unique(hominid_orbits$genus), 
#       title = "Genera", col=unique(genus.col), pch=16, ncol=2, 
#       cex=0.70)  ###Work on this.

eigval<-hom.pca$eig
eigval
#[1] 5.572011e-01 1.432159e-01
hom.scores<-hom.pca$x
plot(hom.scores[,1],hom.scores[,2],ylim=c(-0.125,0.125),
     col=genus.col,pch=16, xlab="PC1 (55.7%)",ylab="PC2 (14.3%)",
     main="Hominid Orbit PCA")
legend("bottomright", legend=unique(hominid_orbits$genus), 
       title = "Genera", col=unique(genus.col), pch=16, ncol=2, 
       cex=0.70)
lda<-LDA(hom.pca,fac=hominid_orbits$genus)
cvtable<-lda$CV.tab
cvtable
#classified
#actual    Gorilla Homo Pan Pongo
#Gorilla      46    6  51    12
#Homo          2   30   9     1
#Pan          29    5 122    10
#Pongo         3    0  22    48

correct<-lda$CV.correct
correct
#[1] 0.6212121

plot_CV(lda)

#Q7- It looks like it is challenging to distinguish
#    hominid genera based on their orbit shape. However,
#    a few stick out to me in the graph and table. For
#    instance, humans are often distinguishable (30). Although
#    it appears that gorillas have a great amount of
#    individuals that were able to be placed (46), nearly half
#    of the inndividuals were placed incorrectly (34). Pan seemed
#    to do fairly well, with over half of their individuals
#    being placed in Pan (122); however, many were also placed in
#    Gorilla (51).

################### Exercise 2 ######################

hominid.2<-c("Gorilla","Homo","Pan","Pongo","Austral",
             "Paranth")
hominid.2_orbits<-subset(orbits,orbits$genus%in%hominid.2)
hominid.2_orbits<-droplevels(hominid.2_orbits)
quant.2<-hominid.2_orbits[5:58]
array.2<-arrayspecs(quant.2,p,k)
hom.out.2<-Out(array.2,fac=as.factor(hominid.2_orbits$genus))

hom.efa.2<-efourier(hom.out.2,9)

harmonic.co<-hom.efa.2$coe
hom.pca.2<-PCA(hom.efa.2)
eigvec<-hom.pca.2$rotation

mean.cent<-function(x){x<-x-mean(x)} #Creates a mean-centering function.
mc.harmonic.matrix<-apply(harmonic.co,2,mean.cent)

scores<-mc.harmonic.matrix%*%eigvec

pca.scores<-hom.pca.2$x
mean.pca.scores<-apply(pca.scores,2,mean)
