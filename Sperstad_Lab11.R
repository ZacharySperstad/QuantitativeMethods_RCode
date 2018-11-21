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

eh<-c("Austral","Paranth")
eh_orbits<-subset(orbits,orbits$genus%in%eh)
eh_orbits<-droplevels(eh_orbits)
eh_quant<-eh_orbits[5:58]
eh_array<-arrayspecs(eh_quant,p,k)
eh_out<-Out(eh_array,
            fac=as.factor(eh_orbits$genus))

eh_efa<-efourier(eh_out,9)

eh_coe<-eh_efa$coe
mean.cent<-function(x){x<-x-mean(x)} #Creates a mean-centering function.
eh_mchm<-apply(eh_coe,2,mean.cent) #Creates a mean-centered harmonic mean matrix from the exinct hominids.

extant_eigvec<-hom.pca$rotation
eh_scores<-eh_mchm%*%extant_eigvec
row.names(eh_scores)<-eh_orbits$genus

mc.hom.scores<-aggregate(hom.scores,
                         list(hominid_orbits$genus),mean)
mc.hom.scores<-mc.hom.scores[2:37]
row.names(mc.hom.scores)<-unique(hominid_orbits$genus)
all.taxa.scores<-rbind(eh_scores,mc.hom.scores)

hom.dist<-dist(all.taxa.scores,method='euclidean')
hom.dist
#           Austral    Paranth    Gorilla       Homo      Pongo
#Paranth 0.06542974                                            
#Gorilla 0.01678213 0.05340095                                 
#Homo    0.08858739 0.03599023 0.07909607                      
#Pongo   0.02840043 0.03929741 0.01806773 0.06727589           
#Pan     0.05213873 0.02760727 0.04534931 0.05050676 0.03402707

#Q8-  It looks like Paranth falls most closely to Pan (d = ~0.028) while
#     Atralipithicus seems to fall most closely to Gorilla (d = ~0.017).

all_individual_scores<-rbind(eh_scores,hom.scores)

desired_taxa<-c("Austral","Paranth","Pan","Pongo","Homo","Gorilla")
speciescol2<-brewer.pal(n=6, name='Set2')
hominids4color<-subset(orbits,orbits$genus%in%desired_taxa)
hominids4color$genus<-droplevels(hominids4color$genus)
hominids4color$genus
hominids4color$genus<-as.factor(hominids4color$genus)
hom.col<-speciescol2[hominids4color$genus]

plot(all_individual_scores[,1],all_individual_scores[,2]
     ,ylim=c(-0.125,0.125),col=hom.col,pch=16, xlab="PC1 ()",
     ylab="PC2 ()",main="Hominid Orbit PCA")
legend("bottomright", legend=unique(hominids4color$genus), 
       title = "Genera", col=unique(hom.col), pch=16, ncol=2, 
       cex=0.70)

require(Momocs)
exta.coe<-hom.efa$coe
exta.mean<-mshapes(hom.efa,FUN=mean,fac=hominid_orbits$genus)
exti.mean<-mshapes(eh_efa,FUN=mean,fac=eh_orbits$genus)

tps_grid(exti.mean$shp$Austral,exta.mean$shp$Gorilla)

tps_grid(exti.mean$shp$Paranth,exta.mean$shp$Pan)
