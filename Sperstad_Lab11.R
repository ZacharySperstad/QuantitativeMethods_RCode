#####################################################
###      Lab 11: Elliptical Fourier Analyses      ###
#####################################################
# Zachary Sperstad

################### Exercise 1 ######################
orbits<-read.csv("Orb2D_upload.csv") #Reads in orbit data.
hominids<-c("Gorilla","Homo","Pan","Pongo") #Creates a vector of generic names of 
                                            #taxa of interest.
hominid_orbits<-subset(orbits,orbits$genus%in%hominids) #Subsets the full dataset to
                                                        #only include the the desired
                                                        #taxa.
hominid_orbits<-droplevels(hominid_orbits) #Drops unwanted levels from data matrix.
quant<-hominid_orbits[5:58] #Makes a new matrix with only quantitative data.
p<-54/2 #Denotes the number of 2D landmarks included in the analysis.
k<-2 #Denotes the number of dimensions used in the subsequent analyses.

require(geomorph) #Allows geomorph functions to be performed.
require(Rcpp) #Allows Rcpp functions to be performed.
require(rgl) #Allows rgl functions to be performed.
array<-arrayspecs(quant,p,k) #Creates an array with primate landmark data.
require(Momocs) #Allows Momocs functions to be performed.
hom.out<-Out(array,fac=as.factor(hominid_orbits$genus)) #Makes an outline using the
                                                        #landmarks.
hom.efa<-efourier(hom.out) #Performs an EFA on the extant hominid data.
#99% 
#9 

#Q1- 9 harmonics are needed to ereach 99% of the harmonic
#      power.
#Q2- Because we are using 9 harmonics, we will have 36
#    coefficients, since each harmonic has four coefficents
#    (9x4=36).
#Q3- Adding harmonics is basically adding shapes until 99% of shape
#    variation in shape is explained.

require(RColorBrewer) #Allows RColorBrewer functions to be performed.
speciescol1<-brewer.pal(n=4, name='Set1') #Creates a vector with 4 colors.
hominid_orbits$genus<-as.factor(hominid_orbits$genus) #Makes generic names factors.
genus.col<-speciescol1[hominid_orbits$genus] #Applies specific colors to each genera.

hom.pca<-PCA(hom.efa) #Performs a PCA on the EFA output.
plot(hom.pca,col=genus.col) #Plots PCA results
#legend(legend=unique(hominid_orbits$genus), 
#       title = "Genera", col=unique(genus.col), pch=16, ncol=2, 
#       cex=0.70)  #I couldn't get this gosh-darn code to work.

eigval<-hom.pca$eig #Extracts the eigenvalues from the PCA object.
eigval #Shows eigenvalue values.
#[1] 5.572011e-01 1.432159e-01
hom.scores<-hom.pca$x #Extracts scores from the PCA object.
plot(hom.scores[,1],hom.scores[,2],ylim=c(-0.125,0.125),
     col=genus.col,pch=16, xlab="PC1 (55.7%)",ylab="PC2 (14.3%)",
     main="Hominid Orbit PCA") #Plots the scores of the PCA of extant hominids.
legend("bottomright", legend=unique(hominid_orbits$genus), 
       title = "Genera", col=unique(genus.col), pch=16, ncol=2, 
       cex=0.70) #Adds a legend to the plot.

#Q4- From the first plot, it can be seen that PC1 explains 55.7% of the variation
#    while 14.3% of the variation is explained by PC2.
#Q5- Moving from left to right on PC1, it can be seen that the orbit changes from 
#    circular to an oval. On PC2 as you move from top to bottom, it seems like the
#    orbit is rotated clockwise.
#Q6- It is hard to see, but I don't notice any clustering.

lda<-LDA(hom.pca,fac=hominid_orbits$genus) #Performs a linear discriminant analysis
                                           #on the extant hominid PCA output.
cvtable<-lda$CV.tab #Creates a cross-validation table from the LDA output.
cvtable #Shows the LDA table.
#classified
#actual    Gorilla Homo Pan Pongo
#Gorilla      46    6  51    12
#Homo          2   30   9     1
#Pan          29    5 122    10
#Pongo         3    0  22    48

correct<-lda$CV.correct #Extracts the proportion of correct placements.
correct #Shows the proportion of correction placements.
#[1] 0.6212121

plot_CV(lda) #Plots the LDA results.

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

eh<-c("Austral","Paranth") #Makes a vector of extinct hominid generic names.
eh_orbits<-subset(orbits,orbits$genus%in%eh) #Extracts the extinct taxa data.
eh_orbits<-droplevels(eh_orbits) #Drops the unwanted levels.
eh_quant<-eh_orbits[5:58] #Creates a subset of only quantitative data.
eh_array<-arrayspecs(eh_quant,p,k) #Creates an array of the landmarks of extinct
                                   #taxa.
eh_out<-Out(eh_array,
            fac=as.factor(eh_orbits$genus)) #Outlines the landmarks of each taxon.
eh_efa<-efourier(eh_out,9) #Performs an EFA on the extinct taxa landmarks
eh_coe<-eh_efa$coe #Extracts the coefficient from the EFA on the extinct taxa.
mean.cent<-function(x){x<-x-mean(x)} #Creates a mean-centering function.
eh_mchm<-apply(eh_coe,2,mean.cent) #Creates a mean-centered harmonic mean matrix from the exinct hominids.

extant_eigvec<-hom.pca$rotation #Extracts the eigenvectors from the extant hominid
                                #PCA.
eh_scores<-eh_mchm%*%extant_eigvec #Multiples the matrices of extant hominid
                                   #eigenvectors with the mean-centered harmonic
                                   #means of the extinct hominids.
row.names(eh_scores)<-eh_orbits$genus #Applies row names to the extinct taxa score
                                      #matrix.
mc.hom.scores<-aggregate(hom.scores,
                         list(hominid_orbits$genus),mean) #Finds means for each
                                                          #genus of extant taxa.
mc.hom.scores<-mc.hom.scores[2:37] #Creates a matrix of just the scores of extant
                                   #taxa.
row.names(mc.hom.scores)<-unique(hominid_orbits$genus) #Applies row names to the
                                                       #matrix of scores of extant
                                                       #taxa.
all.taxa.scores<-rbind(eh_scores,mc.hom.scores) #Binds the matrices of scores for
                                                #extinct and extant taxa.
hom.dist<-dist(all.taxa.scores,method='euclidean') #Finds the Euclidean distances
                                                   #between the scores of all taxa
                                                   #included in the analysis.
hom.dist #Shows the results of computing the Euclidean distances between genera.
#           Austral    Paranth    Gorilla       Homo      Pongo
#Paranth 0.06542974                                            
#Gorilla 0.01678213 0.05340095                                 
#Homo    0.08858739 0.03599023 0.07909607                      
#Pongo   0.02840043 0.03929741 0.01806773 0.06727589           
#Pan     0.05213873 0.02760727 0.04534931 0.05050676 0.03402707

#Q8-  It looks like Paranth falls most closely to Pan (d = ~0.028) while
#     Atralipithicus seems to fall most closely to Gorilla (d = ~0.017).

all_individual_scores<-rbind(eh_scores,hom.scores) #Combines the scores of all taxa.

desired_taxa<-c("Austral","Paranth","Pan","Pongo","Homo","Gorilla") #Creates a vector
                                                                    #of desired taxa.
speciescol2<-brewer.pal(n=6, name='Set2') #Creates a vector with 6 colors.
hominids4color<-subset(orbits,orbits$genus%in%desired_taxa) #Subsets the original
                                                            #dataset to only include
                                                            #the desired taxa.
hominids4color$genus<-droplevels(hominids4color$genus) #Drops unwanted levels.
hominids4color$genus #Checks to make sure unwanted levels are dropped.
hominids4color$genus<-as.factor(hominids4color$genus) #Makes the generic names factors.
hom.col<-speciescol2[hominids4color$genus] #Applies specific colors to different genera.

plot(all_individual_scores[,1],all_individual_scores[,2]
     ,ylim=c(-0.125,0.125),col=hom.col,pch=16, xlab="PC1 ()",
     ylab="PC2 ()",main="Hominid Orbit PCA") #Plots hominid scores (extinct and
                                             #extant).
legend("bottomright", legend=unique(hominids4color$genus), 
       title = "Genera", col=unique(hom.col), pch=16, ncol=2, 
       cex=0.70) #Adds a legend to the plot of hominids scores.

exta.coe<-hom.efa$coe #Extracts the coefficients from the EFA on extant hominids.
exta.mean<-mshapes(hom.efa,FUN=mean,fac=hominid_orbits$genus) #Finds the mean shape
                                                              #for each extant genus.
exti.mean<-mshapes(eh_efa,FUN=mean,fac=eh_orbits$genus) #Finds the mean shape for
                                                        #each extinct genus.

tps_grid(exti.mean$shp$Austral,exta.mean$shp$Gorilla) #Makes a TPS grid showing the
                                                      #the transformation from
                                                      #Australopithecus to Gorilla.
tps_grid(exti.mean$shp$Paranth,exta.mean$shp$Pan) #Makes a TPS grid showing the
                                                  #transformation from Paranthropus
                                                  #to Pan.
