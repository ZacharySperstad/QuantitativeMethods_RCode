#########################################################
###   Lab 09: CVA, MANOVA, and Hotelling's T^2 Test   ###
#########################################################
# Zachary Sperstad

###                     Exercise 1                    ###

master<-read.csv("master.csv") #Reads in master.csv file.
master$spp<-substr(master$cat,1,2) #Creates a new variable 
                                  #for a taxonomic abbreviation.
pp.master<-subset(master, master$spp=="Pp") #Subsets the data so only 
                                            #bonobos are included.
pp.sex<-substr(pp.master$cat,4,4) #Creates a new variable for the sex 
                                  #of each individual.
quant<-pp.master[,2:196] #Makes a new data set with only quantitative 
                         #data.

#Q1-  The cat IDs for the specimen with the unknown sexes are Ppnu067k, 
#     Ppnu068k, Ppnu069k.

require(geomorph) #Allows gemorph functions to be performed.
p<-195/3 #Denotes how many 3D points we will have.
k<-3 #Denotes number of dimensions.
array<-arrayspecs(quant,p,k) #Creates an array of the quanitative data.
gpa<-gpagen(array) #Performes the GPA analysis.
shape<-gpa$coords #Extracts the coordinates from the GPA output.

mfshape<-coords.subset(shape,pp.sex) #Subsets the coordinates  of the GPA
                                     #by sex.
ushape<-mfshape$u #Creates a new object with only the coordinates of
                  #individuals with unknown sexes.
#install.packages("abind") #Pound-ed to avoid continually installing.
require("abind") #Allows abind functions to be performed.
mfshape<-abind(mfshape[1:2])  #Combines multidimensional arrays.
mf.sex<-subset(pp.sex,pp.sex%in%c('m','f')) #Subsets data by sex.

require(Morpho) #Allows Morpho functions to be performed.

cva<-CVA(mfshape,group=mf.sex,cv=T,weighting=T) #Performs a CVA on the 
                                                #male and female coordinates.
var<-cva$Var #Extracts the canonical root information from the CVA output.
var #Shows canonical root.
#     Canonical root
#[1,]       1.419409

#Q2-  Our analysis only generated 1 canonical variate. We will only ever
#     have 1 minus the number of variables/groups (which ever is smaller)
#     variates. Since we have two groups (i.e. males and females), we will
#     only have one variate.
#Q3-  Since we only have one variate, it explains all of the variance
#     explained in this analysis (100%).

cva #Shows cross-validation table for males and females.
#cross-validated classification results in frequencies
#   f  m
#f 17  7
#m  8 13

#cross-validated classification result in %
#       f      m
#f 70.833 29.167
#m 38.095 61.905

#overall classification accuracy: 66.66667 %
#Kappa statistic: 0.32836

#Q4- Males were identified correctly 61.9% of the time. Females were
#    identified correctly 70.8% of the time.

u.aligned<-aperm(ushape) #Transpose an array by permuting its dimensions.
dim(u.aligned)<-c(3,195) #Changes the dimensionality of the data for
                         #the individuals with unknown sexes.
mean.cent<-function(x){x<-x-mean(x)} #Creates a mean-centering function.
u.matrix<-apply(u.aligned,2,mean.cent) #Mean centers the coordinates of
                                       #individuals with unknown sexes.
uscores<-u.matrix%*%cva$CV #Multiplies the mean-centered coordinates of the 
                           #individuals with unknown sexes by the matrix
                           #of canonical variates.

#Q5-

scores<-cva$CVscores #Extracts scores from the CVA output.
boxplot(scores~mf.sex, col=c("chartreuse", "cyan"), xlab="Sex", 
        ylab="Canonical Variate Scores") #Creates a boxplot showing the 
                                         #distribution of male and female 
                                         #scores from the CVA.

#Wow...I spelled chartreuse correctly on the first try.

points(rep(1.5,3),uscores,col="purple",cex=2,pch=19) #Plots unknown scores 
                                                     #as points on graph.

female<-subset(scores,mf.sex=='f') #Subsets score matrix to only include 
                                   #female scores.
male<-subset(scores,mf.sex=='m') #Subsets scores matrix to only include 
                                 #male scores.

mean.male<-mean(male) #Calculates the mean of male scores.
mean.female<-mean(female) #Calculates the mean of female scores.

male.unknown.mdist<-uscores-mean.male #Finds distances between unknowns 
                                      #and mean male scores.
male.unknown.mdist #Shows distances.
#            CV 1
#[1,]  0.05107178
#[2,] -2.03393000
#[3,] -1.75221422

female.unknown.mdist<-uscores-mean.female #Finds distance between unknowns 
                                          #and the mean of female scores.
female.unknown.mdist #Shows distances.
#CV 1
#[1,] 2.3854921
#[2,] 0.3004903
#[3,] 0.5822061

#Q6- Based on the plot and distances I just calculated, I would assign the
#    first unknown specimen to the male catagory and the other two to the
#    female catergory. The first is certainly out of the distribution of
#    the female range, while the other two are within the female AND male
#    range. We should be fairly confident in this assignment. Secondly,
#    the two other specimen are very close to the female average and are
#    relatively close to the end of one of the tails of the male distribution.
#    We should be fairly consident assigning these two individuals to the
#    female category.


###                              Exercise 2                           ###

master$spp<-substr(master$cat,1,3) #Creates a variable in the master data 
                                   #matrix with taxonomic abbreviations.   
h.spp<-c("Haa","Hal","Hau","Hma","Hmf","Hss") #Creates a vector with the 
                                              #taxonomic abbreviations for 
                                              #the desired species.
h.data<-subset(master,master$spp%in%h.spp) #Creates a data matrix with only 
                                           #the hylobatids.
h.quat<-h.data[2:196] #Creates a dataset for the hylobatids with only 
                      #quantitative data.

p<-195/3 #Denotes how many 3D landmarks are in the dataset.
k<-3 #Denotes how many dimensions to be cosidered.

h.array<-arrayspecs(h.quat,p,k) #Creates an array from the hylobatid data.
h.gpa<-gpagen(h.array) #Performs a GPA.
h.coords<-h.gpa$coords #Extracts the coordinates from the GPA.

h.cva<-CVA(h.coords,group=h.data$spp,
           cv=T,weighting=T) #Performs the cva analysis on the hylobatid 
                             #date.
h.cva #Shows cross-validation table for CVA.
# cross-validated classification results in frequencies
#    Haa Hal Hau Hma Hmf Hss
#Haa  17   2   1   3   2   0
#Hal   2   8   0   5   2   0
#Hau   2   1   5   4   1   0
#Hma   4   7   4  21   3   0
#Hmf   2   2   2   5   7   0
#Hss   0   1   0   0   0  37

#cross-validated classification result in %
#        Haa     Hal     Hau     Hma     Hmf     Hss
#Haa 68.0000  8.0000  4.0000 12.0000  8.0000  0.0000
#Hal 11.7647 47.0588  0.0000 29.4118 11.7647  0.0000
#Hau 15.3846  7.6923 38.4615 30.7692  7.6923  0.0000
#Hma 10.2564 17.9487 10.2564 53.8462  7.6923  0.0000
#Hmf 11.1111 11.1111 11.1111 27.7778 38.8889  0.0000
#Hss  0.0000  2.6316  0.0000  0.0000  0.0000 97.3684


#overall classification accuracy: 63.33333 %
#Kappa statistic: 0.54555

#Q7-  It looks like Hss is distinct from the other hylobatids (97.4%
#     of assignments were correctly), while Hmf (38.9%) and and Hau 
#     (38.5%) were assigned correctly the least amount of times (i.e.
#     difficult to diagnose).

h.scores<-h.cva$CVscores #Extracts CV scores from the CVA analysis.

require(RColorBrewer) #Allows RColorBrewer functions to be performed.
h.col<-brewer.pal(n=6, name='Set1') #Creates a vector with 6 colors.
h.data$spp<-as.factor(h.data$spp) #Makes species variable a factor.
species.col<-h.col[h.data$spp] #Assigns specific colors to specific species.

h.cva$Var #Gives a summary of the canonical variates obtained.
#Canonical roots % Variance Cumulative %
#[1,]      123.627051  76.583691     76.58369
#[2,]       12.510708   7.750053     84.33374
#[3,]       10.211811   6.325947     90.65969
#[4,]        8.581967   5.316302     95.97599
#[5,]        6.495850   4.024007    100.00000

plot(h.scores[,1],h.scores[,2], xlim=c(-10,20),ylim=c(-15,15),
     col=species.col,pch=19,xlab="CV1 (76.6%)", ylab="CV2 (7.8%)", 
     main="Hylobatid CVA") #Plots the scores from the hylobatid CVA.
legend("bottomright", legend=unique(h.data$spp), 
       title = "Species", col=unique(species.col), pch=16, ncol=2, 
       cex=0.70) #Adds a sweet legend to my graph showing which color 
                 #goes with which species.

h.manova<-manova(h.scores~h.data$spp) #Performs MANOVA with the CVA scores 
                                      #as the response variable and species 
                                      #identity as the predictor variable.
summary(h.manova) #Shows MANOVA results.
#              Df Pillai approx F num Df den Df    Pr(>F)    
#h.data$spp     5  4.591   323.28     25    720 < 2.2e-16 ***
#Residuals    144                                            
---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  
#Q8-  I am reporting Pillai's Trace. I did this because the sample size
#     was very low and we have an unequal number of males and females.
#Q9-  From our MANOVA, it was demonstrated that the shape of the skulls
#     of hylobatid subspecies differed significantly based on subspecies
#     identity. This was supported by a p-value < 2.2e-16. Although this
#     tells us there is a significant difference within the group, it
#     does not tell us where this difference is. We will need to perform
#     additional tests to see which groups differ.
  
. #I had to put this period here or my code wouldn't run.
haa.scores<-subset(h.scores,h.data$spp=='Haa') #Creates a taxon specific subset of scores.
hal.scores<-subset(h.scores,h.data$spp=='Hal') #Creates a taxon specific subset of scores.
hau.scores<-subset(h.scores,h.data$spp=='Hau') #Creates a taxon specific subset of scores.
hma.scores<-subset(h.scores,h.data$spp=='Hma') #Creates a taxon specific subset of scores.
hmf.scores<-subset(h.scores,h.data$spp=='Hmf') #Creates a taxon specific subset of scores.
hss.scores<-subset(h.scores,h.data$spp=='Hss') #Creates a taxon specific subset of scores.

require(Hotelling) #Allows Hotelling Functions to be perfomed.

haa.hal<-hotelling.test(haa.scores,hal.scores) #Performs Hotelling Test on a subset of the taxa.
haa.hau<-hotelling.test(haa.scores,hau.scores) #Performs Hotelling Test on a subset of the taxa.
haa.hma<-hotelling.test(haa.scores,hma.scores) #Performs Hotelling Test on a subset of the taxa.
haa.hmf<-hotelling.test(haa.scores,hmf.scores) #Performs Hotelling Test on a subset of the taxa.
haa.hss<-hotelling.test(haa.scores,hss.scores) #Performs Hotelling Test on a subset of the taxa.
hal.hau<-hotelling.test(hal.scores,hau.scores) #Performs Hotelling Test on a subset of the taxa.
hal.hma<-hotelling.test(hal.scores,hma.scores) #Performs Hotelling Test on a subset of the taxa.
hal.hmf<-hotelling.test(hal.scores,hmf.scores) #Performs Hotelling Test on a subset of the taxa.
hal.hss<-hotelling.test(hal.scores,hss.scores) #Performs Hotelling Test on a subset of the taxa.
hau.hma<-hotelling.test(hau.scores,hma.scores) #Performs Hotelling Test on a subset of the taxa.
hau.hmf<-hotelling.test(hau.scores,hmf.scores) #Performs Hotelling Test on a subset of the taxa.
hau.hss<-hotelling.test(hau.scores,hss.scores) #Performs Hotelling Test on a subset of the taxa.
hma.hmf<-hotelling.test(hma.scores,hmf.scores) #Performs Hotelling Test on a subset of the taxa.
hma.hss<-hotelling.test(hma.scores,hss.scores) #Performs Hotelling Test on a subset of the taxa.
hmf.hss<-hotelling.test(hmf.scores,hss.scores) #Performs Hotelling Test on a subset of the taxa.

haa.hal #p = 0
haa.hau #p = 0
haa.hma #p = 0
haa.hmf #p = 0
haa.hss #p = 0
hal.hau #p = 0
hal.hma #p = 0
hal.hmf #p = 0
hal.hss #p = 0
hau.hma #p = 0
hau.hmf #p = 0
hau.hss #p = 0
hma.hmf #p = 0
hma.hss #p = 0
hmf.hss #p = 0

#Q10- From this test, it looks like all subspecies differ significantly
#     from each other, which is surprising, given the CVA scatterplot
#     obtained makes it appear as though Hal, Hma, Hau, and Hmf overlap.