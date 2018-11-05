#########################################################
###   Lab 09: CVA, MANOVA, and Hotelling's T^2 Test   ###
#########################################################
# Zachary Sperstad

###                     Exercise 1                    ###

master<-read.csv("master.csv") #Reads in master.csv file.
master$spp<-substr(master$cat,1,2) #Creates a new variable for a taxonomic abbreviation.
pp.master<-subset(master, master$spp=="Pp") #Subsets the data so only bonobos are included.
pp.sex<-substr(pp.master$cat,4,4) #Creates a new variable for the sex of each individual.
quant<-pp.master[,2:196] #Makes a new data set with only quantitative data.

#Q1-  The cat IDs for the specimen with the unknown sexes are Ppnu067k, 
#     Ppnu068k, Ppnu069k.

require(geomorph) #Allows gemorph functions to be performed.
p<-195/3 #Denotes how many 3D points we will have.
k<-3 #Denotes number of dimensions.
array<-arrayspecs(quant,p,k) #Creates an array of the quanitative data.
gpa<-gpagen(array) #Performes the GPA analysis.
shape<-gpa$coords

mfshape<-coords.subset(shape,pp.sex)
ushape<-mfshape$u
install.packages("abind")
require("abind")
mfshape<-abind(mfshape[1:2])
mf.sex<-subset(pp.sex,pp.sex%in%c('m','f'))

require(Morpho)

cva<-CVA(mfshape,group=mf.sex,cv=T,weighting=T)
var<-cva$Var
var
#     Canonical root
#[1,]       1.419409

#Q2-  Our analysis only generated 1 canonical variate. We will only ever
#     have 1 minus the number of variables/groups (which ever is smaller)
#     variates. Since we have two groups (i.e. males and females), we will
#     only have one variate.
#Q3-  Since we only have one variate, it explains all of the variance
#     explained in this analysis (100%).

cva
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

u.aligned<-aperm(ushape)
dim(u.aligned)<-c(3,195)
mean.cent<-function(x){x<-x-mean(x)}
u.matrix<-apply(u.aligned,2,mean.cent)
uscores<-u.matrix%*%cva$CV

#Q5-

scores<-cva$CVscores
boxplot(scores~mf.sex, col=c("chartreuse", "cyan"), xlab="Sex", 
        ylab="Canonical Variate Scores")

#Wow...I spelled chartreuse correctly on the first try.

points(rep(1.5,3),uscores,col="purple",cex=2,pch=19)

female<-subset(scores,mf.sex=='f')
male<-subset(scores,mf.sex=='m')

mean.male<-mean(male)
mean.female<-mean(female)

male.unknown.mdist<-uscores-mean.male
male.unknown.mdist
#            CV 1
#[1,]  0.05107178
#[2,] -2.03393000
#[3,] -1.75221422

female.unknown.mdist<-uscores-mean.female
female.unknown.mdist
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

master$spp<-substr(master$cat,1,3) #Creates a variable in the master data matrix with taxonomic abbreviations.   
h.spp<-c("Haa","Hal","Hau","Hma","Hmf","Hss") #Creates a vector with the taxonomic abbreviations for the desired species.
h.data<-subset(master,master$spp%in%h.spp) #Creates a data matrix with only the hylobatids.
h.quat<-h.data[2:196] #Creates a dataset for the hylobatids with only quantitative data.

p<-195/3 #Denotes how many 3D landmarks are in the dataset.
k<-3 #Denotes how many dimensions to be cosidered.
h.array<-arrayspecs(h.quat,p,k) #Creates an array from the hylobatid data.
h.gpa<-gpagen(h.array) #Performs a GPA.
h.coords<-h.gpa$coords #Extracts the coordinates from the GPA.

h.cva<-CVA(h.coords,group=h.data$spp,
           cv=T,weighting=T) #Performs the cva analysis on the hylobatid date.
h.cva
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

h.scores<-h.cva$CVscores

require(RColorBrewer)
h.col<-brewer.pal(n=6, name='Set1') #Creates a vector with 6 colors.
h.data$spp<-as.factor(h.data$spp) #Makes species variable a factor.
species.col<-h.col[h.data$spp] #Assigns specific colors to specific species.

plot(h.scores[,1],h.scores[,2], xlim=c(-10,20),ylim=c(-15,15),
     col=species.col,pch=19)
legend("bottomright", legend=unique(h.data$spp), title = "Species", 
       col=unique(species.col), pch=16, ncol=2, cex=0.70)

h.manova<-manova(h.scores~h.data$spp)
summary(h.manova)
#              Df Pillai approx F num Df den Df    Pr(>F)    
#h.data$spp     5  4.591   323.28     25    720 < 2.2e-16 ***
#Residuals    144                                            
---
#  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#Q8-
#Q9-
  
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

haa.hal
haa.hau
haa.hma
haa.hmf
haa.hss
hal.hau
hal.hma
hal.hmf
hal.hss
hau.hma
hau.hmf
hau.hss
hma.hmf
hma.hss
hmf.hss

#Q10-