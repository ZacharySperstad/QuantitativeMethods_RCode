#########################################################
###   Lab 09: CVA, MANOVA, and Hotelling's T^2 Test   ###
#########################################################
# Zachary Sperstad

###                     Exercise 1                    ###

master<-read.csv("master.csv")
master$spp<-substr(master$cat,1,2)
pp.master<-subset(master, master$spp=="Pp")
pp.sex<-substr(pp.master$cat,4,4)
quant<-pp.master[,2:196]

#Q1-  The cat IDs for the specimen with the unknown sexes are Ppnu067k, 
#     Ppnu068k, Ppnu069k.

require(geomorph)
p<-195/3
k<-3
array<-arrayspecs(quant,p,k)
gpa<-gpagen(array)
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

master$spp<-substr(master$cat,1,3)
h.spp<-c("Haa","Hal","Hau","Hma","Hmf","Hss")
h.data<-subset(master,master$spp%in%h.spp)
h.quat<-h.data[2:196]

p<-195/3
k<-3
h.array<-arrayspecs(h.quat,p,k)
h.gpa<-gpagen(h.array)
h.coords<-h.gpa$coords

h.cva<-CVA(h.coords,group=h.data$spp,
                   cv=T,weighting=T)
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

h.manova<-manova(h.scores~h.spp)
