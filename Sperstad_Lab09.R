#########################################################
###   Lab 09: CVA, MANOVA, and Hotelling's T^2 Test   ###
#########################################################
# Zachary Sperstad

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
