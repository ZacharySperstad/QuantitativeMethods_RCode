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

