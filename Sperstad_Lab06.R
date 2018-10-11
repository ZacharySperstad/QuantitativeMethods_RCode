#######################################################################
###     Lab 06: Procrustean Analyses, Allometry, & Tangent Space    ###
#######################################################################
# Zachary Sperstad

############################ Exercise 1 ###############################

#Q1

master<-read.csv("master.csv")
master$genus<-substr(master$cat,1,2)
master$sex<-substr(master$cat,4,4)
cat<-master[,1]
gs<-master[,197]
sex<-master[,198]

primate<-master[2:196]
p<-195/3
k<-3

require(geomorph)

array<-arrayspecs(primate, p, k)
dim(array)
# [1]  65   3 547; There are 65 landmarks, these are 3D landmarks, and we have
# 547 individuals.

gpa<-gpagen(array)
coords<-gpa$coords

aperm()