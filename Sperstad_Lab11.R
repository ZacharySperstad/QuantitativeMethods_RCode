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
