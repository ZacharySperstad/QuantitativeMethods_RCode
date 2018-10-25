####################################################
###            Midterm Take-Home Exam            ###
####################################################
# Zachary Sperstad

master<-read.csv("master.csv")
master$spp<-substr(master$cat,1,2)
hominine<-c("Gg","Ho","Pp","Pt")
hdata<-subset(master, master$spp%in%hominine)
hdata$sex<-substr(hdata$cat,4,4)
