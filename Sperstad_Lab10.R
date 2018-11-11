#####################################################
### Lab 10: Cluster Analyses and Categorical Data ###
#####################################################
# Zachary Sperstad

################### Exercise 1 ######################
master<-read.csv("master.csv") #Reads in master.csv file.
master$spp<-substr(master$cat,1,2) #Creates a variable for species identity.
master$sex<-substr(master$cat,4,4) #Creates a variabe for sex identity.
quant<-master[2:196] #Creates a subset of the data with only quantitative variables.

require(geomorph) #Allows geomorph functions to be performed.
k<-3 #Denotes that we are using 3D landmarks.
p<-195/3 #Denotes how many 3D landmarks we are using.
array<-arrayspecs(quant,p,k) #Creates an array for the 3D landmarks.
gpa<-gpagen(array) #Performs a GPA on the array.
coords<-gpa$coords #Extracts the landmarks coordinates from the GPA.
ts<-plotTangentSpace(coords) #Plots coordinates into tangent space.
scores<-ts$pc.scores #Extracts the scores from our coordinates in tangent space.

mean.cran<-aggregate(master[,2:195],list(master$spp),mean)
rownames(mean.cran)<-mean.cran$Group.1 #Makes row names the species identity.
cran.dist<-dist(mean.cran[,-1],method='euclidean')
comp.clust3<-hclust(cran.dist,method='average')
plot(comp.clust3,lwd=2,xlab=NA,main='Average')

#Q1-  

require(ape) #Allows ape functions to be performed.
nj.clust<-nj(cran.dist) #Creates a neighbor joining tree.
plot(nj.clust,main='Neighbor Joining Tree') #Plots the neighbor joining tree.

#Q2-

################### Exercise 2 ######################

sexspp<-cbind(master$spp,master$sex)
table<-matrix(ncol=3,nrow=12) #Creates an empty matrix with 3 columns and 12 rows.

for(i in 1:12){
  for(j in 1:3){
    sub<-subset(sexspp,sexspp[,1]==unique(master$spp)[[i]]&sexspp[,2]==
                  unique(master$sex)[[j]])
    table[i,j]<-nrow(sub)
  }
}

apply(table,2,as.numeric) #Makes columns numerical variables.
rownames(table)<-mean.cran$Group.1 #Applies species identities as row names.
sex<-c("f","m","u") #Creates a vector with sex identity abbreviations.
colnames(table)<-sex #Assigns sex identities as column names.

require(ca) #Allows ca functions to be performed.
sexspp.ca<-ca(table)
