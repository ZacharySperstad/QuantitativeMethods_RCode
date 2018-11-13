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

mean.cran<-aggregate(master[,2:195],
                     list(master$spp),mean) #Finds the mean value for each variable 
                                            #for each species and makes a new table 
                                            #with only one value per variable per 
                                            #species
rownames(mean.cran)<-mean.cran$Group.1 #Makes row names the species identity.
cran.dist<-dist(mean.cran[,-1],method='euclidean') #Finds the Euclideans distances 
                                                   #between each mean variable values.
comp.clust3<-hclust(cran.dist,method='average') #Use the average linkage algorithm
                                                #to cluster species (i.e. see which
                                                #are most similar).
plot(comp.clust3,lwd=2,xlab=NA,main='Average') #Those our average linkage tree.

#Q1-  From this analysis, I see that Hs, Cm, Pb, Ha, and Hm form one large
#     cluster while Ho, Gb, Gg, Ph, Po, Pp, and Pt form another. Additionally,
#     there are smaller clusters that form, such as a gorilla cluster (Gb and Gg)
#     as well as a Pan cluster (Pp and Pt). Wow.

require(ape) #Allows ape functions to be performed.
nj.clust<-nj(cran.dist) #Creates a neighbor joining tree.
plot(nj.clust,main='Neighbor Joining Tree') #Plots the neighbor joining tree.

#Q2- They are different. First, we once again have two large clusters. However,
#    humans form a polytomy with the two large clusters (i.e. they cannot be
#    assigned to either group). In this tree gorillas once again for sister-
#    relationships, but no bonobos and chimpanzees (Pan species) and put into
#    different groups! Ha, Hm, Cm, Pb, and Hs once again belong to a "clade",
#    but Pb and Cm swap places. With all this being said, Gb, Gg, Po, Ph, and Pt
#    form a "clade" once again (with the exception of missing humans), however,
#    the relationship of taxa within this clade is even more jumbled that the
#    aforementioned clade.


################### Exercise 2 ######################

sexspp<-cbind(master$spp,master$sex) #Creates a new object with the species identity 
table<-matrix(ncol=3,nrow=12) #Creates an empty matrix with 3 columns and 12 rows.

for(i in 1:12){
  for(j in 1:3){
    sub<-subset(sexspp,sexspp[,1]==unique(master$spp)[[i]]&sexspp[,2]==
                  unique(master$sex)[[j]])
    table[i,j]<-nrow(sub)
  }
} #This for loop is looking through the sexspp object and essentially tallying
  #how many individuals are a given species AND sex.

#Q3-  Like I said above, this loop is basically looking through the sexspp object
#     for unique combinations of sex and species identities. Then, it begins
#     searching through the master dataset for individuals fitting these various
#     descriptions and counting them. The count of each individual fitting one
#     of the identity combinations is added to the object "table".

apply(table,2,as.numeric) #Makes columns numerical variables.
rownames(table)<-mean.cran$Group.1 #Applies species identities as row names.
sex<-c("f","m","u") #Creates a vector with sex identity abbreviations.
colnames(table)<-sex #Assigns sex identities as column names.

require(ca) #Allows ca functions to be performed.
sexspp.ca<-ca(table) #Performs a correspondence analysis.
summary(sexspp.ca) #Gives a summary of the correspondence analysis.
#Principal inertias (eigenvalues):
#dim    value      %   cum%   scree plot               
#1      0.114761  83.6  83.6  *********************    
#2      0.022445  16.4 100.0  ****                     
#-------- -----                                 
#Total: 0.137206 100.0                                 

#Rows:
#     name   mass  qlt  inr     k=1 cor ctr    k=2 cor ctr  
#1  |   Cm |   16 1000    6 |   205 874   6 |   78 126   4 |
#2  |   Gb |   27 1000   15 |   214 620  11 |  168 380  34 |
#3  |   Gg |  130 1000   74 |   216 592  53 |  179 408 186 |
#4  |   Ha |  101 1000   51 |   -68  67   4 |  256 933 293 |
#5  |   Hm |  104 1000   58 |  -276 996  69 |  -18   4   2 |
#6  |   Ho |  101 1000   28 |   194 993  33 |  -16   7   1 |
#7  |   Hs |   69 1000  587 | -1076 999 701 |  -27   1   2 |
#8  |   Pb |   31 1000   15 |   212 684  12 |  144 316  29 |
#9  |   Ph |   20 1000    6 |   202 926   7 |   57  74   3 |
#10 |   Po |  104 1000   30 |   190 929  33 |  -53  71  13 |
#11 |   Pp |   88 1000   16 |  -146 852  16 |  -61 148  14 |
#12 |   Pt |  208 1000  114 |   173 398  54 | -212 602 419 |
  
#Columns:
#    name   mass  qlt  inr     k=1 cor ctr    k=2 cor ctr  
#1 |    f |  475 1000   88 |    27  30   3 | -157 970 522 |
#2 |    m |  488 1000  114 |   103 330  45 |  147 670 467 |
#3 |    u |   37 1000  798 | -1729 998 952 |   84   2  12 |

#Q4- Pt has the greatest mass.

plot(sexspp.ca,mass=T,arrows=c(T,T), 
     map='rowgreen') #Plots the MCA of sex and species identity of the primates
                     #in the primate dataset.
                                     
#Q5- Hs, Hm, Pp, and Ha all appear to be close to u, meaning the have individuals
#    with unknown sexes. I doubled-checked this in the sexspp matrix and it is true
#    (pretty neat). There seems to be an abundance of females in Pt and those with
#    identifiable sexes in Pp also tend to be females. The Ha individuals that
#    were able to be sexed were predominately male, as were the gorillas (Gg and Gb)
#    and Pb (which I believe is the blue monkey).

################### Exercise 2 ######################

data("HairEyeColor") #Reads in HairEyeColor dataset.
hec.mca<-mjca(HairEyeColor) #Performs an MCA on the HairEyeColor dataset.
summary(hec.mca) #Shows the results of the MCA.
#dim    value      %   cum%   scree plot               
#1      0.054579  65.6  65.6  **********************   
#2      0.006263   7.5  73.1  ***                      
#3      0.000871   1.0  74.1                           
#-------- -----                                 
#Total: 0.083229                                       


#Columns:
#           name   mass  qlt  inr    k=1 cor ctr    k=2 cor ctr  
#1  | Hair:Black |   61  738  117 | -310 613 107 | -140 125 191 |
#2  | Hair:Blond |   72  745  126 |  527 741 364 |  -41   5  19 |
#3  | Hair:Brown |  161  691   71 |  -97 648  27 |   25  43  16 |
#4  |   Hair:Red |   40  667  119 |  -84 111   5 |  187 556 224 |
#5  |   Eye:Blue |  121  754  100 |  337 742 252 |  -42  12  35 |
#6  |  Eye:Brown |  124  724   96 | -295 688 198 |  -68  37  92 |
#7  |  Eye:Green |   36  677  121 |   86  95   5 |  214 582 262 |
#8  |  Eye:Hazel |   52  739  114 | -141 453  19 |  112 286 105 |
#9  | Sex:Female |  176  588   64 |   57 456  11 |  -31 132  27 |
#10 |   Sex:Male |  157  588   72 |  -64 456  12 |   35 132  30 |

#Q6- 65.6% of the total inertia in the dataset is explained by the first factor.

plot(hec.mca,mass=T,arrows=c(T,T), map='rowgreen') #Plots the results of the MCA.

#Q7- Okay...Here we go. People with black hair tend to have brown eyes (quadrent 3).
#    Females tend to have blue eyes and blonde hair (quadrent 4). Males tend to
#    have brown or red hair as well as hazel or green eyes. Ignoring sex, we could
#    also say, in general, people with blue eyes tend to have blonde hair (or vice
#    versa) with thos with brown or red hair tend to have green or hazel eyes (or
#    vice versa). Because black hair and brown eyes seem to be equal distances
#    away from females and males, I am going to say that these characters are
#    evenly split between males and females.
