#######################################################################
###           Lab 05: Generalized Procrustes Analysis               ###
#######################################################################
# Zachary Sperstad

############################ Exercise 1 ###############################

master<-read.csv("master.csv") #Reads in master.csv files.

#install.packages("geomorph") #Installs geomorph; Pounded to avoid reading in.
require(geomorph) #Allows geomorph fuctions to be run.

#Q1- The dimensions of the quantitative data in the master.csv file are 547x195.
#Q2- The dimensions of the final array will be 65x3x547.

row.names(master)<-master[,1] #Makes individual IDs the row names.
master<-master[2:196] #Makes new data matrix with only quantitative data.
k<-3 #Specifies the number of dimensions
p<-195/3 #Specifies the lumber of landmarks.
array<-arrayspecs(master, p, k) #Creates an array of the primate data.
dim(array) #Displays the dimensions of the primate array.
#[1]  65   3 547


############################ Exercise 2 ###############################

#Q3- The external appearance or physical condition of something. That is a tough question.
#Q4- The center of a geometric object. You can calculate centroid size by using the Euclidean distance
#    formula.

gpa<-gpagen(array) #Performes the GPA on the primate array.

#Q5- The geometric shapes made by the landmarks were translated, scaled, and rotated. For translation,
#    a centroid is calculated for each geometric shape (using the Euclidean distance formula). These
#    centroids are then laid ontop of each other. Next, each shape is scaled, so that the sizes are
#    roughly equivalent. Rotation is the process of turning the objects until you minimize the distances
#    between homologous landmarks.
#Q6- Ordinary Procrustes Analysis (OPA) are different that General Procrustes (GPA) distances in that 
#    OPA only uses 2 individuals and is therefore not iterative (i.e. there is only 1 good solution).

coords<-gpa$coords #Creates an object with the GPA coordinates from the primate data.
plotAllSpecimens(coords) #Makes 3D plot of primate GPA coordinates.

#Q7- This 3D plot is a skull. Since this is plotting all specimen, I guess it is also showing
#    me the variation in these points! Wowzie!

Centroids<-gpa$Csize #Makes a object of centroid sizes for individuals.
master$CSize<-Centroids #Adds centroid size as a column in the master.csv data matrix.
min(master$CSize) #21.436; Cmsf007k has the smallest centroid size value for all primates in this dataset.
max(master$CSize) #70.31287; Gggm053k has the largest centroid size value for all primates in this dataset.


#Q8- Cmsf007k has the smallest centroid size (21.436) while Gggm053k has the
#    largest (70.31287). This is not surprising to me. Centroid sizes are
#    used in GPA to remove isometric size. Since adult gorillas will natural
#    have larger heads than the blue monkey, it is not surprising that
#    the size of the gorilla skull will need to be scaled down much more
#    than the blue monkey.

master<-read.csv("master.csv") #Rereads master.csv into R.

Gorilla<-c('Gg', 'Gb') #Makes a vector with genus and species acronyms for gorilla species.
master$spp<-substr(master$cat,1,2) #Adds a column with genus and species acronym.
Gorilla_Data<-subset(master, master$spp%in%Gorilla) #Makes a smaller dataset with only gorilla species

row.names(Gorilla_Data)<-Gorilla_Data[,1] #Makes individual ID the row names.
Gorilla_Data<-Gorilla_Data[2:196] #Makes gorilla dataset only retain quantitative data.
k<-3 #Specifies the number of dimensions.
p<-195/3 #Specifies the number of landmarks.
Gorilla_array<-arrayspecs(Gorilla_Data, p, k) #Creates an array from the gorilla data.
dim(Gorilla_array) #Displays the dimensions of the gorilla array.
#[1] 65  3 86
Gorilla_gpa<-gpagen(Gorilla_array) #Performs the GPA on the Gorrila array.
Gorilla_coords<-Gorilla_gpa$coords #Creates an object with the GPA coordinates from the gorilla data.
plotAllSpecimens(Gorilla_coords) #Makes a 3D plot of Gorilla GPA coordinates. 
Gorilla_Centroids<-Gorilla_gpa$Csize #Makes an object of centroid values for gorillas.
Gorilla_Data$CSize<-Gorilla_Centroids #Adds centroid sizes values for each individual as a column.
min(Gorilla_Data$CSize) #50.63584; Gggf089k has the minimum centroid size value for gorillas.
max(Gorilla_Data$CSize) #70.31287; Gggm053k has the maximum centroid size value for gorillas.

#Wow! That was super neat. It looks like individuals of Gorilla gorilla gorilla
#have both the smallest centroid sizes and the largest! What in tarnation is
#going on here?! In all seriousness, this seems like it would be a really
#cool and useful mathematical tool to use when studying fish. For instance,
#Tetraodontiformes contains species with crazy morphologies! There are tiny
#puffers and huge Molas (Mola mola), which would make these types of analyses
#very interesting!



