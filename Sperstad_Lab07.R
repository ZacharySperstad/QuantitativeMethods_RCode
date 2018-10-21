#######################################################################
###               Lab 07: Thin Plate Splines & PLS Redux            ###
#######################################################################
# Zachary Sperstad

############################ Exercise 1 ###############################

master<-read.csv("master.csv")
master$spp<-substr(master$cat,1,2)
master$sex<-substr(master$cat,4,4)
hp<-c("Ho", "Pt", "Pp")
hpdata<-subset(master, master$spp%in%hp)
hpdata<-droplevels(hpdata)

require(geomorph)
qdata<-hpdata[2:196]
p<-195/3
k<-3
array<-arrayspecs(qdata, p, k)

gpa<-gpagen(array)
shape<-gpa$coords
shape<-shape*-1

species<-hpdata[,197]
species.shape<-coords.subset(shape,group=species)

group.mean<-lapply(species.shape,mshape)
Ho<-as.matrix(group.mean$`Ho`)[1:65,1:3]
Pp<-as.matrix(group.mean$`Pp`)[1:65,1:3]
Pt<-as.matrix(group.mean$`Pt`)[1:65,1:3]

plotRefToTarget(Ho, Pp, method='TPS')

#Q1-  In this plot, I see two thin plate spline (TPS) deformation grids: 
#     one that shows mean landmarks for the human skull either looking
#     down onto the skull or a view from underneath the skull and the second
#     grid shows the mean landmarks of a human skull facing off to the right.
#     The cells in the grid are warped to show how you would have to move
#     the landmarks in order to match the location of the mean bonobo skull
#     landmarks. It looks like the front of the face and back of the head
#     must be compressed for human skulls to look more like bonobo skulls.

plotRefToTarget(Pp, Ho, method='TPS')

#Q2-  This would show the landmarks of the bonobo and how these landmarks
#     would have to change to match the coordinates of the human mean skull
#     landmarks. In order to have a bonobo skull match a human skull, the
#     posterior of the skull must be expanded.

plotRefToTarget(Ho, Pt, method='points')
plotRefToTarget(Ho, Pt, method='vector')

plotRefToTarget(Pp, Pt, method='points')
plotRefToTarget(Pp, Pt, method='vector')


#Q3- 

ts<-plotTangentSpace(shape)
ts.shapes<-ts$pc.shapes

plotRefToTarget(ts.shapes[1,],ts.shapes[2,], method="TPS")

#Q4-


############################ Exercise 2 ###############################

quant<-master[,2:196]
sub<-quant[,c(25:27,34:36,40:45,58:69,79:96,118:120,127:129,133:138,142:156)]

p.2<-69/3
k.2<-3
array.2<-arrayspecs(sub, p.2, k.2)

gpa.2<-gpagen(array.2)
shape.2<-gpa.2$coords
shape.2<-shape.2*-1

group.id<-c("A","A","A","A","A","A","A","A","A","A","A","A","A","A","A",
            "A","A","A","B","B","B","B","B")

int.test<-integration.test(shape.2,partition.gp=group.id)
eigval<-int.test$svd$d/sum(int.test$svd$d)

#Q5-  There are 15 eigenvalues. There are this many eigenvalues because
#     _____________.
#Q6-  The first two eigenvalues explain 83.7% of the variance.

eigvec1<-int.test$left.pls.vectors
dim(eigvec1)
#[1] 54 15
eigvec2<-int.test$right.pls.vectors
dim(eigvec2)
#[1] 15 15

#Q7-  The dimensions of these two eigenvectors are 54x15 and 15x15, so they
#     aren't the same dimension.  They are not the same because __________.

scores1<-int.test$XScores
dim(scores1)
#[1] 547  15
scores2<-int.test$YScores
dim(scores2)
#[1] 547  15

require(RColorBrewer)

speciescol1<-brewer.pal(n=9, name='Set1') #Creates a vector with 9 colors.
speciescol2<-brewer.pal(n=3, name='Set1') #Creates a vector with 3 colors.
ccol<-c(speciescol1, speciescol2) #Combines color vectors to make a vector of 12 colors.
master$spp<-as.factor(master$spp) #Makes species variable a factor.
species.col<-ccol[master$spp] #Assigns specific colors to specific species.

par(mfrow=c(1,2))
plot(scores1[,1],scores2[,1],pch=19,xlab='Face PC1',ylab='Base PC1', 
     col=species.col)
plot(scores1[,2],scores2[,2],pch=19,xlab='Face PC2',ylab='Base PC2',
     col=species.col)

### I need a goddamn legend...

#Q8-  


