#Zachary Sperstad

### Part 1: Canonical Variates Analysis (CVA)

require(geomorph)
require(Morpho)

master<-read.csv("master.csv")
spp<-substr(master$cat,1,2)
quant<-master[2:196]
p<-195/3
k<-3
array<-arrayspecs(quant, p, k)

gpa<-gpagen(array)
shape<-gpa$coords

cva<-CVA(shape,group=spp,cv=T,weighting=T)
var<-cva$Var
cv<-cva$CV
dim(cv)
#[1] 195  11

g.means<-cva$groupmeans
dists<-cva$Dist
dists

plot(dists$GroupdistEuclid, dists$GroupdistMaha) #Look at the differenc in scale.

cva
scores<-cva$CVscores

sppcol<-brewer.pal(n=12,'Set3')[as.factor(spp)]
plot(scores[,1],scores[,2],xlim=range(scores[,2]),col=sppcol,pch=19,
     xlab='CV1 (41.6%)',ylab='CV2 (24.8%)',cex.lab=1.5,cex=1.25)
legend('bottomleft',legend=unique(spp),col=unique(sppcol),ncol=3,pch=19, cex=0.75)


### Part 2: MANOVA

ts<-plotTangentSpace(shape)
pcscores<-ts$pc.scores
manova<-manova(pcscores~spp)
summary(manova)

### Part 3: Hotelling's T^2 Test

install.packages("Hotelling")
require(Hotelling)
eigval<-ts$sdev #You need to do something this to make it correct.


eigval<-cva$
sum(eigval[1:18])
hotelling.test()

