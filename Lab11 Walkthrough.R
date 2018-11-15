wtab<-read.csv('MosquitoWings.csv')
View(wtab)

require(geomorph)
array<-arrayspecs(wtab[-1],100,2)
dim(array)
#[1] 100   2 126

install.packages('Momocs')
require(Momocs)
mos.out<-Out(array,fac=as.factor(wtab$Species))
class(mos.out)

coo_plot(mos.out[1],main=wtab$Species[1],lwd=2,col='lightblue',
         points=T, pch= 19,cex=1)
mos.1<-efourier(mos.out[1],50)
plot(cumsum(harm_pow(mos.1)[-1]),type='o',ylab='Cumulative Hermonic Power',
     xlab='Harmonic Rank',pch=19,lwd=2)
mos.efa<-efourier(mos.out)
mos.efa

mos.pca<-PCA(mos.efa)
eigval<-mos.pca$eig

eigvec<-mos.pca$rotation
scores<-mos.pca$x

plot(scores[,1],scores[,2],ylim=range(scores[,1]),pch=19,
     xlab='PC1 (72.8%)',ylab='PC2 (10.1%)')
plot(mos.pca,cex=2,col='cornflowerblue')

which(scores[,1]==min(scores[,1]))
#shp3 
#3 
which(scores[,1]==max(scores[,1]))
#shp42 
#42 
which(scores[,2]==min(scores[,2]))
#shp124 
#124 
which(scores[,2]==max(scores[,2]))
#shp24 
#24 

par(mfrow=c(1,2))
tps_grid(mos.out[3],mos.out[42],shp.lwd=rep(2,2),shp.border=c('red','blue'), 
         legend.text = c('Mosquito 3 - Min(PC1)','Mosquito 42 - Max(PC1)'))
tps_grid(mos.out[124],mos.out[24],shp.lwd=rep(2,2),shp.border=c('red','blue'), 
         legend.text = c('Mosquito 124 - Min(PC2)','Mosquito 24 - Max(PC2)'))
         