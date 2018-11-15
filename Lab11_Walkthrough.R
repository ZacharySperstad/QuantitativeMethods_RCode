wtab<-read.csv('MosquitoWings.csv')
View(wtab)

require(geomorph)
array<-arrayspecs(wtab[-1],100,2)
dim(array)
#[1] 100   2 126

mos.out<-out(array,fac=as.factor(wtab$Species))
??Out
