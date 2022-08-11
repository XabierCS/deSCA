

setwd('/Users/xabier/Desktop/ipsych/reviseThesis/paper2/packagesR/deSCA/')
devtools::document()
require(deSCA)


optimizeClusterMales2<- function(dataFrame,minPtsRange=c(5,8,10,12,15,20,25,30,35,40,50,100),interClusProb=0.5){
  
  ## New
  dataFrame$mX2<-round(dataFrame$mX,2)
  dataFrame$mY2<-round(dataFrame$mY,2)
  dataFrameOld<-dataFrame
  sample2<- dataFrame[!duplicated(dataFrame[,c('mX2','mY2')]),]
  dataFrame<-data.frame(mX=sample2$mX2, mY=sample2$mY2)
  
  
  
  cls<-c()
  for (x in minPtsRange){
    print(paste0(as.character(length(minPtsRange)-which(minPtsRange==x)),' Optimizations left...'))
    cl1 <- hdbscan(dataFrame, minPts = x)
    cls<-c(cls,list(cl1))
  }
  
  
  clsQC1<-c()
  clsQC2<-data.frame()
  
  for (x in 1:length(minPtsRange)){
    if (sum(c('1','2','3') %in% names(table(cls[[x]]$cluster)))==3 & length(table(cls[[x]]$cluster))<5){ 
      clsQC1<-c(clsQC1,list(cls[[x]]))
      clsQC2<-rbind(clsQC2,data.frame(minPts=minPtsRange[x],prob=mean(cls[[x]]$membership_prob)))
    }
    best <- clsQC2[which.max(clsQC2$prob),]
    bestPos<-which(minPtsRange==best$minPts)}
  
  f1<-cls[[bestPos]]
  dataFrame$cluster<-f1$cluster
  dataFrame$prob<-f1$membership_prob
  sample2 <- subset(dataFrame,dataFrame$cluster!=0)
  
  chrX <- data.frame(aggregate(sample2[,'mX'], list(sample2$cluster), mean))
  names(chrX)<-c('group','mX')
  XXY<-as.integer(chrX$group[which.max(chrX$mX)])
  
  chrY <- data.frame(aggregate(sample2[,'mY'], list(sample2$cluster), mean))
  names(chrY)<-c('group','mY')
  XYY<-as.integer(chrY$group[which.max(chrY$mY)])
  
  clusters<-c(1,2,3)
  clusters
  clusters <- clusters[-which(clusters ==XXY)]
  clusters <- clusters[-which(clusters ==XYY)]
  XY<-clusters
  dataFrame$clusKario<-0
  dataFrame$clusKario[dataFrame$cluster==XY]<-1
  dataFrame$clusKario[dataFrame$cluster==XXY]<-2
  dataFrame$clusKario[dataFrame$cluster==XYY]<-3
  
  
  
  ## Compute samples between clusters
  dataFrame$clusKario2<-dataFrame$clusKario
  dataFrame$id<-1:dim(dataFrame)[1]
  sample2<-subset(dataFrame,dataFrame$clusKario==1) # subset XY Cluster
  
  lm3<-(lm(mY ~ mX, data = sample2)) # Compute regression line 
  b<-as.numeric(lm3$coefficients[1])
  mx<-as.numeric(lm3$coefficients[2])
  
  # Calculate samples over & under the line
  predY<-(sample2$mX*mx)+(b)
  dev<-sample2$mY-predY
  sample2$xyyDev<-dev
  
  #XYY
  mydataInterXYY<-subset(sample2,sample2$prob<interClusProb & sample2$xyyDev>0)
  #mydataInterXYY_indx<-which(sample2$prob<interClusProb & sample2$xyyDev>0)
  #dim(mydataInterXYY)
  #XXY
  mydataInterXXY<-subset(sample2,sample2$prob<interClusProb & sample2$xyyDev<0)
  #mydataInterXXY_indx<-which(sample2$prob<interClusProb & sample2$xyyDev<0)
  
  #dim(mydataInterXXY)
  
  # Plot results
  #plot(sample1$mX,sample1$mY,pch=20,col=optimizeClusterMales$clusKario+1)
  #abline(c(b,mx))
  #points(mydataInterXYY$mX,mydataInterXYY$mY,col='yellow',pch=20)
  #points(mydataInterXXY$mX,mydataInterXXY$mY,col='yellow',pch=20)
  dataFrame$clusKario3<-dataFrame$clusKario2
  
  dataFrame$clusKario3[dataFrame$id %in% mydataInterXYY$id]<-(-1)
  dataFrame$clusKario3[dataFrame$id %in% mydataInterXXY$id]<-(-1)
  
  
  ## new
  names(dataFrameOld)[which(names(dataFrameOld)=='mX')]<-'mX0'
  names(dataFrameOld)[which(names(dataFrameOld)=='mY')]<-'mY0'
  names(dataFrameOld)[which(names(dataFrameOld)=='mX2')]<-'mX'
  names(dataFrameOld)[which(names(dataFrameOld)=='mY2')]<-'mY'

  d2<-merge(dataFrameOld,dataFrame,by=c('mX','mY'))
  d3<-merge(dataFrameOld,mydataInterXYY,by=c('mX','mY'))
  d4<-merge(dataFrameOld,mydataInterXXY,by=c('mX','mY'))
  
  
  #obj1<-list(f1,dataFrame,mydataInterXYY,mydataInterXXY)
  obj1<-list(f1,d2,d3,d4)
  return(obj1)
  
}

sample1 <- simulMales(50000,0,0,0.07,0.07,0.35,0.35)


plot(sample1[,c('mX','mY')],pch=20,col=factor(sample1$kariotype),xlim=c(-0.4,0.6),ylim=c(-0.4,0.65))


dataFrame2<-data.frame(mX=sample1$mX, mY=sample1$mY)

#f1 <- optimizeClusterMales(dataFrame2, interClusProb = 0.5)
f1 <- optimizeClusterMales2(dataFrame2, interClusProb = 0.5)
head(f1[[2]][c(1:5)]) # Print cluster output
plotMaleClusters(dataFrame2,f1)
dim(dataFrame2)
start_time<- Sys.time()
## Simulation Males ####

## mean ####
start_time<- Sys.time()
r=data.frame()
for (x in seq(0.35,0.55,0.04)){
  sample1 <- simulMales(50000,0,0,0.08,0.08,x,x)
  dataFrame2<-data.frame(mX=sample1$mX, mY=sample1$mY)
  f1 <- optimizeClusterMales2(dataFrame2, interClusProb = 0.5)
  vi=sum(f1[[2]]$clusKario3!=1)
  xxy=sum(f1[[2]]$clusKario3==3)
  xyy=sum(f1[[2]]$clusKario3==2)
  xy=sum(f1[[2]]$clusKario3==1)
  r0=data.frame(p=x,visual=vi,XXY=xxy,XYY=xyy,XY=xy)
  r=rbind(r,r0)
  #plotMaleClusters(dataFrame2,f1)
  
  print(x)
}

r2=r
r$XXY<-100*(r$XXY/50)
r$XYY<-100*(r$XYY/50)

par(mar = c(4, 4, 4, 4) + 0.3)    
plot(r[,c(1,3)],pch=20,ylab='True positive rate',xlab='mean LRR dif.',main="Male",ylim=c(min(r[,c(3,4)]),100))
lines(r[,c(1,3)],pch=20 ,col='blue')
points(r[,c(1,4)],pch=20)
lines(r[,c(1,4)],pch=20 ,col='darkgreen')

par(new = TRUE)         
plot(r[,c(1,2)],pch=20,type='l' ,col='grey', axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(r$visual)))      
mtext("Visual inspection", side = 4, line = 3)



## sd ####
start_time<- Sys.time()
r=data.frame()
for (x in c(seq(0.03,0.10,0.009),0.1)){
  sample1 <- simulMales(50000,0,0,x,x,0.45,0.45)
  dataFrame2<-data.frame(mX=sample1$mX, mY=sample1$mY)
  f1 <- optimizeClusterMales2(dataFrame2, interClusProb = 0.5)
  vi=sum(f1[[2]]$clusKario3!=1)
  xxy=sum(f1[[2]]$clusKario3==3)
  xyy=sum(f1[[2]]$clusKario3==2)
  xy=sum(f1[[2]]$clusKario3==1)
  r0=data.frame(p=x,visual=vi,XXY=xxy,XYY=xyy,XY=xy)
  r=rbind(r,r0)
  #plotMaleClusters(dataFrame2,f1)
  
  print(x)
}

r2=r
r$XXY<-100*(r$XXY/50)
r$XYY<-100*(r$XYY/50)

par(mar = c(4, 4, 4, 4) + 0.3)              
plot(r[,c(1,3)],pch=20,ylab='True positive rate',xlab='sd LRR dif.',main="Male",ylim=c(min(r[,c(3,4)]),100))
lines(r[,c(1,3)],pch=20 ,col='blue')
points(r[,c(1,4)],pch=20)
lines(r[,c(1,4)],pch=20 ,col='darkgreen')

par(new = TRUE)         
plot(r[,c(1,2)],pch=20,type='l' ,col='grey', axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(r$visual)))      
mtext("Visual inspection", side = 4, line = 3)


stop_time<- Sys.time()
print(stop_time-start_time)





## Simulation Females ####
#mean ####
start_time<- Sys.time()
f=data.frame()
for (x in seq(0.35,0.62,0.04)){
  sample1 <- simulFemales(50000,0,0,0.08,0.08,x,x)
  dataFrame2<-data.frame(mX=sample1$mX, mA=sample1$mA)
  f1 <- optimizeClusterFemales(dataFrame2, interClusProb = 0.5)
  vi=sum(f1[[2]]$clusKario3!=1)
  xxx=sum(f1[[2]]$clusKario3==3)
  x0=sum(f1[[2]]$clusKario3==2)
  xx=sum(f1[[2]]$clusKario3==1)
  r0=data.frame(p=x,visual=vi,XXX=xxx,X0=x0,Xx=xx)
  f=rbind(f,r0)
  #plotFemaleClusters(dataFrame2,f1)
  
  print(x)
}

stop_time<- Sys.time()
print(stop_time-start_time)


f2=f
f$XXX<-100*(f$XXX/50)
f$X0<-100*(f$X0/50)

par(mar = c(4, 4, 4, 4) + 0.3)              
plot(f[,c(1,3)],pch=20,ylab='True positive rate',xlab='mean LRR dif.',main="Female",ylim=c(min(f[,c(3,4)]),100) )
lines(f[,c(1,3)],pch=20 ,col='blue')
points(f[,c(1,4)],pch=20)
lines(f[,c(1,4)],pch=20 ,col='darkgreen')

par(new = TRUE)         
plot(f[,c(1,2)],pch=20,type='l' ,col='grey', axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(f$visual)))      
mtext("Visual inspection", side = 4, line = 3)


#sd ####
start_time<- Sys.time()
f=data.frame()
for (x in c(seq(0.03,0.10,0.009),0.1)){
  sample1 <- simulFemales(50000,0,0,x,x,0.45,0.45)
  dataFrame2<-data.frame(mX=sample1$mX, mA=sample1$mA)
  f1 <- optimizeClusterFemales(dataFrame2, interClusProb = 0.5)
  vi=sum(f1[[2]]$clusKario3!=1)
  xxx=sum(f1[[2]]$clusKario3==3)
  x0=sum(f1[[2]]$clusKario3==2)
  xx=sum(f1[[2]]$clusKario3==1)
  r0=data.frame(p=x,visual=vi,XXX=xxx,X0=x0,Xx=xx)
  f=rbind(f,r0)
  #plotFemaleClusters(dataFrame2,f1)
  
  print(x)
}

stop_time<- Sys.time()
print(stop_time-start_time)

f
f2=f
f$XXX<-100*(f$XXX/50)
f$X0<-100*(f$X0/50)

par(mar = c(4, 4, 4, 4) + 0.3)              
plot(f[,c(1,3)],pch=20,ylab='True positive rate',xlab='sd LRR dif.',main="Female",ylim=c(min(f[,c(3,4)]),100))
lines(f[,c(1,3)],pch=20 ,col='blue')
points(f[,c(1,4)],pch=20)
lines(f[,c(1,4)],pch=20 ,col='darkgreen')

par(new = TRUE)         
plot(f[,c(1,2)],pch=20,type='l' ,col='grey', axes = FALSE, xlab = "", ylab = "")
axis(side = 4, at = pretty(range(f$visual)))      
mtext("Visual inspection", side = 4, line = 3)


stop_time<- Sys.time()
print(stop_time-start_time)

## Main figures  ####




 ## trash ####
devtools::document()
require(deSCA)
# git token: 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
#  Males ####
sample1 <- simulMales(30000,0,0,0.09,0.09,0.45,0.45)

sample1 <- simulMales(100000,0,0,0.06,0.06,0.45,0.45)
start_time<- Sys.time()
dataFrame2<-sample1[,c(1,2)]
#dataFrame<-sample1[,c(1,2)]
dim(dataFrame2)
#plot(dataFrame2,pch=20,col='blue')

f2 <- optimizeClusterMales(dataFrame2, interClusProb = 0.7)
stop_start<- Sys.time()
print(stop_start-start_time)
plotMaleClusters(dataFrame2,f2)

pcs <- prcomp(f2[[2]][,c(1:2)], center = TRUE,scale. = TRUE)
f3<-cbind(f2[[2]],(pcs[[5]]))
f3_0<-subset(f3,f3$clusKario3==0)
f3_1<-subset(f3,f3$clusKario3==1)
f3_2<-subset(f3,f3$clusKario3==2)
f3_3<-subset(f3,f3$clusKario3==3)
hist(f3_1$PC2,xlim=c(-6,6),freq = F, main='Males', xlab='PC1',
     col=adjustcolor("#4F76C3", alpha.f = 0.8),10,ylim=c(0,0.8))
hist(f3_0$PC2,add=T,freq = F,col=adjustcolor("#adadad", alpha.f = 0.8),20)
hist(f3_2$PC2,add=T,freq = F,col=adjustcolor("#37C27B", alpha.f = 0.8),5)
hist(f3_3$PC2,add=T,freq = F,col=adjustcolor("#D80913", alpha.f = 0.8),5)


# Females ####
sample1 <- simulFemales(30000,0,0,0.08,0.08,0.45,0.45)
dataFrame2<-sample1[,c(1,2)]
#dataFrame<-sample1[,c(1,2)]

plot(dataFrame2,pch=20,col='blue')

f2 <- optimizeClusterFemales(dataFrame2, interClusProb = 0.7)
plotFemaleClusters(dataFrame2,f2)

pcs <- prcomp(f2[[2]][,c(1:2)], center = TRUE,scale. = TRUE)
f3<-cbind(f2[[2]],(pcs[[5]]))
f3_0<-subset(f3,f3$clusKario3==0)
f3_1<-subset(f3,f3$clusKario3==1)
f3_2<-subset(f3,f3$clusKario3==2)
f3_3<-subset(f3,f3$clusKario3==3)
hist(f3_1$PC2,xlim=c(-6,6),freq = F, main='Females', xlab='PC1',
     col=adjustcolor("#4F76C3", alpha.f = 0.8),10,ylim=c(0,0.8))
hist(f3_0$PC2,add=T,freq = F,col=adjustcolor("#adadad", alpha.f = 0.8),30)
hist(f3_2$PC2,add=T,freq = F,col=adjustcolor("#37C27B", alpha.f = 0.8))
hist(f3_3$PC2,add=T,freq = F,col=adjustcolor("#D80913", alpha.f = 0.8))
