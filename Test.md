setwd("Q:\\Dropbox\\Uni\\Kandidat\\DTU Fall 2017\\02450 - Machine Learning\\Projects")

data<-read.csv("2017-09-22 VAC1506 data.csv",sep=",")
attributes<-read.csv("2017-09-22 VAC1506 data - attributes.csv",sep=",")

head(data)
head(attributes)

sum(is.na(data)) # 71 NAs in the data

summary(data) # basic summary statistics for all
summary(data)[,c(3,5,6,11,17,20,21)] # selected representative statistics

d.names<-names(data)

library("ggplot2")
library("reshape2")
library("plyr")

# View single
IFU.data<-data[,1:4]
IFU.long<-melt(IFU.data,id=c("Group","Animal"))

ggplot(IFU.long)+geom_boxplot(aes(variable,value))

# View all
ggdata<-data
head(ggdata)

ggSer<-melt(ggdata,id=d.names[-(6:10)]) # Serum IgG levels
ggplot(ggSer)+geom_boxplot(aes(variable,value))

ggOcu<-melt(ggdata,id=d.names[-(11:15)]) # Ocular IgG levels
ggplot(ggOcu)+geom_boxplot(aes(variable,value))

ggNeu<-melt(ggdata,id=d.names[-(16:18)]) # Neutralizing assay
ggplot(ggNeu)+geom_boxplot(aes(variable,value))

ggAI<-melt(ggdata,id=d.names[-(19:20)]) # Antibody avidity
ggplot(ggAI)+geom_boxplot(aes(variable,value))

ggCD<-melt(ggdata,id=d.names[-(21:22)]) # T-cell types
ggplot(ggCD)+geom_boxplot(aes(variable,value))

# Clear
ggplot(data)+geom_boxplot(aes(x=as.factor(data[,"Group"]),y=data[,"Clear"]))

# View  all correlations
plot(data)

# Log serum to baseline and Serum max to account for different schedule
Ser.Max<-apply(ggdata[,6:10],1,max,na.rm=T)
Ser.Max.Log<-log10(Ser.Max)
Ser.Base<-ggdata[,6]
FC.Ser.Max<-Ser.Max/Ser.Base
logFC.Ser<-log10(FC.Ser.Max)

plot(logFC.Ser)

# Log ocular to baseline
ggdata.ocu<-ggdata
ggdata.ocu[,11:15]<-ggdata.ocu[,11:15]/ggdata.ocu[,11] # divide with baseline
ggdata.ocu[,11:15]<-log10(ggdata.ocu[,11:15]) # log10 transform

ggOcu.log<-melt(ggdata.ocu,id=d.names[-(11:15)]) # Ocular IgG levels
ggplot(ggOcu.log)+geom_boxplot(aes(variable,value))


# Normality tests ---------------------------------------------------------
graphics.off()
par(mar=c(1,1,1,1))
par(mfrow=c(5,5))

for(i in 1:22){
  qqnorm(ggdata[,i],main=d.names[i])
  qqline(ggdata[,i])
}

# Shapiro test for normality
for(i in 1:22){
  test<-shapiro.test(ggdata[,i])
  print(test)
}


# Log of all
ggdata.norm<-ggdata
ggdata.norm<-log10(ggdata.norm) # problem with 0 becomes -Inf

# Substitute with the half of the lowest
ggdata.no.zero<-ggdata

ggdata.no.zero[ggdata.no.zero==0]<-max(ggdata.no.zero,na.rm=T) # make all 0 to max
ggdata.min<-apply(ggdata.no.zero,2,min,na.rm=T) # make a vector of the min (which is not 0)

ggdata.sub.zero<-ggdata

for(i in 1:22){
 ggdata.sub.zero[,i][ggdata.sub.zero[,i]==0]<-ggdata.min[i]/2 # substitute all 0 with ½ min for the specific column
}

ggdata.sub.zero.log<-log10(ggdata.sub.zero) # log10 transform

for(i in 1:22){
  qqnorm(ggdata.sub.zero.log[,i],main=d.names[i])
  qqline(ggdata.sub.zero.log[,i])
}

# Remove control group and test for normality
for(i in 1:22){
  qqnorm(ggdata.sub.zero.log[-(26:30),i],main=d.names[i])
  qqline(ggdata.sub.zero.log[-(26:30),i])
}


# Fold Change AI ----------------------------------------------------------

FC.AI<-data[,20]/data[,19]
plot(FC.AI)

bin.clear<-data[,5]
bin.clear[bin.clear<median(data[,5])]<-0
bin.clear[bin.clear>median(data[,5])]<-1

bin.IFU.D3<-data[,3]
bin.IFU.D3[bin.IFU.D3<median(data[,3])]<-0
bin.IFU.D3[bin.IFU.D3>median(data[,3])]<-1


# New.Data.frame ----------------------------------------------------------
ndata<-cbind(data,Ser.Max,Ser.Max.Log,logFC.Ser,FC.AI,bin.clear,bin.IFU.D3)

names(ndata)


# PCA analysis ------------------------------------------------------------
PCA.data<-data[,3:22] # exclude Group and Animal from the PCA

means<-colMeans(PCA.data,na.rm=T) # find means

# subtract the column means from each row. Transpose result since apply returns a matrix corresponding to the transposed datfinal
datzeromean<- t(apply(PCA.data,1,'-',means))

# check that column means are now zero (or are numerically different from zero by something on the order of 10^-14)
colMeans(datzeromean,na.rm=T)

# divide with the variance to standardize
variance<-apply(PCA.data,2,var,na.rm=T)
standard.data<-datzeromean/variance

# get the SVD decomposition of the standardized data
svdres <- svd(standard.data) # problem with NA

# One way to deal with NA is to place the points at the non-informative part, i.e. origin, by setting the NA to mean/variance
standard.corrected.data<-standard.data

for(i in 1:20){
  standard.corrected.data[,i][is.na(standard.data[,i])]<-means[i]/variance[i] # substitute all NA with means/var for the specific column
}

# Transpose the data so we have each attribute as row
standard.corrected.data<-t(standard.corrected.data)

# get the SVD decomposition of the standardized NA corrected data
svdres <- svd(standard.corrected.data)

# extract the singular values from the result list, svdres
singularvals <- svdres$d

# calculate the variance explained by each PC
pcvariance <- singularvals^2/sum(singularvals^2)

# plot the cumulative proportion of variance explained by the PCs
graphics.off()
plot(cumsum(pcvariance), main="Data variance explained by PCs", xlab="Number of PCs included in variance sum", ylab="Proportion of variance explained")

# Plot the PC versus each other
y <- as.numeric(as.factor(d.names[3:22]))

Z <- svdres$u%*%diag(svdres$d)

# extract the two principal components i and j
i <- 1; j <- 2;
pcx <- Z[,i]
pcy <- Z[,j]

plot(c(min(pcx), max(pcx)), c(min(pcy), max(pcy)), xlab="PC A", ylab="PC B", main="PCA Analysis", type="n")

# plot points for each pc in separate colors
cols <- colors()
for(i in sort(unique(y))){
  points(pcx[y==i], pcy[y==i], col=cols[(i+1)*10])
}

# get the order that classes were gone through and plotted in in for loop
sorted <- sort(d.names[3:22], index.return=TRUE)
# add legend
legend("topright", legend=d.names[3:22][sorted$ix], fill = cols[10*(1:5)])
