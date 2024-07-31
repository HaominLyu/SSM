#!/usr/bin/Rscript
# Usage: Rscript this_script.r SML outputname
# SML (Similarity Matrix of LTRs) is output of ClusterProByLTRBlastFull.pl


library(pvclust)    # Clustering from correlation using pvclust
library(MASS)       # MDS clustering


########## Final Clustering
data1=read.table("Etef.SML.90.To.95.proportion",header=T,row.names=1)
set.seed(123)
res.pv1 = parPvclust(cl=NULL, data1, method.hclust = "average", 
                     method.dist = "correlation", nboot = 1000, 
                     iseed = NULL)
clusters1 <- pvpick(res.pv1)
d1 = dist(data1) # euclidean distances between the rows
fit1 = isoMDS(d1, k=3) # k is the number of dim
x1 = fit1$points[,1]
y1 = fit1$points[,2]
z1 = fit1$points[,3]


data2=read.table("Etef.SML.98.To.100.proportion",header=T,row.names=1)
set.seed(123)
res.pv2 = parPvclust(cl=NULL, data2, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)
clusters2 <- pvpick(res.pv2)
d2 = dist(data2) # euclidean distances between the rows
fit2 = isoMDS(d2, k=3) # k is the number of dim
x2 = fit2$points[,1]
y2 = fit2$points[,2]
z2 = fit2$points[,3]


#### Drawing figure
par(mfcol=c(2,2))
plot(res.pv2, print.pv=T, hang = -1, cex = 0.8, main = "pvclust (98,100)", print.num=TRUE)
pvrect(res.pv2,alpha=0.95)
plot(res.pv1, print.pv=T, hang = -1, cex = 0.8, main = "pvclust (90,95)", print.num=TRUE)
pvrect(res.pv1,alpha=0.95)
plot(y2, x2, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (1,2)", pch=19)
text(y2, x2, labels = row.names(data2), cex=.7)
plot(y1, x1, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (1,2)", pch=19)
text(y1, x1, labels = row.names(data1), cex=.7)


#######################################
data1=read.table("Fananassa.SML.90.To.95.proportion",header=T,row.names=1)
set.seed(123)
res.pv1 = parPvclust(cl=NULL, data1, method.hclust = "average", 
                     method.dist = "correlation", nboot = 1000, 
                     iseed = NULL)
clusters1 <- pvpick(res.pv1)
d1 = dist(data1) # euclidean distances between the rows
fit1 = isoMDS(d1, k=3) # k is the number of dim
x1 = fit1$points[,1]
y1 = fit1$points[,2]
z1 = fit1$points[,3]


par(mfcol=c(2,2))
plot(res.pv1, print.pv=T, hang = -1, cex = 0.8, main = "pvclust (90,95)", print.num=TRUE)
pvrect(res.pv1,alpha=0.99)
plot(y1, x1, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (1,2)", pch=19)
text(y1, x1, labels = row.names(data2), cex=.7)
plot(z1, y1, xlab="Coordinate 2", ylab="Coordinate 3",  main="NomMetric MDS (2,3)", pch=19)
text(z1, y1, labels = row.names(data1), cex=.7)

#######################################

data1=read.table("6sp9gStrawberry.SML.96.To.100.proportion",header=T,row.names=1)
set.seed(123)
res.pv1 = parPvclust(cl=NULL, data1, method.hclust = "average", 
                     method.dist = "correlation", nboot = 1000, 
                     iseed = NULL)
clusters1 <- pvpick(res.pv1)
d1 = dist(data1) # euclidean distances between the rows
fit1 = isoMDS(d1, k=3) # k is the number of dim
x1 = fit1$points[,1]
y1 = fit1$points[,2]
z1 = fit1$points[,3]


data2=read.table("6sp9gStrawberry.SML.90.To.95.proportion",header=T,row.names=1)
set.seed(123)
res.pv2 = parPvclust(cl=NULL, data2, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)
clusters2 <- pvpick(res.pv2)
d2 = dist(data2) # euclidean distances between the rows
fit2 = isoMDS(d2, k=3) # k is the number of dim
x2 = fit2$points[,1]
y2 = fit2$points[,2]
z2 = fit2$points[,3]


data3=read.table("6sp9gStrawberry.SML.85.To.89.proportion",header=T,row.names=1)
set.seed(123)
res.pv3 = parPvclust(cl=NULL, data3, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)
clusters3 <- pvpick(res.pv3)
d3 = dist(data3) # euclidean distances between the rows
fit3 = isoMDS(d3, k=3) # k is the number of dim
x3 = fit3$points[,1]
y3 = fit3$points[,2]
z3 = fit3$points[,3]



#### Drawing figure
par(mfcol=c(2,2))
plot(res.pv1, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (95,100)")
pvrect(res.pv1,alpha=0.99)
plot(res.pv2, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (90,95)")
pvrect(res.pv2,alpha=0.99)
plot(res.pv3, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (85,89)")
pvrect(res.pv3,alpha=0.99)


##############################################
data1=read.table("Strawberry6G.SML.98.To.100.proportion",header=T,row.names=1)
set.seed(123)
res.pv1 = parPvclust(cl=NULL, data1, method.hclust = "average", 
                     method.dist = "correlation", nboot = 1000, 
                     iseed = NULL)
clusters1 <- pvpick(res.pv1)
d1 = dist(data1) # euclidean distances between the rows
fit1 = isoMDS(d1, k=3) # k is the number of dim
x1 = fit1$points[,1]
y1 = fit1$points[,2]
z1 = fit1$points[,3]


data2=read.table("Strawberry6G.SML.95.To.98.proportion",header=T,row.names=1)
set.seed(123)
res.pv2 = parPvclust(cl=NULL, data2, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)
clusters2 <- pvpick(res.pv2)
d2 = dist(data2) # euclidean distances between the rows
fit2 = isoMDS(d2, k=3) # k is the number of dim
x2 = fit2$points[,1]
y2 = fit2$points[,2]
z2 = fit2$points[,3]


data3=read.table("Strawberry6G.SML.92.To.95.proportion",header=T,row.names=1)
set.seed(123)
res.pv3 = parPvclust(cl=NULL, data3, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)
clusters3 <- pvpick(res.pv3)
d3 = dist(data3) # euclidean distances between the rows
fit3 = isoMDS(d3, k=3) # k is the number of dim
x3 = fit3$points[,1]
y3 = fit3$points[,2]
z3 = fit3$points[,3]


data4=read.table("Strawberry6G.SML.89.To.92.proportion",header=T,row.names=1)
set.seed(123)
res.pv4 = parPvclust(cl=NULL, data4, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)
clusters4 <- pvpick(res.pv4)
d4 = dist(data4) # euclidean distances between the rows
fit4 = isoMDS(d4, k=3) # k is the number of dim
x4 = fit4$points[,1]
y4 = fit4$points[,2]
z4 = fit4$points[,3]


data5=read.table("Strawberry6G.SML.85.To.89.proportion",header=T,row.names=1)
set.seed(123)
res.pv5 = parPvclust(cl=NULL, data5, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)
clusters5 <- pvpick(res.pv5)
d5 = dist(data5) # euclidean distances between the rows
fit5 = isoMDS(d5, k=3) # k is the number of dim
x5 = fit5$points[,1]
y5 = fit5$points[,2]
z5 = fit5$points[,3]


data6=read.table("Strawberry6G.SML.70.To.85.proportion",header=T,row.names=1)
set.seed(123)
res.pv6 = parPvclust(cl=NULL, data6, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)
clusters6 <- pvpick(res.pv6)
d6 = dist(data6) # euclidean distances between the rows
fit6 = isoMDS(d6, k=3) # k is the number of dim
x6 = fit6$points[,1]
y6 = fit6$points[,2]
z6 = fit6$points[,3]



#### Drawing figure
par(mfcol=c(3,6))
plot(res.pv1, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (98,100)")
pvrect(res.pv1,alpha=0.99)
plot(res.pv2, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (95,98)")
pvrect(res.pv2,alpha=0.99)
plot(res.pv3, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (92,95)")
pvrect(res.pv2,alpha=0.99)

plot(x1, y1, xlab="Coordinate 1", ylab="Coordinate 2",  main="Nonmetric MDS (1,2)", pch=19)
text(x1, y1, labels = row.names(data1), cex=.7)
plot(x2, y2, xlab="Coordinate 1", ylab="Coordinate 2",  main="Nonmetric MDS (1,2)", pch=19)
text(x2, y2, labels = row.names(data2), cex=.7)
plot(x3, y3, xlab="Coordinate 1", ylab="Coordinate 2",  main="Nonmetric MDS (1,2)", pch=19)
text(x3, y3, labels = row.names(data3), cex=.7)

plot(y1, z1, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (2,3)", pch=19)
text(y1, z1, labels = row.names(data1), cex=.7)
plot(y2, z2, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (2,3)", pch=19)
text(y2, z2, labels = row.names(data2), cex=.7)
plot(y3, z3, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (2,3)", pch=19)
text(y3, z3, labels = row.names(data3), cex=.7)



plot(res.pv4, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (89,92)")
pvrect(res.pv4,alpha=0.99)
plot(res.pv5, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (85,89)")
pvrect(res.pv5,alpha=0.99)
plot(res.pv6, print.pv=T,print.num=TRUE, hang = -1, cex = 0.8, main = "pvclust (70,85)")
pvrect(res.pv6,alpha=0.99)

plot(x4, y4, xlab="Coordinate 1", ylab="Coordinate 2",  main="Nonmetric MDS (1,2)", pch=19)
text(x4, y4, labels = row.names(data4), cex=.7)
plot(x5, y5, xlab="Coordinate 1", ylab="Coordinate 2",  main="Nonmetric MDS (1,2)", pch=19)
text(x5, y5, labels = row.names(data5), cex=.7)
plot(x6, y6, xlab="Coordinate 1", ylab="Coordinate 2",  main="Nonmetric MDS (1,2)", pch=19)
text(x6, y6, labels = row.names(data6), cex=.7)

plot(y4, z4, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (2,3)", pch=19)
text(y4, z4, labels = row.names(data4), cex=.7)
plot(y5, z5, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (2,3)", pch=19)
text(y5, z5, labels = row.names(data5), cex=.7)
plot(y6, z6, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (2,3)", pch=19)
text(y6, z6, labels = row.names(data6), cex=.7)



