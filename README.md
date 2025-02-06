# Serial of Similarity Matrix (SSM) for subgenome partition: Deciphering octoploid strawberry evolution
![image](https://github.com/user-attachments/assets/928c2537-9b5b-49ce-8923-26f87cfe278c)

# Workflow of SSM
#1.Pairwise comparison of LTR-RT sequences for octoploid strawberry (_Fragaria × ananassa_)
#We recommended using LTR_retriever pipeline. The commands are listed as below. Merge multiple genome in a single file to explore the clusters of their chromosomes together. If you want to use your own LTR identified using other pipelines, please name your LTR elements with the format "chromosome name:position", such as Chr1:1234..5678. And then generate the blast tabular format file for the SSM.

#2.Reconstruct the similarity matrix among all the chromosomes
#The SSM.ctl file contains the pairs of similarity intervals and chromosome names for matrix.
#Note: the similarity intervals could be adjusted according to your own genome. Mainly considering the Ks peaks between subgenomes and genome-wide LTR-RT content. For most genome, we recommended our pre-set intervals. You can set smaller intervals to increase the resolutions only if your genome are much larger than highland cotton.


# 3.Clustering using R packages
#Run in Run using Fragaria x ananassa as an example:

library(pvclust)

library(MASS)


#Read similarity matrix

data1=read.table("SML.90.To.95.proportion",header=T,row.names=1)


#Clustering from correlation using pvclust

set.seed(123)

res.pv1 = parPvclust(cl=NULL, data1, method.hclust = "average", method.dist = "correlation", nboot = 1000, iseed = NULL)

clusters1 <- pvpick(res.pv1)


#MDS clustering

d1 = dist(data1) # euclidean distances between the rows

fit1 = isoMDS(d1, k=3) # k is the number of dim

x1 = fit1$points[,1]

y1 = fit1$points[,2]

z1 = fit1$points[,3]
