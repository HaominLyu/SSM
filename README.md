# Serial of Similarity Matrix (SSM) for subgenome partition: Deciphering octoploid strawberry evolution
Polyploidization has been recognized as a major force in plant evolution. With the continuous progress in sequencing technologies and genome assembly algorithms, high-quality chromosome-level assemblies of polyploid genomes have become increasingly attainable. However, accurately delineating these assemblies into subgenomes remains a challenging task, especially in case where known diploid ancestors are absent. In this study, we introduce a novel approach that leverages long terminal repeat retrotransposons (LTR-RTs) coupled with the Serial Similarity Matrix (SSM) method to assign genome assemblies to subgenomes, particularly beneficial for those without known diploid progenitor genomes. The SSM method helps identify subgenome-specific LTRs and facilitates the inference of the timing of allopolyploidization events. We then applied the SSM method to the octoploid strawberry genome. Our analysis revealed three allopolyploidization events in the evolutionary trajectory of the octoploid strawberry genome, shedding light on the evolutionary process of the origin of the octoploid strawberry genome and enhancing our understanding of allopolyploidization in this complex species.

# 1.Pairwise comparison of LTR-RT sequences for octoploid strawberry (_Fragaria × ananassa_)
#We recommended using LTR_retriever pipeline. The commands are listed as below. Merge multiple genome in a single file to explore the clusters of their chromosomes together.

  nohup gt suffixerator -db Fananassa.genome.fa -indexname Fananassa.genome.fa -tis -suf -lcp -des -ssp -sds -dna &

  nohup ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.75 Fananassa.genome.fa > Fananassa.ltrfinder.scn &

  nohup gt -j 48 ltrharvest -index Fananassa.genome.fa -similar 75 -seqids yes -vic 10 -seed 20 -motif tgca -motifmis 1 -minlenltr 100 -maxlenltr 7000 -mintsd 2 -maxtsd 6 > Fananassa.ltrharvest.scn &

  nohup LTR_retriever -genome Fananassa.genome.fa -inharvest Fananassa.ltrharvest.scn -infinder Fananassa.ltrfinder.scn -threads 64 2>> error.log &

# 2.Reconstruct the similarity matrix among all the chromosomes
#The LTR-RT sequences were filtered for the matrix according to the interval of sequence similarity.

  cat Fananassa.genome.fa.out.*.LAI.LTR.ava.out > Fananassa.genome.fa.out.LAI.LTR.ava.out &

  perl ../ClosestLTRMatch.pl Fananassa.genome.fa.out.LAI.LTR.ava.out > Fananassa.genome.fa.out.LAI.LTR.ava.out.closest &

  perl ../ClusterProByLTRBlastFull.pl *.genome.fa.out.LAI.LTR.ava.out.closest *.ctl

# 3.Clustering using R packages
#Run in Run using Fragaria x ananassa as an example:

library(pvclust)

library(MASS)


#Read similarity matrix

data1=read.table("Fananassa.SML.90.To.95.proportion",header=T,row.names=1)


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


#Drawing clustering figures

par(mfcol=c(2,2))

plot(res.pv1, print.pv=T, hang = -1, cex = 0.8, main = "pvclust (90,95)", print.num=TRUE)

pvrect(res.pv1,alpha=0.99)

plot(y1, x1, xlab="Coordinate 1", ylab="Coordinate 2",  main="NomMetric MDS (1,2)", pch=19)

text(y1, x1, labels = row.names(data2), cex=.7)

plot(z1, y1, xlab="Coordinate 2", ylab="Coordinate 3",  main="NomMetric MDS (2,3)", pch=19)

text(z1, y1, labels = row.names(data1), cex=.7)
