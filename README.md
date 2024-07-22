# Serial of Similarity Matrix (SSM) for subgenome partition: Deciphering octoploid strawberry evolution
Polyploidization has been recognized as a major force in plant evolution. With the continuous progress in sequencing technologies and genome assembly algorithms, high-quality chromosome-level assemblies of polyploid genomes have become increasingly attainable. However, accurately delineating these assemblies into subgenomes remains a challenging task, especially in case where known diploid ancestors are absent. In this study, we introduce a novel approach that leverages long terminal repeat retrotransposons (LTR-RTs) coupled with the Serial Similarity Matrix (SSM) method to assign genome assemblies to subgenomes, particularly beneficial for those without known diploid progenitor genomes. The SSM method helps identify subgenome-specific LTRs and facilitates the inference of the timing of allopolyploidization events. We then applied the SSM method to the octoploid strawberry genome. Our analysis revealed three allopolyploidization events in the evolutionary trajectory of the octoploid strawberry genome, shedding light on the evolutionary process of the origin of the octoploid strawberry genome and enhancing our understanding of allopolyploidization in this complex species.

# 1.Pairwise comparison of LTR-RT sequences for octoploid strawberry (_Fragaria × ananassa_)
#We recommended using LTR_retriever pipeline. The commands are listed as below. Merge multiple genome in a single file to explore the clusters of their chromosomes together.
  nohup gt suffixerator -db Fananassa.genome.fa -indexname Fananassa.genome.fa -tis -suf -lcp -des -ssp -sds -dna &
  nohup ltr_finder -D 15000 -d 1000 -L 7000 -l 100 -p 20 -C -M 0.75 Fananassa.genome.fa > Fananassa.ltrfinder.scn &
  nohup gt -j 48 ltrharvest -index Fananassa.genome.fa -similar 75 -seqids yes -vic 10 -seed 20 -motif tgca -motifmis 1 -minlenltr 100 -maxlenltr 7000 -mintsd 2 -maxtsd 6 > Fananassa.ltrharvest.scn &
  nohup LTR_retriever -genome Fananassa.genome.fa -inharvest Fananassa.ltrharvest.scn -infinder Fananassa.ltrfinder.scn -threads 64 2>> error.log &

# 2.Reconstruct the similarity matrix among all the chromosomes
#The LTR-RT sequences were filtered for the matrix according to the interval of sequence similarity
  cat Fananassa.genome.fa.out.*.LAI.LTR.ava.out > Fananassa.genome.fa.out.LAI.LTR.ava.out &
  perl ../ClosestLTRMatch.pl Fananassa.genome.fa.out.LAI.LTR.ava.out > Fananassa.genome.fa.out.LAI.LTR.ava.out.closest &
  perl ../ClusterProByLTRBlastFull.pl *.genome.fa.out.LAI.LTR.ava.out.closest *.ctl

# 3.Clustering using R packages
  
