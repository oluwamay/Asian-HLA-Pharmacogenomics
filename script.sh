#!/bin/bash
#Download the complete sample list for all populations
wget https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Global_genome_structure_project/complete_1000_genomes_sample_list_.tsv
#Download data set for 4 asian population, they have been converted to bed, bim and fam files from vcf.
wget https://github.com/HackBio-Internship/public_datasets/blob/main/Asia_HLA_Distribution/binary_plink_file/asia.bed.gz
wget https://github.com/HackBio-Internship/public_datasets/blob/main/Asia_HLA_Distribution/binary_plink_file/asia.bim
wget https://github.com/HackBio-Internship/public_datasets/raw/main/Asia_HLA_Distribution/binary_plink_file/asia.fam
mv asia.bed.gz asia.bed

#Using PLINK we will perform the following PCA analysis
#generate Eigenvalues
plink --bfile asia --pca 
#This generates two values .eigen.vec and .eigen.val which will be plotted on R

###Codes used on R
# Read
pca1 <- read.table("plink.eigenvec",sep=" ",header=F)
#This command will read the eigen.vec file as a table and output it to pca1
#Plot the values of PC1 and PC2 which are in column 3 and 4
plot(data=pca1, V3~V4)
#ggplot
install.packages("ggplots2")
pca1 <- read.table("plink.eigenvec",sep=" ",header=F)
library("ggplot2")
#Preliminary plot
ggplot(data=pca1, aes(V3,V4)) + geom_point()
#set complete genome as PCA2 and create a metadata table
metadata <- read.table("complete_1000_genomes_sample_list_.tsv", sep ="\t", header = TRUE)
#Merge PCA1 with PCA2
merge_data <- merge(x= pca1,y = metadata, by.x = "V2", by.y = "Sample.name", all = F )

#Create a coloured plot with ggplot

 ggplot(data=merge_data, aes(V3,V4,color = Population.name)) + geom_point() + xlab("Principal Component 1 (PC1)"
                                                              + ylab("Principal Component 2 (PC2)" + ggtitle("PCA of selected Asian Populations")

#Multidimensional Scaling Analysis
##Create a LD pruned set of markers
plink --bfile asia --indep-pairwise 1000 10 0.05 --out prune1 

#Calculate identity by descent score on pruned markers
plink --bfile asia --extract prune1.prune.in --genome --out ibs1

#Cluster individuals into homogenous groups
plink --bfile asia --read-genome ibs1.genome --cluster --ppc 1e-3 --cc --mds-plot 2 --out strat1

#Set strat1.mds as mds data
mdsdata <- read.table("strat1.mds", header = TRUE)

#Merge mdsdata with metadata
merged_mds <- merge(x = mdsdata, y = metadata, by.x = "FID", by.y = "Sample.name", all = F)

#Create scatter plot for merged_mds colored by population code
ggplot(data=merged_mds, aes(C1, C2, color = population.code)) + geom_point() + ggtitle("Multidimensional Scaling Analysis")

plink --bfile asia --chr 7 --from-kb 66 --to-kb 159078 --make-bed --out asia_c7
plink --bfile asia --chr 14 --from-kb 19041  --to-kb 107273 --make-bed --out asia_c14
plink --bfile asia --chr 21 --from-kb 9442 --to-kb 48094 --make-bed --out asia_c21

#Perform linkage disequiibrium

plink --bfile asia_c7 --r2 --out asia_c7.ld 
plink --bfile asia_c14 --r2 --out asia_c14.ld 
plink --bfile asia_c21 --r2 --out asia_c21.ld 

