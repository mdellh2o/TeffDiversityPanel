########################################################################################
# This script perform Phenotypic data analysis
# Date: 11/24/2020
# Authors: Matteo Dell'Acqua  and Leonardo Caproni
########################################################################################
rm(list=ls())
options(stringsAsFactors = F)

#PRELIMINARY
maindir <- "/glab/projects/teff"
setwd(maindir)

library(ggplot2)
library(ggcorrplot)
library(ggfortify)
library(patchwork)
library(tidyverse)
library(ggsci)

load(file="metadata/metadata.diversity.analysis.Rdata")
#bio<-read.csv("input/passport.data.qualitative.traits.bioclim.csv")
#summary(bio)

load(file="metadata/passport.data.qualitative.traits.bioclim_OK.Rdata")
#summary(passbio)

#make PC of biolcim
#names(passbio)

dftmp<-merge(out, passbio, by.x='ID', by.y='ID', all.x=TRUE)

biovar<-dftmp[,c(1, grep("^bio", colnames(dftmp)))]
#drop NAs
biovar<-na.omit(biovar)
dim(biovar)

#head(biovar2) 
biovar2 <- biovar[,-1]
rownames(biovar2) <- biovar[,1]
colnames(biovar2) <- sub("_27","", colnames(biovar2))

#make a df with same number of rows, this allows to plot CLusters on the PCA
tmpdf1<- merge(biovar,out, by.x='ID', by.y = 'ID', all.x = TRUE )
#out[1:10,1:10]

pcbio<-prcomp(biovar2, scale. =T, center=T)
str(pcbio)

pcbioplot <- autoplot(pcbio, 
                      scale. = TRUE, 
                      data = tmpdf1, 
                      colour = 'Cluster', 
                      loadings = TRUE,
                      loadings.colour = 'black', 
                      loadings.label = TRUE, 
                      loadings.label.size = 4
                      ) +
                ggsci::scale_fill_nejm() +
                ggsci::scale_color_nejm()+
                ggplot2::theme_light()
  

ggsave(pcbioplot, file="figures/PCA.bioclim.geno.6.clusters.jpg", height = 5, width = 6)

pcout<-pcbio$x[,1:5]
pcout<-data.frame(ID=biovar[,1], pcout)
colnames(pcout)<-paste0(colnames(pcout), "_bio")
#head(pcout)

#bring it back to passbio
passbio<-merge(passbio, pcout, by.x="ID", by.y="ID_bio", all.x=T)
#head(passbio)

#make a couple of plots
passbio1 <- subset(passbio, ID %in% out$ID)
colnames(passbio1) <- sub("_27","", colnames(passbio)) #get rid of "_27"
rawbio<-passbio1[,grep("^bio", colnames(passbio))]
rawpc<-passbio1[,grep("bio$", colnames(passbio))]

corrbioPC<-ggcorrplot(cor(rawbio, rawpc, use="pairwise.complete"), 
                      method = "circle",
                      outline.col = "gray",
                      ggtheme = ggplot2::theme_light) 
                      #ggsci::scale_color_nejm()
  
ggsave(corrbioPC, file="figures/corrBIO-PCs_OK.jpg", height=3, width =7.5)


#time to import phenotype BLUPs
#read phenotypes
ph<-read.delim("phenotypes/output_AsReml/all.BLUPs.teff.txt")

###add first 3 PCs for ph (only combined measures, separated by metric traits and farmer traitas)
comb<-merge(out, passbio1, by.x="ID", by.y="ID", all.x = TRUE)
head(comb)
#names(comb)
comb<-comb[,-c(11:24)]

comb1<-merge(comb, ph, by="ID", all.x=T)
head(comb1)
#dim(comb1)

names(comb1)
f<-c("AEZ_type", "Cluster","REGION.x","AEZ31", "NV", "NVc", 
     "Plc", "Gfpc", "PCMdate", "SBP", "Ptype", "SC", "PCHdate", "Edatec", "Hdatec", "Mdatec")

comb1[f]<-lapply(comb1[f], as.factor)
lapply(comb1, class)

#get out names by clusters
byclust<-split(comb1, comb1[,"Cluster"])
idclust<-lapply(byclust, function(x) x[,1])  

#remove all localtion-specific traits

#comb.comb<-comb1[, -grep("Ak|Ad", colnames(comb1))]
comb.comb<-comb1
# names(comb.comb)
# dim(comb.comb)

#Scale bioclimatic variables
for(i in c(29,30,33,34,35,36,37,38,39)){
comb.comb[,i] <- (comb.comb[,i]/10)
}

for(j in c(31,32)){
comb.comb[,j] <- (comb.comb[,j]/100)
}
#head(comb.comb[,c(29:39)])

#########################################################################
#Non-parametric Kruskas-Wallis Test
#########################################################################

names(comb.comb)
formulae2 <- lapply(colnames(comb.comb)[8:ncol(comb.comb)], function(x) as.formula(paste0(x, " ~ Cluster")))

res.kw <- lapply(formulae2, function(x) kruskal.test(x, data = comb.comb))
names(res.kw) <- format(formulae2)

pvals.kw <- lapply(formulae2, function(x) kruskal.test(x, data = comb.comb)[["p.value"]])
names(pvals.kw) <- format(formulae2)

#res.kw
#pvals.kw
#length(pvals.kw)

###Bonferroni correction
alpha=0.05

thr<-alpha/length(pvals.kw)

hits.kw<-which(pvals.kw<thr)
length(pvals.kw[hits.kw])
sig_kw <- t(as.data.frame(pvals.kw[hits.kw]))

#sig_kw

write.csv(sig_kw, file="metadata/Kruskas.Wallis.phenotypes.after.Bonferroni.csv", row.names=TRUE)

############################################################################################
#export data 
############################################################################################

save.image(file="metadata/pheno.bio.analysis.step.1.Rdata")
#load("pheno.bio.analysis.step.1.Rdata")

#newnames
#colnames(bioclim)<-newnames
pheno<-comb.comb[,c(1,13:ncol(comb.comb))]
pheno[1:10, 1:25]
names(comb.comb)

dim(pheno)
#drop single environment data
pheno <- pheno[,c(-grep(".Akaki", colnames(pheno)))]
pheno <- pheno[,c(-grep(".Adet", colnames(pheno)))]
pheno <- pheno[,c(-grep("OrderNorm", colnames(pheno)))]
#drop combined phenotypes with low HEb2
#he.thr=0.4

####Write table for GWAS
pheno$ID <- str_replace(pheno$ID, "T_", "T") #ID without underscore to match vcf to perform GWAS

write.table(pheno, file = "phenotypes/phenotypes_all.txt", quote = FALSE, sep="\t", row.names = FALSE)
save.image(file="metadata/metadata.pheno.analysis.Rdata")
