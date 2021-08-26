#!/usr/bin/R

#GWAS
rm(list=ls())
options(stringsAsFactors = F)
maindir <- "/glab/projects/teff"
setwd(maindir)

library(rMVP)
library(stringr)
library(dplyr)
library(tidyverse)
library(data.table)
library(ape)
library("qvalue")

phe <-"phenotypes/phenotypes_all.txt"
vcf <- "/glab/projects/teff/snps_HaplotypeCaller/no_wild/teff.snps.noWild.FINAL.AF001.vcf" 
vcf.pruned <- "/glab/projects/teff/snps_HaplotypeCaller/no_wild/plink/pruned_100_10_03/pruned_100_10_03.vcf" 

# # make an output GWAS directory
# GWAS <-"GWAS"
# dir.create(GWAS)

#2.1) Prepare MVP datasets for full SNPs data
MVP.Data(fileVCF=vcf, 
         filePhe=phe, 
         sep.phe="\t", 
         fileKin=F, 
         filePC=F,
         out="GWAS/mvp.out")
genotype <- attach.big.matrix("GWAS/mvp.out.geno.desc")
phenotype <- read.table("GWAS/mvp.out.phe",head=TRUE)
map <- read.table("GWAS/mvp.out.geno.map" , head = TRUE)
#K <- MVP.K.VanRaden(genotype)
#CV <- MVP.PCA(genotype, pcs.keep = 30)
#map[1:5,1:5]

map$CHROM <- str_replace(map$CHROM, "lcl\\|", "") # ouble \\ to interpret | symbol literally

chromOrder = c("1A", "1B","2A","2B","3A","3B","4A","4B",
               "5A","5B","6A","6B","7A","7B","8A","8B",
               "9A","9B","10A","10B","U")

map <- map[order(match(map$CHROM, chromOrder)),]


#2.2) Prepare MVP datasets for pruned SNPs data.
MVP.Data(fileVCF=vcf.pruned, 
         filePhe=phe, sep.phe="\t", 
         fileKin=FALSE, 
         filePC=FALSE, 
         out="GWAS/mvp.pruned")

genotype_pruned <- attach.big.matrix("GWAS/mvp.pruned.geno.desc")
K.p <- MVP.K.VanRaden(genotype_pruned)
CV.p <- MVP.PCA(genotype_pruned, pcs.keep = 30)
#MVP.PCAplot(CV.p, plot3D = F)


results.bon <-"GWAS/results.bonferroni"
dir.create(results.bon)
setwd(results.bon)

#Define threshold
n.markers <- 12153 #n of markers for GWAS
alpha <- 0.05 #Bonferroni alpha
hap.blocks <- 2100 #n of non-associated markers r2=0.3 in the set
thr <- (alpha/hap.blocks)*n.markers #0.512

#set a vector for errors
errs<-c()

#Run for loop MVP
for (j in c(10,15)){ 
  
  dirname<-paste0("PC.", j)
  dir.create(paste0("./",dirname))
  setwd(dirname)
  
  for(i in c(2:ncol(phenotype))){    ##(2:5,50:52,68:70,116:118,161:165,197:201
    tryCatch({ #this function is used to handle errors (in any)
      tmpMVP <- 
        MVP(
          phe=phenotype[,c(1,i)], 
          geno=genotype,
          map=map,
          K=K.p,
          #CV.GLM=CV.p[,c(1:j)],
          #CV.MLM=CV.p[,c(1:j)],
          CV.FarmCPU=CV.p[,c(1:j)],
          #nPC.FarmCPU=2, # "If pcs have been added in covariate files, PLEASE DO NOT assign value to nPC.GLM, nPC.MLM, nPC.FarmCPU"
          #nPC.GLM=2,
          #nPC.MLM=2,
          #perc=1,
          #priority="speed",
          ncpus=3,
          vc.method="EMMA",
          maxLoop=10,
          method.bin="FaST-LMM",
          permutation.threshold=F,
          file.output=T,
          #permutation.rep=500,
          threshold=thr,
          method=c("FarmCPU"),
          file="jpg",
          dpi=300,
          col = c("#667fce",
                  "#4faa6f") ### M plot colors
        )
      save(tmpMVP, file=paste0(colnames(phenotype)[i],".rMVP.bonferroni.Rdata"))
      gc()
      
    }, error=function(e){errs<-c(errs, i)}) #trycatch ends #for i loop ends
  }#for i
  #go back one level
  setwd("..")
} # for jloop ends



#3.3) Add FDR threshold
wd_fdr <- "/glab/projects/teff/GWAS/results.bonferroni/PC.10"
setwd(wd_fdr)

files=list.files(path= wd_fdr, pattern="*FarmCPU.csv", full.names = FALSE) 
pattern="*FarmCPU.csv"
pat=".FarmCPU.csv"

outlist <- list()
minsig<-list()
thrfdr<-0.05

for(i in 1:length(files)){ #i=1
  tmp <- read.csv(files[i], header=T)
  
  colnames(tmp)[8]<-"pvalue"
  tmp$qvalue<-qvalue(tmp[,8])$qvalue

  hit<-tmp[which(tmp$qvalue<thrfdr),]
  
  
  if(nrow(hit)>0){
    #store mininym sig
    minsig[[i]]<-max(hit$pvalue)
    names(minsig)[i]<-files[i]

    out<-data.frame(trait=files[i], hit)
    outlist[[i]] <- out
  }#if hit
}

minsig
minsig_df <- do.call(rbind, minsig)

sigFDR<-do.call(rbind, outlist)
#dim(sigFDR)
#table(sigFDR[,1])
