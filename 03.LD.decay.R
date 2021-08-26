########################################################################################
# This script perform LD analyses 
# Date: 11/24/2020
# Authors: Matteo Dell'Acqua  and Leonardo Caproni
########################################################################################

#load packages, set working directory, import .vcf and .csv files and create genind object and data frame based on first script (Fig2a_SI7_Phylogenetic_tree.R)

rm(list=ls())
maindir<-"/glab/projects/teff"
LDdir <-"/glab/projects/teff/LD"

setwd(LDdir)
options(stringsAsFactors = F)

library(tidyverse)
library(data.table)
library(genetics)
library(LDheatmap)
library(zoo)
library(RColorBrewer)

#A reduced SNP subset excluding SNPs of minor allele frequency <5% is used
#convert VCF t0 hmp and filter to MAF 5%
#Import hapmap file 
data<-read.delim("teff.snps.noWild.FINAL.AF005.hmp.txt")
#trim hapmap
data<-data[,-5:-11]
rownames(data)<-data[,1]
data<-data[,-1:-2]

#revert to mbp
data[,2]<-data[,2]/1000000

#converd to matrix for faster substitutions
dat1<-as.matrix(data)
# dat2 <- dat1[,-c(1:3)]
# dat2[1:10,1:10]

dat1[1:10,1:10]
dat.chr.pos <- dat1[,1:2]
#substitue characters
dat1<-sub("A", "A/A", dat1[,3:ncol(dat1)])
dat1<-sub("C", "C/C", dat1[,3:ncol(dat1)])
dat1<-sub("G", "G/G", dat1[,3:ncol(dat1)])
dat1<-sub("T", "T/T", dat1[,3:ncol(dat1)])
dat1<-sub("R", "A/G", dat1[,3:ncol(dat1)])
dat1<-sub("Y", "C/T", dat1[,3:ncol(dat1)])
dat1<-sub("S", "G/C", dat1[,3:ncol(dat1)])
dat1<-sub("W", "A/T", dat1[,3:ncol(dat1)])
dat1<-sub("K", "G/T", dat1[,3:ncol(dat1)])
dat1<-sub("M", "A/C", dat1[,3:ncol(dat1)])
dat1<-sub("N", "<NA>", dat1[,3:ncol(dat1)])

dat2 <- cbind(dat.chr.pos,dat1)
dat2[1:10,1:10]

data<-data.frame(dat2) #revert to data.frame
data$chrom <- str_replace(data$chrom, "LCL\\|", "")

###################################################
# adjust chromosome names
###################################################
data$chrom[data$chrom == "1A"] <- 1
data$chrom[data$chrom == "1B"] <- 2
data$chrom[data$chrom == "2A"] <- 3
data$chrom[data$chrom == "2B"] <- 4
data$chrom[data$chrom == "3A"] <- 5
data$chrom[data$chrom == "3B"] <- 6
data$chrom[data$chrom == "4A"] <- 7
data$chrom[data$chrom == "4B"] <- 8
data$chrom[data$chrom == "5A"] <- 9
data$chrom[data$chrom == "5B"] <- 10
data$chrom[data$chrom == "6A"] <- 11
data$chrom[data$chrom == "6B"] <- 12
data$chrom[data$chrom == "7A"] <- 13
data$chrom[data$chrom == "7B"] <- 14
data$chrom[data$chrom == "8A"] <- 15
data$chrom[data$chrom == "8B"] <- 16
data$chrom[data$chrom == "9A"] <- 17
data$chrom[data$chrom == "9B"] <- 18
data$chrom[data$chrom == "10A"] <- 19
data$chrom[data$chrom == "10B"] <- 20

data<-data[order(data$chrom, data$pos),]

data[1:10,1:7]

data<-t(data) #transposing

data[1:10,1:7]


gen<-makeGenotypes(data) #Make genotypes. 
save(gen, file="../metadata/LDgenotypes.Rdata")

#remove snp position
snpspos<-data.frame(t(data[1:2,]))
snpspos[,1]<-as.numeric(as.character(snpspos[,1]))##################################################################
snpspos[,2]<-as.numeric(as.character(snpspos[,2]))

save(snpspos, file="../metadata/snpspos.Rdata")

setwd(LDdir)
#create heatmaps of pairwise LD for all SNPs, iterating per chromosome. Long process to run on a personal computer
snpchr<-split(snpspos, snpspos[,1])

for (c in 1:length(snpchr)){
  print(c)
  tmp<-gen[,colnames(gen) %in% rownames(snpchr[[c]])] #extract from the data frame the SNP data from each chromosome
  #remove all non-genotype columns for now
  zz<-lapply(tmp, class)
  notgen<-grep("character", zz)
  if(length(notgen)>0){
    tmp<-tmp[,-notgen]
  }
  #remove also from snpchr
  snpchr[[c]]<-snpchr[[c]][rownames(snpchr[[c]]) %in% colnames(tmp),]
  print(paste("SNP removed", length(notgen), "out of", ncol(tmp)))
  pdf(paste("CHR", c, "LDheatmap.pdf", sep = "."), width = 9, height = 7)
  tmp<-LDheatmap(tmp, genetic.distances=snpchr[[c]][,2], distance="physical",color=heat.colors(20))
  dev.off()
  save(tmp, file=paste("CHR", c, "LDheatmap.Rdata", sep = "."))
}

######################
#restart from here
setwd(LDdir)
options(stringsAsFactors = F)

#print out easier to handle files
load("../metadata/snpspos.RData")

for (i in 1:20){
  print(i)
  load(paste("CHR", i, "LDheatmap.Rdata", sep = "."))
  tmpld<-as.matrix(tmp$LDmatrix)
  
  #keep positions and rearrange the dataframe / create a position data frame
  tmpos<-snpspos[rownames(snpspos) %in% colnames(tmpld),]
  tmpos$marker<-rownames(tmpos)
  colnames(tmpos)<-c("chr","position","marker")
  
  #collapse the LD matrix
  m<-cbind(which(!is.na(tmpld),arr.ind = TRUE),na.omit(as.vector(tmpld)))
  colnames(m)<-c("m1","m2","r2")
  tmpos<-tmpos[order(tmpos[,2]),]
  tmpos$idx<-1:nrow(tmpos)
  
  #add position information. For each pairwise comparison this will yield marker names, positions and r2
  newpos<-merge(m, tmpos, by.x="m1", by.y="idx", all.x=T)
  newpos1<-merge(newpos, tmpos, by.x="m2", by.y="idx", all.x=T)
  #drop chromosome data
  newpos1<-subset(newpos1, select=-c(chr.x,chr.y))
  
  
  dist<-abs(newpos1$position.x-newpos1$position.y) #calculate distance as the difference in Mbp between two markers
  newpos1<-cbind(newpos1,dist) #add distance vector to matrix
  #rearrange columns in matrix
  m<-newpos1[,c("m1","position.x","m2","position.y","r2","dist")]
  #order based on dist
  m<-m[order(m[,"dist"]),]
  save(m, file=paste("collapsed_matrix_LD_chr",i,"Rdata", sep="."))
  
}

#set interpolation function based on Hill and Weir (1988) equation based on Marroni (2011) script available at https://fabiomarroni.wordpress.com/2011/08/09/estimate-decay-of-linkage-disequilibrium-with-distance/
Er<-function(C_,d){
  length(d)
  res<-((10+C_*d)/((2+C_*d)*(11+C_*d)))*(1+((3+C_*d)*(12+12*C_*d+(C_*d)^2))/(n*(2+C_*d)*(11+C_*d)))
  return(res)
}

#set additional parameters
n=366
exp <- expression(italic(paste(displaystyle(r^2))))

flist<-list()
datashd<-list()

#obtain fpoints by chromosomes
for (j in 1:20){
  print(j)
  load(paste("collapsed_matrix_LD_chr",j,"Rdata", sep="."))
  m1<-m
  m1<-m1[order(m1$dist),]
  #drop variables that will not be used
  m1<-m1[,c("r2","dist")]
  ld<-m1[,"r2"]
  d<-m1[,"dist"]
  nlm<-nls(ld~Er(C_,d[order(d)]),start=c(C_=0.1),control=nls.control(maxiter=100))
  C_<-summary(nlm)$parameters[1]
  fpoints<-Er(C_,d[order(d)])
  tmp<-data.frame(d,fpoints,j)
  flist[[j]]<-data.frame(d,fpoints,j)
  ldd<-0.1
  xpos<-tmp[which((tmp[,2]-ldd)<=0)[1],1]
  datashd[[j]]<-c(max(tmp[,2]), xpos)
}

fpointsall<-do.call(rbind, flist)
# fpointsall[,3]<-as.factor(fpointsall[,3])
datashd<-do.call(rbind, datashd)

colnames(fpointsall)<-c("distance","r2","Chromosome")

fpointsall$distanceKb<-fpointsall$distance*1000

#Export maximum LD and LD decay distance per chromosome to .csv file
colnames(datashd)<-c("max", "Mb")
datashd<- as.data.frame(datashd)
datashd$Kb = datashd$Mb * 1000
write.csv(datashd,"SupplementaryInformation.SX.csv")
datashd[1:10,1:3]
#now upgrade the basic plots using sliding windows

#set parameters:
#pairwise r2 is averaged for all surrounding markers within a +/-'win' of each SNP. In the paper, the 'win' parameter is set to 5 Mbp.
win<-5 
#resulting LD are plotted against physical position, averaging values over a sliding window of 'swin'% of each chromosome's markers. In the paper, 'swin' parameter is set to 5%.
swin<-floor(length(unique(m[,1]))*.05) 

#import centromeres' positions
plotter<-list()
centromeres<-read.csv("centromeres.csv")
centromeres$chromosome<-as.factor(centromeres$chromosome)
colnames(centromeres)[3]<-"Chromosome"

#estimate average LD per SNP
for (i in 1:20){
  load(paste("collapsed_matrix_LD_chr",i,"Rdata", sep="."))
  #remove long distances (the interval is on each side of the current SNP)
  m1<-m[which(m$dist<win),] 
  #sort by pos
  m1<-m1[order(m1$position.x),]
  #binding x and y positions
  m2<-m1[,-c(1:2)]
  colnames(m2)[1:2]<-c("m","position")
  m1<-m1[,-c(3:4)]
  colnames(m1)[1:2]<-c("m","position")
  m1<-rbind(m1,m2)
  #get mean for each unique SNP
  m1<-aggregate(m1['r2'], by=m1['position'], mean) #mean of r2 by position for SNPs
  #go through a sliding window
  #noting SNP positions
  md<-rollapply(m1$r2,width=swin, function(x) mean(x), partial=T) #mean of r2 by position for SNPs that are close by
  #create df for plotting
  md<-cbind(m1$position, md, i)
  plot(md[,1],md[,2], type="l", xlab="Mb", ylim=c(0,0.5))
  points(x=centromeres[i,1], y=0, col="red")
  plotter[[i]]<-md
}

#organize data frame
plotter<-as.data.frame(do.call(rbind, plotter))
colnames(plotter)<-c("position","mean_r2","Chromosome")

#change chr names again to plot decay along chrs
plotter$Chromosome[plotter$Chromosome == "1"] <- "1A"
plotter$Chromosome[plotter$Chromosome == "2"] <- "1B"
plotter$Chromosome[plotter$Chromosome == "3"] <- "2A"
plotter$Chromosome[plotter$Chromosome == "4"] <- "2B"
plotter$Chromosome[plotter$Chromosome == "5"] <- "3A"
plotter$Chromosome[plotter$Chromosome == "6"] <- "3B"
plotter$Chromosome[plotter$Chromosome == "7"] <- "4A"
plotter$Chromosome[plotter$Chromosome == "8"] <- "4B"
plotter$Chromosome[plotter$Chromosome == "9"] <- "5A"
plotter$Chromosome[plotter$Chromosome == "10"] <- "5B"
plotter$Chromosome[plotter$Chromosome == "11"] <- "6A"
plotter$Chromosome[plotter$Chromosome == "12"] <- "6B"
plotter$Chromosome[plotter$Chromosome == "13"] <- "7A"
plotter$Chromosome[plotter$Chromosome == "14"] <- "7B"
plotter$Chromosome[plotter$Chromosome == "15"] <- "8A"
plotter$Chromosome[plotter$Chromosome == "16"] <- "8B"
plotter$Chromosome[plotter$Chromosome == "17"] <- "9A"
plotter$Chromosome[plotter$Chromosome == "18"] <- "9B"
plotter$Chromosome[plotter$Chromosome == "19"] <- "10A"
plotter$Chromosome[plotter$Chromosome == "20"] <- "10B"

#change chr names again to plot decay
fpointsall$Chromosome[fpointsall$Chromosome == "1"] <- "1A"
fpointsall$Chromosome[fpointsall$Chromosome == "2"] <- "1B"
fpointsall$Chromosome[fpointsall$Chromosome == "3"] <- "2A"
fpointsall$Chromosome[fpointsall$Chromosome == "4"] <- "2B"
fpointsall$Chromosome[fpointsall$Chromosome == "5"] <- "3A"
fpointsall$Chromosome[fpointsall$Chromosome == "6"] <- "3B"
fpointsall$Chromosome[fpointsall$Chromosome == "7"] <- "4A"
fpointsall$Chromosome[fpointsall$Chromosome == "8"] <- "4B"
fpointsall$Chromosome[fpointsall$Chromosome == "9"] <- "5A"
fpointsall$Chromosome[fpointsall$Chromosome == "10"] <- "5B"
fpointsall$Chromosome[fpointsall$Chromosome == "11"] <- "6A"
fpointsall$Chromosome[fpointsall$Chromosome == "12"] <- "6B"
fpointsall$Chromosome[fpointsall$Chromosome == "13"] <- "7A"
fpointsall$Chromosome[fpointsall$Chromosome == "14"] <- "7B"
fpointsall$Chromosome[fpointsall$Chromosome == "15"] <- "8A"
fpointsall$Chromosome[fpointsall$Chromosome == "16"] <- "8B"
fpointsall$Chromosome[fpointsall$Chromosome == "17"] <- "9A"
fpointsall$Chromosome[fpointsall$Chromosome == "18"] <- "9B"
fpointsall$Chromosome[fpointsall$Chromosome == "19"] <- "10A"
fpointsall$Chromosome[fpointsall$Chromosome == "20"] <- "10B"

## now convert Chrs as factors
# plotter[,3]<-as.factor(plotter[,3])
# fpointsall[,3]<-as.factor(fpointsall[,3])


#Define graph parameters to represent each chromosome for both subgenomes
cols<-c("#a55bd1",
        "#72be46",
        "#cf4da3",
        "#4a8a35",
        "#596dd0",
        "#d0a733",
        "#8a539b",
        "#60c185",
        "#d83f64",
        "#43c2c6",
        "#d95633",
        "#5e90cb",
        "#ce8a4e",
        "#c18fd9",
        "#a8ae59",
        "#9e4867",
        "#388661",
        "#e2849b",
        "#756e29",
        "#ad563c")

## Define graph parameters to represent each chromosome pairs (homeologus)
cols10 <- c("#b74454",
            "#a8db57",
            "#a9a240",
            "#6261b8",
            "#61ca83",
            "#e7553f",
            "#4e8d38",
            "#b255a1",
            "#d88f47",
            "#92c8a8")

datashd$cols<-cols

#Export to .pdf file
library(ggplot2)

#split subgenomes to plot indipendently A and B
fpointsA <- filter(fpointsall, grepl("A", fpointsall$Chromosome, fixed = FALSE))
fpointsB <- filter(fpointsall, grepl("B", fpointsall$Chromosome, fixed = FALSE))

plotterA  <- filter(plotter, grepl("A", plotter$Chromosome, fixed = FALSE))
plotterB  <- filter(plotter, grepl("B", plotter$Chromosome, fixed = FALSE))

centromeresA <- filter(centromeres, grepl("A", centromeres$Chromosome, fixed = FALSE))
centromeresB <- filter(centromeres, grepl("B", centromeres$Chromosome, fixed = FALSE))


decay <- ggplot(fpointsall) + 
        geom_line(aes(y=r2, x=distance, col=fct_relevel(Chromosome,'1A','1B','2A','2B','3A','3B','4A','4B','5A','5B','6A','6B','7A','7B','8A','8B','9A','9B','10A','10B')),size=1) + 
        labs(col = "Chromosomes") +
        scale_color_manual(values=cols) +
        xlab("Mb") + ylab(expression(italic(r)*{}^2)) +
        coord_cartesian(ylim=c(0.0,0.40),xlim=c(0,3)) + #set x and y axes limit
        theme_classic(base_size = 14)+
        theme(legend.position="none")

decayA <- ggplot(fpointsA) + 
  geom_line(aes(y=r2, x=distance, col=fct_relevel(Chromosome,'1A','2A','3A','4A','5A','6A','7A','8A','9A','10A')),size=1) + 
  labs(col = "Chromosomes") +
  scale_color_manual(values=cols10) +
  xlab("Mb") + ylab(expression(italic(r)*{}^2)) +
  coord_cartesian(ylim=c(0.0,0.40),xlim=c(0,3)) + #set x and y axes limit
  theme_classic(base_size = 14) +
  theme(legend.position="none")

decayB <- ggplot(fpointsB) + 
  geom_line(aes(y=r2, x=distance, col=fct_relevel(Chromosome,'1B','2B','3B','4B','5B','6B','7B','8B','9B','10B')),size=1) + 
  labs(col = "Chromosomes") +
  scale_color_manual(values=cols10) +
  xlab("Mb") + ylab(expression(italic(r)*{}^2)) +
  coord_cartesian(ylim=c(0.0,0.40),xlim=c(0,3)) + #set x and y axes limit
  theme_classic(base_size = 14) +
  theme(legend.position="none")



############
LD.along.chrs<- ggplot(plotter, aes(x=position)) + 
        geom_line(aes(y=mean_r2, col=fct_relevel(Chromosome,'1A','1B','2A','2B','3A','3B','4A','4B','5A','5B','6A','6B','7A','7B','8A','8B','9A','9B','10A','10B')), 
                  size=0.7 ) + labs(col = "Chromosomes") +
        geom_point(data=centromeres, mapping=aes(x=x.position, y=y.position),shape=17, color="black", size=1.5)+
        xlab("Mb") + ylab(expression(italic(r)*{}^2)) +
        scale_color_manual(values=cols) +
        theme_classic(base_size = 12) +
        theme(legend.justification=c(1,0), legend.position="right", strip.text.y = element_blank()) +
        scale_x_continuous(breaks=seq(0,320,20), sec.axis = dup_axis()) +
        facet_grid(fct_relevel(Chromosome,'1A','1B','2A','2B','3A','3B','4A','4B','5A','5B','6A','6B','7A','7B','8A','8B','9A','9B','10A','10B') ~ .)
        
#dev.off()

LD.along.chrs.A<- ggplot(plotterA, aes(x=position)) + 
  geom_line(aes(y=mean_r2, col=fct_relevel(Chromosome,'1A','2A','3A','4A','5A','6A','7A','8A','9A','10A')),size=0.7) +
  labs(col = "Chromosomes") +
  geom_point(data=centromeresA, mapping=aes(x=x.position, y=y.position),shape=17, color="black", size=1.5)+
  xlab("Mb") + ylab(expression(italic(r)*{}^2)) +
  scale_color_manual(values=cols10) +
  theme_classic(base_size = 12) +
  theme(legend.justification=c(1,0), legend.position="right", strip.text.y = element_blank()) +
  scale_x_continuous(breaks=seq(0,320,20), sec.axis = dup_axis()) +
  facet_grid(fct_relevel(Chromosome,'1A','2A','3A','4A','5A','6A','7A','8A','9A','10A') ~ .)
  
LD.along.chrs.B<- ggplot(plotterB, aes(x=position)) + 
  geom_line(aes(y=mean_r2, col=fct_relevel(Chromosome,'1B','2B','3B','4B','5B','6B','7B','8B','9B','10B')),size=0.7) +
  labs(col = "Chromosomes") +
  geom_point(data=centromeresB, mapping=aes(x=x.position, y=y.position),shape=17, color="black", size=1.5)+
  xlab("Mb") + ylab(expression(italic(r)*{}^2)) +
  scale_color_manual(values=cols10) +
  theme_classic(base_size = 12) +
  theme(legend.justification=c(1,0), legend.position="right", strip.text.y = element_blank()) +
  scale_x_continuous(breaks=seq(0,320,20), sec.axis = dup_axis()) +
  facet_grid(fct_relevel(Chromosome,'1B','2B','3B','4B','5B','6B','7B','8B','9B','10B') ~ .)

############# Patchwork
# library(patchwork)
# LD.sub.A <- LD.along.chrs.A + inset_element(decayA,left = 0.5, bottom = 0, right = 1, top = 0.5, on_top = TRUE, align_to = 'full')


############ Save images
ggsave(LD.along.chrs, filename = "../figures/LD_chrs_teff.pdf",height=13, width =8)
ggsave(LD.along.chrs.A, filename = "../figures/LD_chrs_teff.subgenome.A.pdf",height=8, width =6)
ggsave(LD.along.chrs.B, filename = "../figures/LD_chrs_teff.subgenome.B.pdf",height=8, width =6)

ggsave(decay, filename = "../figures/LD_decay_teff.pdf", height=5, width =5)
ggsave(decayA, filename = "../figures/LD_decay_teff.subgenome.A.pdf", height=5, width =5)
ggsave(decayB, filename = "../figures/LD_decay_teff.subgenome.B.pdf", height=5, width =5)


##################################################
# Now prepare data for candidate gene analysis
##################################################

# split fpoints by chromosome
bychr<-split(fpointsall, fpointsall[,3])

# set a r2 thrshold 
thr<-0.3

# find position closer to r2 = thr
tmp<-bychr[[1]]
idx<-which(tmp$r2 <= thr)[1]
outdf<-tmp[idx,]

for (c in 2:length(bychr)){ #c=1
  tmp<-bychr[[c]]
  idx<-which(tmp$r2 <= thr)[1]
  outdf<-rbind(outdf, tmp[idx,])
}#for c

decay.window.bychr.03<-outdf [,-c(1,2)]
write.csv(decay.window.bychr.03, file="decay.window.by.chr.03.csv", row.names=F, quote=F)

#Let's do it again
# change to a less stringent threshold (in terms of window size)
thr<-0.2

tmp<-bychr[[1]]
idx<-which(tmp$r2 <= thr)[1]
outdf<-tmp[idx,]

for (c in 2:length(bychr)){ #c=1
  tmp<-bychr[[c]]
  idx<-which(tmp$r2 <= thr)[1]
  outdf<-rbind(outdf, tmp[idx,])
}#for c

decay.window.bychr.02<-outdf [,-c(1,2)]
write.csv(decay.window.bychr.03, file="decay.window.by.chr.02.csv", row.names=F, quote=F)

save(decay.window.bychr.03, decay.window.bychr.02, file="../metadata/ouput.LD.analysis.Rdata")
