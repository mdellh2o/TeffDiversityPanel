########################################################################################
# This script performs plotting for penotypic analysis
# onto the geographic space
# Date: 11/24/2020
# Authors: Matteo Dell'Acqua 
########################################################################################

# load libraries
options(stringsAsFactors = F)
library(corrplot)
library(ggplot2)
library(vioplot)
library(asreml)
library(reshape2)
library(ggcorrplot)
library(patchwork)
library(tidyr)
library(ggExtra)
library(tidyr)
library(ggfortify)
library(ggpubr)
library(ggsci)
#devtools::install_github("jfq3/QsRutils", build_vignettes = TRUE)
library(QsRutils)
library(ggridges)
library(ggalluvial)
library(dplyr)
library(easyalluvial)
library(adegenet)
library(ggsci)
library(UpSetR)
library(ggVennDiagram)

wd <- "/glab/projects/teff"
setwd(wd)

# set colors by gender
colmen<-"#0093dc"
colwom<-"#ac1084"

#get info about varieties
load("metadata/diversity.varieties.Rdata")

#get in metadata about farmers
load("metadata/metadata.diversity.analysis.Rdata")
load("metadata/metadata.pheno.analysis.Rdata")
names(comb.comb)
comb.comb[1:5,1:10]
rownames(comb.comb)<-comb.comb[,1]
colnames(comb.comb)<-sub("PL.y", "PL", colnames(comb.comb))

quants<-colnames(comb.comb)[10:28]
bios<-colnames(comb.comb)[29:52]
traits<-colnames(comb.comb)[53:ncol(comb.comb)]

#remove color appreciation (CA from the phenotypic dataset)
traits<-traits[-grep("^CA", traits)]

#make plots by traits
for (t in traits){ ## t = traits[1]
  tmp<-comb.comb[,c(1, grep(t, colnames(comb.comb)))]
  png(paste0("./traits/", t, ".BLUP.vioplot.png"))
  vioplot(tmp[,2:ncol(tmp)], col='lightgrey', ylab=t, las=2)
  dev.off()
}


#keep only phenotypes and remove normalized phenotypes for the moment
tmp.comb.comb<-comb.comb[,traits]
tmp.comb.comb<-tmp.comb.comb[,-grep("OrderNorm", colnames(-tmp.comb.comb))]
head(tmp.comb.comb)

#correlate among locations
adblup<-tmp.comb.comb[,grep("Adet", colnames(tmp.comb.comb))]
akblup<-tmp.comb.comb[,grep("Akaki", colnames(tmp.comb.comb))]
cradak<-cor(adblup, akblup, use="complete")

akadplot<-ggcorrplot(cradak, title = "Akaki VS Adet")
akadplot
ggsave(akadplot, filename="figures/BLUPs.correlation.by.location.png")

#correlate among  genders
mblup<-tmp.comb.comb[,grep("\\.W$", colnames(tmp.comb.comb))]
wblup<-tmp.comb.comb[,grep("\\.M$", colnames(tmp.comb.comb))]
crmw<-cor(mblup, wblup, use="complete")

mwplot<-ggcorrplot(crmw, title = "Men VS Women") +theme_light()
mwplot
ggsave(mwplot, filename="figures/farmer.BLUPs.correlation.by.gender.png")

#make a general correlation plot
cr1<-cor(tmp.comb.comb, use="complete")
totcorr<-ggcorrplot(cr1, tl.cex=10, tl.srt = 90)

ggsave(totcorr, filename="figures/BLUPs.all.correlation.pdf", height=10, width=10)
ggsave(totcorr, filename="figures/BLUPs.all.correlation.png", height=10, width=10)


#################################
#make a correlation plot on combined values
######################

combined<-tmp.comb.comb[,-grep("Ade|Aka", colnames(tmp.comb.comb))]
head(combined)

cr<-cor(y=combined[,-grep("^OA|^PA", colnames(combined))], x=combined[,grep("^OA$|^PA$", colnames(combined))], use="complete")

allcorplot<-ggcorrplot(cr, outline.col = "gray")
#allcorplot<- allcorplot + geom_hline(yintercept=16.5, linetype="dashed") + geom_vline(xintercept=16.5, linetype="dashed")
allcorplot <- allcorplot + #labs(tag="A") +  
            theme_light() + xlab("") +ylab("")
allcorplot

ggsave(allcorplot, filename="figures/FIGURE.BLUPs.correlation.pdf", height=6, width=3)


#################################
#make scatterplot of farmers scores
######################

#keep only phenotypes and remove normalized phenotypes for the moment
farm.comb.comb<-tmp.comb.comb[,grep("^PA|^OA", colnames(tmp.comb.comb))]
head(farm.comb.comb)

p1<-ggplot(farm.comb.comb, aes(x=OA.Adet.M, y=OA.Adet.W)) + geom_point()+
  theme_light() + xlim(1, 5)+ ylim(1, 5)  + geom_abline(intercept = 0, slope = 1) + geom_point(color='darkred')+ coord_fixed() 
p1

p2<-ggplot(farm.comb.comb, aes(x=OA.Akaki.M, y=OA.Akaki.W)) + geom_point()+
  theme_light()+ xlim(1, 5)+ ylim(1, 5)  + geom_abline(intercept = 0, slope = 1)+ geom_point(color='navyblue')+ coord_fixed()
p2 

gendcombplot<-p1|p2
gendcombplot

ggsave(gendcombplot, file="figures/gender.correlation.plot.pdf", height=4, width=9)

#define ranking by farmers in both locations
perc<-.75
Adperc<-quantile(farm.comb.comb[,"OA.Akaki"], perc, na.rm=T)
Akperc<-quantile(farm.comb.comb[,"OA.Akaki"], perc, na.rm=T)

selAd<-rownames(farm.comb.comb)[which(farm.comb.comb[,"OA.Adet"]>Akperc)]
selAk<-rownames(farm.comb.comb)[which(farm.comb.comb[,"OA.Akaki"]>Adperc)]

length(selAd)
length(selAk)

which(selAd %in% selAk)

#make a pca on phentoypes
tokeep<-which(!names(tmp.comb.comb) %in% names(farm.comb.comb))
met.comb.comb<-tmp.comb.comb[,tokeep]
met.comb.comb<-met.comb.comb[,-grep("Ak|Ad", colnames(met.comb.comb))]
head(met.comb.comb)
traitpc<-prcomp(na.omit(met.comb.comb))
summary(traitpc)

#combine it with the general file
pcout<-traitpc$x[,1:5]
colnames(pcout)<-paste0(colnames(pcout), "_pheno")
head(pcout)

all.meta<-merge(comb.comb, pcout, by="row.names", all.x=T)

#reorganize metafile before plotting
all.meta<-all.meta[,-grep("Order", colnames(all.meta))]
colnames(all.meta)<-sub("Axis", "PC",colnames(all.meta))
rownames(all.meta)<-all.meta[,1]
all.meta<-all.meta[,-c(1:2)]
head(all.meta)

#do some plotting
evaltraits<-c("OA", "OA.Akaki", "OA.Adet", "OA.M", "OA.W")

toplot<-all.meta[,grep("PC[0-9]|Cluster", colnames(all.meta))]
toplot<-cbind(toplot, all.meta[,evaltraits])
head(toplot)

#remove samples with NAs
toplot<-na.omit(toplot)

#in each evalutation column, just put the top percentile
perc<-.9
for (i in evaltraits){ # i = "OA"
  tmp<-quantile(toplot[,i], perc)
  top<-rownames(toplot)[which(toplot[,i] >= tmp)]
  toplot[which(!rownames(toplot) %in% top),i]<-NA
}
head(toplot)


#plot farmer selection on genetic PCA
farmersel.geno<-ggplot(toplot, aes(x=PC1, y=PC2))+
  geom_point(alpha=0.9,size=1.5) + theme_light() + 
  #theme(legend.title = element_blank(), legend.position = "bottom", legend.text=element_text(size=8)) + 
  #scale_color_manual(values = c("gray", colmen, colwom),
  #                   limits = c("1", "2", "3"), 
  #                   labels = c("Teff genotype", "Men choice", "Women choice"))+ 
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")
farmersel.geno<-farmersel.geno + geom_point(data=toplot, aes(x=PC1, y=PC2), alpha=0.9,size=1.5, col = "gray")
farmersel.geno<-farmersel.geno + geom_point(data=subset(toplot, !is.na(toplot$OA.M)), aes(x=PC1, y=PC2), col = colmen, alpha=0.3, size=2)
farmersel.geno<-farmersel.geno + geom_point(data=subset(toplot, !is.na(toplot$OA.W)), aes(x=PC1, y=PC2), col = colwom, alpha=0.3, size=2)
farmersel.geno<-farmersel.geno + labs(title="Genotypic diversity", x="PC 1", y="PC 2")
farmersel.geno

#ggsave(farmersel.geno, file="PCA.farmers.selection.genetic.PCA.pdf", height=4, width=4)


#plot farmer selection on biolcim PCA
farmersel.bio<-ggplot(toplot, aes(x=PC1_bio, y=PC2_bio))+
  geom_point(alpha=0.9,size=1.5) + theme_light() + 
  #theme(legend.title = element_blank(), legend.position = "bottom", legend.text=element_text(size=8)) + 
  #scale_color_manual(values = c("gray", colmen, colwom),
  #                   limits = c("1", "2", "3"), 
  #                   labels = c("Teff genotype", "Men choice", "Women choice")) + 
geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")
farmersel.bio<-farmersel.bio + geom_point(data=toplot, aes(x=PC1_bio, y=PC2_bio), alpha=0.9,size=1.5, col = "gray")
farmersel.bio<-farmersel.bio + geom_point(data=subset(toplot, !is.na(toplot$OA.M)), aes(x=PC1_bio, y=PC2_bio), col = colmen, alpha=0.3, size=2)
farmersel.bio<-farmersel.bio + geom_point(data=subset(toplot, !is.na(toplot$OA.W)), aes(x=PC1_bio, y=PC2_bio), col = colwom, alpha=0.3, size=2)
farmersel.bio<-farmersel.bio + labs(title="Climatic diversity", x="PC 1", y="PC 2")
farmersel.bio

#ggsave(farmersel.bio, file="PCA.farmers.selection.bioclim.PCA.pdf", height=4, width=4)


#plot farmer selection on phenotypic PCA
farmersel.pheno<-ggplot(toplot, aes(x=PC1_pheno, y=PC2_pheno, col=Cluster))+
  geom_point(alpha=0.9,size=1.5) + theme_light() + 
  theme(legend.title = element_blank(), legend.position = "right", legend.text=element_text(size=8)) + 
  scale_color_manual(values = c("gray", colmen, colwom),
                     limits = c("1", "2", "3"), 
                     labels = c("Teff genotype", "Men choice", "Women choice")) + 
  geom_hline(yintercept=0, linetype="dashed") + geom_vline(xintercept=0, linetype="dashed")
farmersel.pheno<-farmersel.pheno + geom_point(data=toplot, aes(x=PC1_pheno, y=PC2_pheno), alpha=0.9,size=1.5, col = "gray")
farmersel.pheno<-farmersel.pheno + geom_point(data=subset(toplot, !is.na(toplot$OA.M)), aes(x=PC1_pheno, y=PC2_pheno), col = colmen, alpha=0.3, size=2)
farmersel.pheno<-farmersel.pheno + geom_point(data=subset(toplot, !is.na(toplot$OA.W)), aes(x=PC1_pheno, y=PC2_pheno), col = colwom, alpha=0.3, size=2)
farmersel.pheno<-farmersel.pheno + labs(title="Phenotypic diversity", x="PC 1", y="PC 2")
farmersel.pheno

#ggsave(farmersel.pheno, file="PCA.farmers.selection.phenotypes.PCA.pdf", height=4, width=4)

#make a composite plot
combotop<-farmersel.geno|farmersel.pheno
combotop
ggsave(combotop, file="figures/PCA.farmers.selection.combined.pdf", height=4, width=12)

#correlate PCs with traits
cortrait<-all.meta[, colnames(all.meta) %in% colnames(combined)]
cortrait<-cortrait[,-grep("W$|M$", colnames(cortrait))]
names(cortrait)

corpc<-all.meta[, grep("^PC[1-2]", colnames(all.meta))]
names(corpc)

cr <- cor(y=corpc, x=cortrait, use="complete")
corplot<-ggcorrplot(cr)  #+  labs(tag="C")
  ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8))
corplot

ggsave(corplot, file="figures/coorelation.PCs.traits.pdf", height=2, width=12)

##combo<-combotop/corplot
#combo
#ggsave(comboplot, file="PCA.farmers.selection.combined.pdf", height=6, width=12)

#############
# check GY and OA distribution among varieties and landraces
############
#get improved lines ids
varieties[,1]<-sub("T", "T_", varieties[,1])
imps<-varieties[which(varieties$Status == "breeding material"),1]
imps

all.meta$status<-"Landrace"
all.meta$status[which(rownames(all.meta) %in% imps)]<-"Breeding material"

#make a boxplot
oap<-ggplot(all.meta, aes(y=OA, x=as.factor(status)))+
  geom_boxplot(fill="gray") + theme_light() + labs(x="")+ guides(fill=F)

dhp<-ggplot(all.meta, aes(y=DM, x=as.factor(status)))+
  geom_boxplot(fill="gray") + theme_light() + labs(x="") + guides(fill=F)

lrvsimp<-oap|dhp
lrvsimp

###########################à
# make a bunch of boxplots
##########################
OA.box <- ggboxplot(all.meta, x = "Cluster", y="OA", color="Cluster") + 
  scale_color_nejm(alpha=1)+ 
  geom_jitter(data=all.meta, aes(color=Cluster, shape=status, size=status), alpha=0.6)+
  scale_shape_manual(values=c(15, 20))+
  scale_size_manual(values=c(2,1.2))+
  theme_light() +theme(legend.position = "none")
OA.box
box.rslt <- with(all.meta, graphics::boxplot(OA ~ Cluster, plot = FALSE))
ptt.rslt <- with(all.meta,pairwise.wilcox.test(OA, Cluster, p.adjust.method="bonferroni"))
ltrs <- make_letter_assignments(ptt.rslt)
x <- c(1:length(ltrs$Letters))
y <- box.rslt$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame(x, y, cbd)
lmts <- get_plot_limits(OA.box)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
#OA.box <- OA.box +geom_text(data = ltr_df, aes(x=x, y=y, label=cbd), nudge_y=y.nudge)
OA.box <- OA.box +geom_text(data = ltr_df, aes(x=x, y=5, label=cbd))
OA.box

pheno.box <- ggboxplot(all.meta, x = "Cluster", y="PL", color="Cluster") + 
  scale_color_nejm(alpha=1)+
  geom_jitter(data=all.meta, aes(color=Cluster, shape=status, size=status), alpha=0.6)+
  scale_shape_manual(values=c(15, 20))+
  scale_size_manual(values=c(2,1.2))+
  theme_light()# +theme(legend.position = "none")
pheno.box
box.rslt <- with(all.meta, graphics::boxplot(PL ~ Cluster, plot = FALSE))
ptt.rslt <- with(all.meta,pairwise.wilcox.test(PL, Cluster, p.adjust.method="bonferroni"))
ltrs <- make_letter_assignments(ptt.rslt)
x <- c(1:length(ltrs$Letters))
y <- box.rslt$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame(x, y, cbd)
lmts <- get_plot_limits(pheno.box)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
#pheno.box <- GFP.box +geom_text(data = ltr_df, aes(x=x, y=y, label=cbd), nudge_y=y.nudge)
pheno.box <- pheno.box +geom_text(data = ltr_df, aes(x=x, 42, label=cbd))
pheno.box

PC_bio.box <- ggboxplot(all.meta, x = "Cluster", y="PC2_bio", color="Cluster") + 
  scale_color_nejm(alpha=1)+
  geom_jitter(data=all.meta, aes(color=Cluster, shape=status, size=status), alpha=0.6)+
  scale_shape_manual(values=c(15, 20))+
  scale_size_manual(values=c(2,1.2))+
  theme_light()# +theme(legend.position = "none")
PC_bio.box
box.rslt <- with(all.meta, graphics::boxplot(PC2_bio ~ Cluster, plot = FALSE))
ptt.rslt <- with(all.meta,pairwise.wilcox.test(PC2_bio, Cluster, p.adjust.method="bonferroni"))
ltrs <- make_letter_assignments(ptt.rslt)
x <- c(1:length(ltrs$Letters))
y <- box.rslt$stats[5, ]
cbd <- ltrs$Letters
ltr_df <- data.frame(x, y, cbd)
lmts <- get_plot_limits(PC_bio.box)
y.range <- lmts$ymax - lmts$ymin
y.nudge <- 0.05 * y.range
#PC_bio.box <- PC_bio.box +geom_text(data = ltr_df, aes(x=x, y=y, label=cbd), nudge_y=y.nudge)
PC_bio.box <- PC_bio.box +geom_text(data = ltr_df, aes(x=x, y=10, label=cbd))
PC_bio.box

boxplots<-OA.box | pheno.box
boxplots
ggsave(boxplots, file="figures/pheno.by.cluster.pdf", height=4, width=12)

#make all the remaining boxplots
bxpltrt<-c("DH", "DM", "PH", "PL", "PBPM", "TT", "CLF", "CDF", "PW", "PY", "GY", 
           "BY", "HI", "GFP", "GFR", "BPR", "OA", "PA")
boxout<-list()
it<-1
  
for (i in bxpltrt){ #i = "DH"
  
  print(i)
  
  tmpbox <- ggboxplot(all.meta, x = "Cluster", y=i, color="Cluster") + 
    scale_color_nejm(alpha=1)+
    geom_jitter(data=all.meta, aes(color=Cluster, shape=status, size=status), alpha=0.6)+
    scale_shape_manual(values=c(15, 20))+
    scale_size_manual(values=c(2,1.2))+
    theme_light()# +theme(legend.position = "none")
  tmpbox
  
  box.rslt <- with(all.meta, graphics::boxplot(formula(paste0(i," ~ Cluster")), plot = FALSE))
  ptt.rslt <- pairwise.wilcox.test(all.meta[,i], all.meta[,"Cluster"], p.adjust.method="bonferroni")
  
  ltrs <- make_letter_assignments(ptt.rslt)
  x <- c(1:length(ltrs$Letters))
  y <- box.rslt$stats[5, ]
  cbd <- ltrs$Letters
  ltr_df <- data.frame(x, y, cbd)

  tmpboxout <- tmpbox + geom_text(data = ltr_df, aes(x=x, y=max(all.meta[,i], na.rm=T)*1.1, label=cbd))
  
  print(tmpboxout)
  
  boxout[[it]]<-tmpboxout
  
  ggsave(tmpboxout, file=paste0("./boxplots/boxplot.", i, ".pdf"), width=6, height=4)
  ggsave(tmpboxout, file=paste0("./boxplots/boxplot.", i, ".jpeg"), width=6, height=4)
  
  it<-it+1
}


###################################
# make alluvial 
###################################

alluvial<-all.meta[,52:ncol(all.meta)]
alluvial$ids<-rownames(all.meta)
alluvial<-na.omit(alluvial)

#remove breeding lines
#alluvial<-alluvial[alluvial$status == "Landrace",]

head(alluvial)

#select columns to work on
selcol<-c("OA", "OA.Akaki.W", "OA.Akaki.M",  "OA.Adet.W", "OA.Adet.M")
alluvialtmp<-alluvial[,selcol]
rng<-1:length(selcol)

#replace values with quantiles
alluvialq<-alluvialtmp
alluvialq[,rng]<-as.data.frame(sapply(alluvialq[,rng], function (x) {
  qx <- quantile(x)
  cut(x, qx, include.lowest = TRUE,
      labels = c("q1", "q2", "q3", "q4"))
}))

alluvialq[,rng]<-sapply(alluvialq[,rng], as.factor)
head(alluvialq)

alluvialq <- alluvialq %>% mutate_if(is.character,as.factor)
sapply(alluvialtmp, class)
sapply(alluvialq, class)

#reorder factors before plotting
#alluvialq[,rng]<-sapply(alluvialq[,rng], function(x) factor(x, levels=rev(levels(x))))

alluvialplot<-ggplot(data = alluvialq,
                aes(axis1 =OA.Adet.W , axis2= OA.Adet.M, axis3= OA.Akaki.M, axis4 = OA.Akaki.W  )) +
                geom_alluvium(aes(fill = OA), width=0.2, knot.pos = 0, reverse=FALSE) +
                geom_stratum(alpha = 0.5, width=0.2, reverse=FALSE) +
                geom_text(stat = "stratum",
                          aes(label = after_stat(stratum)), reverse=FALSE) +
                scale_fill_futurama(alpha=1) +
                scale_x_discrete(limits = c("OA.Adet.M", "OA.Adet.W","OA.Akaki.W", "OA.Akaki.M"),
                                expand = c(0.05, 0.05)) +
                ggtitle("Overall appreciation by gender and location")+
                guides(fill=guide_legend(title="Combined OA")) +
                theme_minimal()
alluvialplot                
  
#################################???
# compare landraces with improved varieties
##################################
tmp<-all.meta[,52:ncol(all.meta)]
head(tmp)
tmp<-na.omit(tmp)

#split landraces and varieties
lrs<-subset(tmp, status=="Landrace")
ivs<-subset(tmp, status=="Breeding material")
dim(lrs)
lrs<-lrs[,-ncol(lrs)]

dim(ivs)
ivs<-ivs[,-ncol(ivs)]

#select traits to include in the venn
lowbest<-c("DH", "DM", "GFP")

#get improved lines means
lrsqnt<-apply(lrs, 2, function(x) quantile(x))
ivsqnt<-apply(ivs, 2, function(x) quantile(x))
ivsqnt

#plot(lrsmax - ivsmax)

#loop through lrs
listout<-list()
for (i in colnames(lrsqnt)){ #i="GFP.Adet"
  tmplr<-lrs[,i]
  tmpiv<-ivsqnt[,i]
  #check if the phenotype is one such phneotypes for which the lower the better
  if(length(grep(paste0(lowbest, collapse="|"), i))>0){
    listout[[i]]<-rownames(lrs)[which(tmplr<=tmpiv["25%"])]
  } else {
    listout[[i]]<-rownames(lrs)[which(tmplr>=tmpiv["75%"])]
  }
}

names(listout)
length(listout)

#subset the traits prior plotting
tmplist<-listout[c("OA", "PA", "GFP", "GY")]
venn <- Venn(tmplist)
data <- process_data(venn)

vennplot<- ggplot()+
  geom_sf(aes(fill=count), data = venn_region(data)) +
  #geom_sf(size = 1, lty = "dashed", color = "grey", data = venn_setedge(data), show.legend = F) +
  geom_sf_text(aes(label = name), data = venn_setlabel(data)) +
  geom_sf_label(aes(label=count), data = venn_region(data), alpha=0.9) +
  scale_fill_gradient(low="white",high = "lightblue") +
  theme(title = "Landraces outperforming improved varities")+
  theme_void()
vennplot

#get the best variety included in all subsets
which(venn["OA"] %in% venn["PA"])


##################
#make a composite plot with patchwork
##################


layout<-"
        ABBBB
        ABBBB
        ACCCC
        ACCCC
        DDDEE
        DDDEE
      "

combined<- allcorplot + combotop + boxplots + 
          alluvialplot + vennplot + 
          plot_layout(design = layout) +
          plot_annotation(tag_levels = 'a') 
combined


ggsave(combined, file="figures/Fig.2.pdf", width = 12, height = 12)
ggsave(combined, file="figures/Fig.2.jpg", width = 12, height = 12)


########## make sankey graph
# alluvialq$id<-rownames(alluvial)
# tmp<-alluvialq[,c("id", "OA", "GY", "PA")]
# 
# data<- table(melt(tmp, id.var="id")[-2])
# data
# head(data)
# 
# links <- data %>% 
#   as.data.frame() %>% 
#   rownames_to_column(var="source") %>% 
#   gather(key="target", value="value", -1) %>%
#   filter(value != 0)




