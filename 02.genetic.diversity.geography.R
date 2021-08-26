########################################################################################
# This script does perform DIVERSITY ANALYSIS PCA & DAPC and maps the diversity 
# onto the geographic space
# Date: 11/24/2020
# Authors:Leonardo Caproni, Sessen Daniel Iohannes and Matteo Dell'Acqua 
########################################################################################

rm(list=ls())
options(stringsAsFactors = F)

#Install packages, set directories and load dataset
library(gplots)
library(RColorBrewer)
library(adegenet)
library(car)
library(tidyverse)
library(stringr)
library(pegas)
library(ape)
library(seqinr)
library(ade4)
library(factoextra)
library(CMplot)
library(ggsci)
library(LDheatmap)
library(RColorBrewer)
library(genetics)
library(EMMREML)
library(compiler)
library(scatterplot3d)
library(hierfstat)
library(poppr)
library(hierfstat)
library(viridis)
library(raster)
library(sf)
library(rgdal)
library(dplyr)
library(fasterize)
library(mapview)
library(ggspatial)
library(gridExtra)
library(ggrepel)
library(patchwork)
library(ggplot2)
library(rgl)
library(magick)
library(vcfR)
library(data.table)

maindir <- "/glab/projects/teff"
setwd(maindir)

############################################################################################
#1) Kinship visualization
############################################################################################

#load Identity By Status Matrix and its IDs computed using plink 1.90
mibs<-read.table(file = "snps_HaplotypeCaller/no_wild/plink/pruned_100_10_03/plink.mibs", sep = " ", header = F)
mibs$V367<-NULL

mIDs<-read.table(file = "snps_HaplotypeCaller/no_wild/plink/pruned_100_10_03/plink.mibs.id", sep = "\t", header = F)
mIDs$V1<- NULL
colnames(mIDs) <- c("IDs")

mIDs <- data.frame(gsub("T", "T_", mIDs$IDs))
colnames(mIDs) <- c("IDs")

##insert IDs
mat <- as.matrix(mibs)
rownames(mat)<-mIDs$"IDs"
colnames(mat)<-mIDs$"IDs"

#PLOT a heatmap
pl <- colorRampPalette(c("#ff3419",
                         "#01529f"))

pdf(paste0("figures/IBS.pairwise.pruned_100_10_03.pdf"))
heatmap.2(mat, col = pl(50), scale = "none", trace="none", labRow = FALSE, labCol = FALSE)
dev.off()

# get order of the dendrogram created using "hclust()"
x<-heatmap.2(mat, col = pl(50), scale = "none", trace="none") #class(x)
order <- colnames(x$carpet)
rm(x)

sort <- str_replace(order, "T_", "T")
sort_bed <- data.frame(fam.id=sort,
                       id=sort)
# sort_bed[1:10,1:2]
write.table(sort_bed, file = "metadata/sort.bed.txt", sep='\t', quote = FALSE, row.names = FALSE, col.names = FALSE)

############################################################################################
#2) DIVERSITY ANALYSIS PCA & DAPC
############################################################################################

## 2.1) Data preprocessing
myvcf <- read.vcfR ("snps_HaplotypeCaller/no_wild/plink/pruned_100_10_03/pruned_100_10_03.vcf", verbose = FALSE)
genind <- vcfR2genind(myvcf)
genind #366 individuals; 2100 
summary(genind@loc.n.all) # all loci are biallelic

genind_scale <- scaleGen(genind, NA.method="mean") # Impute missing data using mean. This genotype dataset is used for PCA

## 2.2 Calculating basic population genetics statistics
div <- summary(genind) 
class(div$Hexp)
div$Hobs <- as.array(div$Hobs)
div_het <- cbind(div$Hexp, div$Hobs)
div_het <- as.data.frame(div_het) #Create dataframe displaying observed and expected heterozygosity

pdf("figures/Heterozygosity.pdf")
plot(div$Hobs,div$Hexp, xlab="Hobs", ylab="Hexp", 
     main="Expected heterozygosity as a function of observed heterozygosity per locus", cex=0.1, pch=1, col="blue")
#abline(0,1,col="red") #Create expected_het vs. observed_het plot
dev.off()

bartlett.test(list(div$Hexp, div$Hobs)) #assess homogeneity of variance
t.test(div$Hexp, div$Hobs, paired=TRUE, var.equal=FALSE) #assess significance for differences between observed and expected heterozigosity


############################################################################################
#3. PCA
############################################################################################

## 3.1 Diversity analysis using PCA
geno.pca <-dudi.pca(genind_scale,cent=TRUE,scale=TRUE,scannf=FALSE, nf=3)
geno.pca_scores <- data.frame(geno.pca$li)# order scores
geno.pca_scores <- data.frame(scale(geno.pca_scores)) #scale scores
#geno.pca_scores [1:2,1:3]


pdf("figures/scree_plot.pdf")
fviz_eig(geno.pca) #visualize scree plot
dev.off()

get_eigenvalue(geno.pca) #Visualize loadings statistics (variance explained etc...)

#plot 2 PCAs
pc12<- ggplot(geno.pca_scores, aes(x= Axis1, y = Axis2)) +
  geom_point(size=1) +
  coord_cartesian() +
  theme_minimal()
ggsave(pc12, file="figures/pc12.geno.jpg", height=8, width=8)


pc13 <- ggplot(geno.pca_scores, aes(Axis1, Axis3)) +
  geom_point(size=1) + 
  coord_cartesian() +
  theme_minimal()
ggsave(pc13, file="figures/pc13.geno.jpg", height=8, width=8)

## DISCRIMINANT ANALYSIS OF PRINCIPAL COMPONENTS
grp <- find.clusters(genind, 
                     max.n.clust=30,
                     n.pca = 365,
                     n.clust = 6
                     ) 
dapc1 <- dapc(genind, grp$grp, n.pca = 365, 
              n.da = 3,
              ) #kept 3 discriminant functions
#dapc1
#scatter(dapc1, cell = 0, pch = 18:23, cstar = 2, mstree = FALSE, lwd = 2, lty = 2, posi.da="topright")

##Plot PCA by DAPC assig
assig <- as.data.frame(grp$grp)
geno <- rownames(assig)
rownames(assig) <- NULL
assig <- cbind(geno,assig)
colnames(assig)<-c("geno", "Cluster")

geno1 <- rownames(geno.pca_scores)
rownames(geno.pca_scores) <- NULL
geno.pca_scores<- cbind(geno1,geno.pca_scores)

#geno.pca_scores[1:4,1:3]

geno.pca_scores_assig <- merge(x= assig, y= geno.pca_scores, by.x = "geno", by.y= "geno1" )
#geno.pca_scores_assig[1:4,1:5]
geno.pca_scores_assig$geno <- str_replace(geno.pca_scores_assig$geno, "T", "T_") #ID with underscore


outpca12 <- ggplot(geno.pca_scores_assig, aes(x= Axis1, y = Axis2, color=Cluster)) +
  #scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  scale_color_nejm() +
  geom_point(size=2) +
  coord_cartesian()+
  labs(fill = "DAPC cluster", 
       x= "PC1 (5.71%)",
       y=  "PC2 (5.27%)", 
       size = "") +
  theme_light()
#labs(size = "")
outpca12 <- outpca12 + guides(color = FALSE)

outpca13 <- ggplot(geno.pca_scores_assig, aes(x= Axis1, y = Axis3, color=Cluster)) +
  #scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  scale_color_nejm() +
  geom_point(size=2) +
  coord_cartesian()+
  labs(fill = "DAPC cluster", size = "") +
  theme_light()+
  labs(
    x= "PC1 (5.71%)",
    y=  "PC3 (3.65%)",
    color= "DAPC cluster",
    size = "")

ggsave(outpca12, file="figures/DAPC12.jpg", height=8, width=8)
ggsave(outpca13, file="figures/DAPC13.jpg", height=8, width=8)

save.image(file="metadata/diversity.analysis.step.1.Rdata")

##Plot PCA by varieties
load(file= "metadata/passport.data.qualitative.traits.bioclim_OK.Rdata")
load(file="metadata/diversity.analysis.step.1.Rdata")

#For PCA by varieties
varieties <- passbio[,c(2,7)]
varieties$TYPE <- as.factor(varieties$TYPE)
levels(varieties$TYPE)[levels(varieties$TYPE)=="breeding material"] <- "breeding material"
levels(varieties$TYPE)[levels(varieties$TYPE)=="modern variety"] <- "breeding material"
colnames(varieties)[1] <- "geno1"
colnames(varieties)[2] <- "Status"
varieties$geno1 <- str_replace(varieties$geno1, "T_", "T")

#For DAPC by varieties
varieties_dapc <- varieties
varieties_dapc$geno1 <- str_replace(varieties_dapc$geno1, "T", "T_")


biostatus <- merge(geno.pca_scores, varieties, by="geno1", all.x=TRUE)
biostatus <- within(biostatus, Status[geno1 == 'T4b'] <- 'landrace')
biostatus <- within(biostatus, Status[geno1 == 'T459'] <- 'breeding material')


dapcstatus <- merge(geno.pca_scores_assig, varieties_dapc, by.x="geno", by.y="geno1", all.x=TRUE)
dapcstatus <- within(dapcstatus, Status[geno== 'T_4b'] <- 'landrace')
dapcstatus <- within(dapcstatus, Status[geno == 'T_459'] <- 'breeding material')

pcavar12 <- ggplot(biostatus, aes(x= Axis1, y = Axis2, color=Status)) +
  #scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  scale_color_nejm() +
  geom_point(size=2) +
  coord_cartesian()+
  labs(x= "PC1 (5.71%)",
       y=  "PC2 (5.27%)", 
       size = "") +
  theme_light()
#labs(size = "")
pcavar12 <- pcavar12 + guides(color=FALSE)

pcavar13 <- ggplot(biostatus, aes(x= Axis1, y = Axis3, color=Status)) +
  #scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  scale_color_nejm() +
  geom_point(size=2) +
  coord_cartesian()+
  labs( x= "PC1 (5.71%)",
       y=  "PC3 (3.65%)", 
       size = "") +
  theme_light()
#labs(size = "")
pcavar13 <- pcavar13 + theme(legend.title = element_text(face = "bold"))


pcavar <- (pcavar12+pcavar13) +
  plot_annotation(
    title = 'PCA: landraces vs. breeding material') +
  plot_annotation(tag_levels = 'A')

dapcvar12 <- ggplot(dapcstatus, aes(x= Axis1, y = Axis2, color=Cluster, shape=Status)) +
  #scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  scale_color_nejm() +
  geom_point(size=2) +
  scale_shape_manual(values = c(17, 16))+
  coord_cartesian()+
  labs(fill = "DAPC cluster", 
       x= "PC1 (5.71%)",
       y=  "PC2 (5.27%)", 
       size = "") +
  theme_light()
#labs(size = "")
dapcvar12 <- dapcvar12 + theme(legend.position = "none")

dapcvar13 <- ggplot(dapcstatus, aes(x= Axis1, y = Axis3, color=Cluster,shape=Status)) +
  #scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  scale_color_nejm() +
  geom_point(size=2) +
  scale_shape_manual(values = c(17, 16))+
  coord_cartesian()+
  labs(fill = "DAPC cluster", size = "") +
  theme_light()+
  labs(
    x= "PC1 (5.71%)",
    y=  "PC3 (3.65%)",
    color= "DAPC cluster",
    size = "")
dapcvar13 <- dapcvar13 + theme(legend.title = element_text(face = "bold"))


dapcvar <- (dapcvar12+dapcvar13) +
  plot_annotation(
    title = 'PCA: DAPC clusters and biological status') +
  plot_annotation(tag_levels = 'A')

ggsave(pcavar, file="figures/PCAVAR.jpg", height=6, width=10)
ggsave(pcavar12, file="figures/PCAVAR12.jpg", height=8, width=8)
ggsave(pcavar13, file="figures/PCAVAR13.jpg", height=8, width=8)

ggsave(dapcvar, file="figures/DAPCVAR.jpg", height=6, width=10)
ggsave(dapcvar12, file="figures/DAPCVAR12.jpg", height=8, width=8)
ggsave(dapcvar13, file="figures/DAPCVAR13.jpg", height=8, width=8)

#save outputs
save.image(file="metadata/diversity.varieties.Rdata")


############################################################################################
#4) Map teff genotypes on Ethiopian map
############################################################################################
####restart from here

#PRELIMINARY
setwd(maindir)
load(file="metadata/diversity.analysis.step.1.Rdata")
#passport.data.qualitative.traits.bioclim.Rdata_OK

aez<- readOGR(dsn = "geographic.data/AEZ_32/Agro_ecology.shp")
#coord<-read.csv("input/acc_coordinates.csv")
#coord<-read.delim("input/acc_coord_UPDATED.txt")
names(passbio)


coord <- passbio[, which(colnames(passbio) %in% c("ID","REGION", "lat","lon", "altitude"))]

#head(coord)
#coord [1:20, 1:3]

#get coord of genotyped accessions & omit NAs
coord <- coord[coord$ID %in% mIDs$IDs,]
coord <- na.omit(coord)

#now include DAPC assicgnation
assig$geno <- str_replace(assig$geno, "T", "T_")
coord <- merge(coord, assig, by.x="ID", by.y="geno")

#Convert foreign object into an sf object
pts<-st_as_sf(coord, coords = c("lon", "lat"))
pts$LON<-coord$lon
pts$LAT<-coord$lat

aez2 <- st_as_sf(aez)
st_crs(pts)<-st_crs(aez2)

pdf("figures/bioclim_distributions.pdf")
plot(pts)
dev.off()

#get intersection bw points and AEZs
inter<-st_intersection(aez2, pts)
#head(inter)

#get relevant zones
dftmp<-data.frame(inter)
hitzones<-unique(dftmp[,"AEZ31"])


#add a factor for plotting
aez2$AEZ<-as.character(aez2$AEZ31)
aez2$AEZ[which(!aez2$AEZ %in% hitzones)]<-"ZNR"
aez2$AEZ<-as.factor(aez2$AEZ)

#include a rougher AEZ definition
aez2$AEZ_type<-as.factor(sub("[0-9]$", "", aez2$AEZ))
aez2$AEZ_type

head(aez2)

#create nice color palette(s)
#nicecols <- c ("#de9342","#b05cc6","#6db643", "#d1417b","#55b684","#cb5742","#4bafd0","#b4ad49","#7178ca", "#577936","#be6b92","#9d7239")

nicecols3<- c ("#657452",
               "#74e451",
               "#5f5848",
               "#d3dd3e",
               "#45813b",
               "#d0dfb2",
               "#6ead2f",
               "#a0aa7e",
               "#a8e88a",
               "#686e2a",
               "#d9d778",
               "#62b866",
               "#9f9b3b"
               )

#make map
outmap<-ggplot(data = aez2) +
  geom_sf(aes(fill=AEZ), alpha=0.4) + 
  geom_sf(data=pts, size=2.5, aes(color=Cluster)) +
  guides(color=FALSE) +
  labs(fill = "Agroecological Zones", size = "") +
  #scale_fill_discrete(name = "Agroecological Zones", labels = c("x", "B", "C")) +
  theme_light() +
  scale_fill_manual(values=c(nicecols3,"white")) +
  #scale_color_manual(values=myCol)+
  #scale_color_viridis(discrete=TRUE) +
  scale_color_nejm() +
  annotation_scale(location = "tr", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering)

#save a .jpg
ggsave(outmap, file="figures/AEZ.map_6.clusters.jpg", height=7, width =9)

#make another type of map
#make map
outmap1<-ggplot(data = aez2) +
  geom_sf(aes(fill=AEZ_type), alpha=0.6) + 
  geom_sf(data=pts, size=2.5, aes(color=Cluster)) +
  guides(color=FALSE) +
  labs(fill = "Agroecological Zone Type", size = "") +
  #scale_fill_discrete(name = "Agroecological Zones", labels = c("x", "B", "C")) +
  theme_light() +
  scale_fill_manual(values=c(gray.colors(length(unique(aez2$AEZ_type))-1),"white")) +
  #scale_color_viridis(discrete=TRUE) +
  #scale_color_manual(values=myCol)+
  scale_color_nejm() +
  annotation_scale(location = "tr", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering)

#save a .jpg
ggsave(outmap1, file="figures/AEZ.type.map_6.clusters.jpg", height=7, width =9)


dftmp$AEZ_type<-sub("[0-9]$", "", dftmp$AEZ31)
#names(dftmp)

dapcAGZ_type <-ggplot(data=dftmp, aes(x=AEZ_type, fill=Cluster)) + 
  geom_bar(stat="count") +
  #scale_fill_manual(values=myCol)+
  #scale_fill_viridis(discrete=TRUE) +
  scale_fill_nejm() +
  labs(x="Agroecological Zone Type", y="count") +
  guides(fill=guide_legend(title="DAPC cluster")) +
  theme_light()


dapcAGZ <-ggplot(data=dftmp, aes(x=AEZ31, fill=Cluster)) + 
  geom_bar(stat="count") +
  #scale_fill_manual(values=myCol)+
  #scale_fill_viridis(discrete=TRUE) +
  scale_fill_nejm() +
  labs(x="Agroecological Zone", y="count") +
  guides(fill=guide_legend(title="DAPC cluster")) +
  theme_light()


#save a .jpg
ggsave(dapcAGZ_type, file="figures/AEZ_type.barplot_6_clusters.jpg", height=3, width =8)
ggsave(dapcAGZ, file="figures/AEZ.barplot_6_clusters.jpg",height=3, width =8)

############################################################################################
# 5. PLOT ITEM 1 PROTOTYPE Diversity analysis
############################################################################################

item <- (outpca12 | outpca13) /
        outmap /
        dapcAGZ +
  plot_layout(heights=c(2,4,1)) +
  plot_annotation(
    title = 'Teff landrace diversity across Ethiopian Agroecological Zones'#,
    # subtitle = 'These 4 plots reveal yet-untold secrets about teff',
    # caption = 'Disclaimer: None of these plots are insightful'
    ) +
  plot_annotation(tag_levels = 'A')

ggsave(item, file="figures/item1_6_CLUSTERS_AGZ.jpg", height=14, width =10)



item1 <- (outpca12 | outpca13) /
  outmap1 /
  dapcAGZ_type +
  plot_layout(heights=c(2,4,1)) +
  plot_annotation(
    title = 'Teff landrace diversity across Ethiopian Agroecological Zone Types'#,
    #subtitle = 'These 4 plots reveal yet-untold secrets about teff',
    #caption = 'Disclaimer: None of these plots are insightful'
    ) +
  plot_annotation(tag_levels = 'A')


ggsave(item1, file="figures/item1_6_CLUSTERS_AGZ_type.jpg", height=14, width =11)

pca.comb <- outpca12 + outpca13
ggsave(pca.comb, file="figures/PCAs_6clusters.jpg", height=6, width =13)


############################################################################################
# 6. export metadata information
############################################################################################
# head(dftmp)
# names(dftmp)

metaout<-dftmp[,c("ID", "AEZ_type", "AEZ31", "REGION", "LON", "LAT")]
out<-merge(metaout, geno.pca_scores_assig, by.x="ID", by.y="geno", all=T)
#out<- out[,c(1:6,8,7,9:ncol(out))]
dim(out)
names(out)
save(out, file="metadata/metadata.diversity.analysis.Rdata")