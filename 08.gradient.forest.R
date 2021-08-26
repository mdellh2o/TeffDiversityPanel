########################################################################################
# This script performs the LANDSCAPE GENOMICS Gradient forest analysis.
# Date: 11/24/2020
# Author: Leonardo Caproni
########################################################################################

#  Load the libraries
library(vegan)
library(adespatial)
library(BiodiversityR)
library(gradientForest)
library(raster)
library(dismo)
library(maptools)
library(measurements)
library(rgdal)
library(maps)
library(mapdata)
library(rJava) 
library(jsonlite)
library(RStoolbox) 
library(data.table)
library(stringr)

# graph libararies
library(ggplot2)
library(ggspatial)
library(ggsflabel)
library(ggfortify)
library(patchwork)
library(RColorBrewer)

rm(list=ls())
maindir <- "/glab/projects/teff"
setwd(maindir)

options(stringsAsFactors = F)

#load working files
load(file="metadata/pheno.bio.analysis.step.1.Rdata")
load(file="metadata/diversity.analysis.step.1.Rdata")

#Calculate Moran's eigenvector maps (MEMs)
##This approach generates a set of uncorrelated spatial variables too

mem <- dbmem(dist(coord), thresh=1.012, MEM.autocor = "positive", silent = FALSE)
mem.keep <- scores(mem)[,1:10] #keep first 10 MEMs

###################################################################################################################
#Gradient forest analysis
###################################################################################################################

#we need same number of samples with bioclimatic variables and MEMs
coord.tmp <- comb.comb[,c("ID","LON","LAT", "bio1", "bio2", "bio3", 
                                 "bio4", "bio5", "bio6", "bio7", "bio8", 
                                 "bio9", "bio10", "bio11", "bio12",
                                 "bio13", "bio14", "bio15", "bio16",
                                 "bio17", "bio18", "bio19", 
                                 "PC1_bio", "PC2_bio", "PC3_bio", "altitude"
                                 )]
coord.tmp <-na.omit(coord.tmp)

###################################################################################################################
## Download climate data
###################################################################################################################

setwd(inputdir)
#get ethiopian extension
eth <- getData("GADM", country="ETH", level=0)
eth1 <- getData("GADM", country="ETH", level=1)

#get historical environmental data
currentEnv=getData('worldclim', var="bio", res=2.5)
#Climate projections with CMIP5, rcp45 and year 2070
futureEnv=getData('CMIP5', var='bio', res=2.5, rcp=45, model='HE', year=70)
names(futureEnv)=names(currentEnv)
# #Climate projections with CMIP5, rcp26 and year 2070
futureEnv_2670=getData('CMIP5', var='bio', res=2.5, rcp=26, model='HE', year=70)
names(futureEnv_2670)=names(currentEnv)
# #Climate projections with CMIP5, rcp60 and year 2070
futureEnv_6070=getData('CMIP5', var='bio', res=2.5, rcp=60, model='HE', year=70)
names(futureEnv_6070)=names(currentEnv)
# #Projections with CMIP5, rcp85 and year 2070
futureEnv_8570=getData('CMIP5', var='bio', res=2.5, rcp=85, model='HE', year=70)
names(futureEnv_8570)=names(currentEnv)

#crop climate data to study area
model.extent<-extent(eth)
modelEnv=crop(currentEnv,model.extent)

modelFutureEnv=crop(futureEnv, model.extent)
modelFutureEnv_2670=crop(futureEnv_2670, model.extent)
modelFutureEnv_6070=crop(futureEnv_6070, model.extent)
modelFutureEnv_8570=crop(futureEnv_8570, model.extent)

#save current data
save(modelEnv, modelFutureEnv, eth, eth1, coord, file="../metadata/distribution.model.step.1.Rdata")

##############################################################################################
#check redundancy in bioclim variables using the current values
##############################################################################################

modelEnv<-stack(modelEnv)
vif <- ensemble.VIF(
  x = modelEnv,
  a = data.frame(sp),
  VIF.max = 10,
  keep = NULL,
  layer.drops = "bio7",
  factors = NULL,
  dummy.vars = NULL
)

tokeep<-names(vif$VIF.final)

###drop colinear variables from environmental datasets
redbc<-which(names(modelEnv) %in% tokeep)
modelEnv<-modelEnv[[redbc]]
modelEnv<-brick(modelEnv)

# RCP 2.6
modelFutureEnv_2670<-stack(modelFutureEnv_2670)
modelFutureEnv_2670<-modelFutureEnv_2670[[redbc]]
modelFutureEnv_2670<-brick(modelFutureEnv_2670)

# RCP 4.5
modelFutureEnv<-stack(modelFutureEnv)
modelFutureEnv<-modelFutureEnv[[redbc]]
modelFutureEnv<-brick(modelFutureEnv)

#RCP 6.0
modelFutureEnv_6070<-stack(modelFutureEnv_6070)
modelFutureEnv_6070<-modelFutureEnv_6070[[redbc]]
modelFutureEnv_6070<-brick(modelFutureEnv_6070)

#RCP 8.5
modelFutureEnv_8570<-stack(modelFutureEnv_8570)
modelFutureEnv_8570<-modelFutureEnv_8570[[redbc]]
modelFutureEnv_8570<-brick(modelFutureEnv_8570)


###################################################################################################################
# create an object that contains only the climate and MEM spatial variables (no lat/lon)

setwd(maindir)

env.gf.mem <- cbind(coord.tmp[,c("ID","bio2","bio3", "bio4", "bio9",
                                 "bio13", "bio14", "bio15",
                                 "bio18", "bio19")], mem.keep)
								 

row.names(env.gf.mem) <- env.gf.mem[,1]
env.gf.mem <- env.gf.mem[,2:ncol(env.gf.mem)]

#In GF, a maximum number of splits can be defined following the developers suggestion
maxLevel.gf.mem <- log2(0.368*nrow(env.gf.mem)/2)


# Prepare a SNP file
snpGF <- read.table("/glab/projects/teff/gradientForest/full.set", header = T, row.names = 1)

# make sure we cot the same order of samples 
ordGF<-coord.tmp[,1]
snpGF<-subset(snpGF, rownames(snpGF) %in% ordGF)
snpGF<-snpGF[order(match(rownames(snpGF), ordGF)), , drop = FALSE]

save(env.gf.mem, snpGF,maxLevel.gf.mem, file="metadata/metadata.beforeGF.Rdata")

####################################################################################################
#RunGF
####################################################################################################

gf.mem <- gradientForest(cbind(env.gf.mem, snpGF), 
                     predictor.vars=colnames(env.gf.mem), 
                     response.vars=colnames(snpGF), 
                     ntree=500, 
                     maxLevel=maxLevel.gf, trace=T, 
                     corr.threshold=0.50)

# save metadata after GF
save.image("metadata/metadata.afterGF.full.set.Rdata")


# restart from here
load(file="metadata.afterGF.full.set.Rdata")

#arrange variables by importance
by.imp.mem <- names(importance(gf.mem))

#Plot overall importance of our predictors
pdf("gradientForest/full.set/GF_Overall importance.pdf")
plot(gf.mem, plot.type="Overall.Importance")
dev.off()

# plot turnover fucytions
pdf("gradientForest/full.set/GF_TurnoverFunctions.pdf")
plot(gf.mem, plot.type = "C", imp.vars = by.imp.mem, common.scale = T, show.species = F,
     cex.axis = 1, cex.lab = 1.2, line.ylab = 1, 
     par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 2, 2, 2), omi = c(0.2, 0.3, 0.2, 0.4)))
dev.off()


#############################################################################################################
#Now is time to make predictions using the gf function
#############################################################################################################
## use the predict function to "predict" the allelic turnover across the whole
## landscape using our fitted model (gf) and the climate values on the landscape

# Predict allelic turnover actual environment
clim.land <- extract(modelEnv, 1:ncell(modelEnv), df = TRUE) #clim.layer.crop
clim.land <- na.omit(clim.land)

pred.mem <- predict(gf.mem, clim.land[,-1]) 

# #For projections with CMIP5, rcp26 and year 2070
clim.land.2670 <- extract(modelFutureEnv_2670, 1:ncell(modelFutureEnv_2670), df = TRUE) #clim.layer.crop
clim.land.2670 <- na.omit(clim.land.2670)

proj.mem.2670 <- predict(gf.mem, clim.land.2670[,-1]) 

# Predict allelic turnover actual environment
clim.land.4570 <- extract(modelFutureEnv, 1:ncell(modelFutureEnv), df = TRUE) 
clim.land.4570 <- na.omit(clim.land.4570)

proj.mem.4570 <- predict(gf.mem, clim.land.4570[,-1]) 

# #For projections with CMIP5, rcp60 and year 2070
clim.land.6070 <- extract(modelFutureEnv_6070, 1:ncell(modelFutureEnv_6070), df = TRUE) #clim.layer.crop
clim.land.6070 <- na.omit(clim.land.6070)

proj.mem.6070 <- predict(gf.mem, clim.land.6070[,-1]) 

# #For projections with CMIP5, rcp85 and year 2070
clim.land.8570 <- extract(modelFutureEnv_8570, 1:ncell(modelFutureEnv_8570), df = TRUE) #clim.layer.crop
clim.land.8570 <- na.omit(clim.land.8570)

proj.mem.8570 <- predict(gf.mem, clim.land.8570[,-1]) 

#############################################################################################################
#Plot allelic turnover
#############################################################################################################
##These predictions then need to be converted to a color scale for mapping. 
##One way is to use principal components analysis (PCA) on the predictions and use 
##the first three axes to define red, green, blue color scales, respectively. 
##After the values are defined in color space, they can be stacked and mapped.

PCs <- prcomp(pred.mem, center=T, scale=F) #For pred.mem
a1 <- PCs$x[, 1]
a2 <- PCs$x[, 2]
a3 <- PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 +a2 -a1
r <- (r - min(r))/(max(r) - min(r)) * 255
g <- (g - min(g))/(max(g) - min(g)) * 255
b <- (b - min(b))/(max(b) - min(b)) * 255
mask<-modelEnv$bio19 #Precipitation of Coldest Quarter
mask[]<-as.numeric(mask[]>0)
rastR <- rastG <- rastB <- mask
rastR[clim.land$ID] <- r
rastG[clim.land$ID] <- g
rastB[clim.land$ID] <- b
rgb.rast <- stack(rastR, rastG, rastB)

# #For projection rcp85 year 2070
PCs.proj.8570 <- prcomp(proj.mem.8570, center=T, scale.=F) #For proj.mem
r.proj.8570 <- PCs.proj.8570$x[, 1]
g.proj.8570 <- PCs.proj.8570$x[, 2]
b.proj.8570 <- PCs.proj.8570$x[, 3]
r.proj.8570 <- (r.proj.8570 - min(r.proj.8570))/(max(r.proj.8570) - min(r.proj.8570)) * 255
g.proj.8570 <- (g.proj.8570 - min(g.proj.8570))/(max(g.proj.8570) - min(g.proj.8570)) * 255
b.proj.8570 <- (b.proj.8570 - min(b.proj.8570))/(max(b.proj.8570) - min(b.proj.8570)) * 255
mask.proj.8570<-modelFutureEnv_8570$bio19 #Precipitation of Coldest Quarter
mask.proj.8570[]<-as.numeric(mask.proj.8570[]>0)
rastR.proj.8570 <- rastG.proj.8570 <- rastB.proj.8570 <- mask.proj.8570
rastR.proj.8570[clim.land.8570$ID] <- r.proj.8570
rastG.proj.8570[clim.land.8570$ID] <- g.proj.8570
rastB.proj.8570[clim.land.8570$ID] <- b.proj.8570
rgb.rast.proj.8570 <- stack(rastR.proj.8570, rastG.proj.8570, rastB.proj.8570)

#############################################################################################################
# DEFINE CROPPIN AREA
### project vulnerability only on the cropping area, AEZ with at least 2 hits

aez<- readOGR(dsn = "input/AEZ_32/Agro_ecology.shp")

ETH_AEZs <- subset(aez, AEZ31 == "H3"|AEZ31 == "H4"|
                     AEZ31 == "M2"|AEZ31 == "M3"|AEZ31 == "M4"|
                     AEZ31 == "SH2"|AEZ31 == "SH3"|AEZ31 == "SH4"|
                     AEZ31 == "SM2"|AEZ31 == "SM3"|AEZ31 == "SM4")

## now create common garden locations refs on maps (WGS84)

pheno.loc <- data.frame(loc=c("Akaki","Adet"),
                        lat= c(8.8354, 11.2756),
                        lon= c(38.8329, 37.4917))          

pts_pheno<-SpatialPointsDataFrame(data=pheno.loc, coords = pheno.loc[,c("lon", "lat")],
                                  proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))  

pts_pheno<-st_as_sf(pts_pheno, coords = c("lon", "lat"))

#############################################################################################################
# Map allelic turnover in geographic space (only on the cropping area)

# convert to sf to draw Ethiopian boarders
eth_ok <- st_as_sf(eth)
eth1_ok <- st_as_sf(eth1)

# plot allelic turnover
rgb.rast <- crop(rgb.rast, extent(ETH_AEZs))
rgb.rast  <- mask(rgb.rast, ETH_AEZs) 

rgb.rast.proj.8570 <- crop(rgb.rast.proj.8570, extent(ETH_AEZs))
rgb.rast.proj.8570  <- mask(rgb.rast.proj.8570, ETH_AEZs) 


#############################################################################################################
# make biplot of the biological space using ggfortify

PCs.gf.mem <-autoplot(PCs, col =rgb(r, g, b, max = 255),
                      loadings = TRUE, loadings.colour = 'gray30',
                      loadings.label = TRUE, loadings.label.size = 2.5, loadings.label.colour='black') +
  geom_vline(linetype = "dashed", xintercept = 0, color="gray10") +
  geom_hline(linetype = "dashed", yintercept = 0, color="gray10") +
  #xlab("") + ylab("") +
  theme_light() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) #+
  # scale_y_continuous(breaks=NULL) +
  # scale_x_continuous(breaks=NULL)

PCs.gf.mem.proj.8570 <-autoplot(PCs.proj.8570, col =rgb(r.proj.8570, g.proj.8570, b.proj.8570, max = 255),
                          loadings = TRUE, loadings.colour = 'gray20',
                          loadings.label = TRUE, loadings.label.size = 2.5, loadings.label.colour='black') +
  geom_vline(linetype = "dashed", xintercept = 0, color="gray10") +
  geom_hline(linetype = "dashed", yintercept = 0, color="gray10") +
  #xlab("") + ylab("") +
  theme_light() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank()) #+
  # scale_y_continuous(breaks=NULL) +
  # scale_x_continuous(breaks=NULL)

###################################################################???
# plot allelic turnover in the geographic space using ggplot

map <- ggRGB(rgb.rast, alpha = 1) +
  geom_sf(data = eth1_ok, fill=NA, color="gray30", alpha=1) +
  #geom_sf(data=pts, size=2, aes(color=Cluster)) +
  geom_sf(data=pts_pheno, shape=23, fill="yellow", color="black", size=2.5, alpha=0.8 ) +
  geom_sf_label_repel(data = pts_pheno, aes(label = loc),nudge_x = +3.2, nudge_y = +1.4, seed = 10, alpha =0.8) +
  theme_light() +
  xlab("longitude") + ylab("latitude")+
  labs(title = "Current") +
  #scale_color_viridis(discrete=TRUE) +
  annotation_scale(location = "tl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering)


map.proj.8570 <- ggRGB(rgb.rast.proj.8570, alpha = 1) +
  geom_sf(data = eth1_ok, fill=NA, color="gray30", alpha=1) +
  #geom_sf(data=pts, size=2, aes(color=Cluster)) +
  geom_sf(data=pts_pheno, shape=23, fill="yellow", color="black", size=2.5, alpha=0.8 ) +
  geom_sf_label_repel(data = pts_pheno, aes(label = loc),nudge_x = +3.2, nudge_y = +1.4, seed = 10, alpha =0.8) +
  theme_light() +
  xlab("longitude") + ylab("latitude")+
  labs(title = "RCP 8.5", tag="d") +
  #scale_color_viridis(discrete=TRUE) +
  annotation_scale(location = "tl", width_hint = 0.2) +
  annotation_north_arrow(location = "tl", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering)


#############################################################################################################
# Estimate genomic offset, so-called "genomic vulnerability"
#############################################################################################################

temp.2670 <- vector("numeric", length = nrow(proj.mem.2670))
for (i in 1:ncol(proj.mem.2670)) {
  temp.2670 <- temp.2670 + (proj.mem.2670[,i]-pred.mem[,i])^2
}

temp.4570 <- vector("numeric", length = nrow(proj.mem.4570))
for (i in 1:ncol(proj.mem.4570)) {
  temp.4570 <- temp.4570 + (proj.mem.4570[,i]-pred.mem[,i])^2
}

temp.6070 <- vector("numeric", length = nrow(proj.mem.6070))
for (i in 1:ncol(proj.mem.6070)) {
  temp.6070 <- temp.6070 + (proj.mem.6070[,i]-pred.mem[,i])^2
}

temp.8570 <- vector("numeric", length = nrow(proj.mem.8570))
for (i in 1:ncol(proj.mem.8570)) {
  temp.8570 <- temp.8570 + (proj.mem.8570[,i]-pred.mem[,i])^2
}


##############################################################
GenVuln.2670 <- data.frame(sqrt(temp.2670))
GenVuln.4570 <- data.frame(sqrt(temp.4570))
GenVuln.6070 <- data.frame(sqrt(temp.6070))
GenVuln.8570 <- data.frame(sqrt(temp.8570))


GenVuln <- cbind(clim.land[,c(1)], GenVuln.2670, GenVuln.4570, GenVuln.6070, GenVuln.8570 )

colnames(GenVuln)[1] <- "cell_ID"
colnames(GenVuln)[2] <- "RCP 2.6"
colnames(GenVuln)[3] <- "RCP 4.5"
colnames(GenVuln)[4] <- "RCP 6.0"
colnames(GenVuln)[5] <- "RCP 8.5"
summary(GenVuln)

#############################################################################################################
clim.land2 <- extract(modelEnv, 1:ncell(modelEnv), df = TRUE)
clim.land2 <- cbind(coordinates(modelEnv), clim.land2)
clim.land2 <- na.omit(clim.land2)
# assign coordinates to each of the pixels
genVuln <- cbind(clim.land2[,c(1,2)], GenVuln)
#names(genVuln)


# plot genomic offset only on the cropping area

coordinates(genVuln) <- ~ x + y
gridded(genVuln) <- TRUE

genVuln.rast.8570 <- raster(genVuln, "RCP 8.5")
genVuln.rast.8570 <- crop(genVuln.rast.8570, extent(ETH_AEZs))
genVuln.rast.8570 <- mask(genVuln.rast.8570, ETH_AEZs)
names(genVuln.rast.8570)<-"Offset"

# RCP 8.5"
#########################################################################################

vuln.8570 <- ggplot() +
  geom_raster(data=genVuln.rast.8570, aes(x=x, y=y, fill=Offset)) + 
  scale_fill_gradientn(colours=rev(brewer.pal(7,"Spectral")), limits=c(0, 0.0058), na.value=NA) +
  xlab("longitude") + ylab("latitude")+
  coord_quickmap() + 
  #coord_sf() +
  theme_light() +
  labs(title="RCP 8.5") +
  geom_sf(data = eth1_ok, fill=NA, color="gray45", alpha=0.8) +
  geom_sf(data=pts_pheno, shape=23, fill="yellow", color="black", size=2.5, alpha=0.8 ) +
  geom_sf_label_repel(data = pts_pheno, aes(label = loc),nudge_x = +3.2, nudge_y = +1.4, seed = 10, alpha =0.8) +  xlab("longitude") + ylab("latitude")+
  #geom_sf(data=pts, size=2, aes(color=Cluster)) + 
  #scale_color_viridis(discrete=TRUE) +
  #labs(title = "Genomic vulnerability with future climate: rcp45, year 2070") +
  annotation_scale(location = "tr", width_hint = 0.2) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.2, "in"), pad_y = unit(0.4, "in"),
                         style = north_arrow_fancy_orienteering)


# save your workspace
save.image(file="gradientForest/full.set/genVuln.crop.area.Rdata")
