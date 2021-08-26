########################################################################################
# This script merges data from passport file and Table S2 from Woldeyohannes et al. (2020)
# Date: 11/24/2020
# Author: Leonardo Caproni
########################################################################################


#1) Load packages and set directories
library(raster)
library(maptools)
library(rgdal)

inputdir<-"/glab/projects/teff/geographic.data"
setwd(inputdir)

#2) Read the passport data
csvfile<-"../phenotypes/passport.info.qualitative.traits.csv"
passport <- read.csv(csvfile,stringsAsFactors = FALSE)


#Read Table S2 from Woldeyohannes et al 2020
wold<-read.delim("../phenotypes/Woldeyohannes.et.al.Table.S2.pheno.txt")

#make sure everything is lowercase
wold$Accession<-tolower(wold$Accession)
passport$ACCESSION<-tolower(passport$ACCESSION)

length(which(wold$Accession %in% passport$ACCESSION))

#remove useless information
head(passport)
subpass<-passport[,1:13]

passnew<-merge(subpass, wold, by.x="ACCESSION", by.y="Accession", all=T)
dim(passnew)
head(passnew)
names(passnew)

#remove bioclimatic information
passnew<-passnew[,-grep("bio|Altitude|ph", colnames(passnew))]
#remove redundant columns
passnew<-passnew[,-which(colnames(passnew) %in% c("Region","Zone" ,"District" ,"Locality"))]
names(passnew)

#overwrite passport
passport<-passnew
passport [1:20, 12:18]

####convert DMS to decimal degrees
# change the degree symbols to somthing that could be understood by R/measurements
passport$DMSlat = sub('-', ' ', passport$Lat)
passport$DMSlat = sub('-', '.', passport$DMSlat)
passport$DMSlat = sub('-N', '', passport$DMSlat)

passport$DMSlon = sub('-', ' ', passport$Long)
passport$DMSlon = sub('-', '.', passport$DMSlon)
passport$DMSlon = sub('-E', '', passport$DMSlon)

passport$DMSlat2 = sub('-', ' ', passport$LAT)
passport$DMSlat2 = sub('-', '.', passport$DMSlat2)
passport$DMSlat2 = sub('-N', '', passport$DMSlat2)

passport$DMSlon2 = sub('-', ' ', passport$LON)
passport$DMSlon2 = sub('-', '.', passport$DMSlon2)
passport$DMSlon2 = sub('-E', '', passport$DMSlon2)

passport [1:20, 33:36]

#merge the two column from wold and passport in decimal degrees
passport$DMSlat <- ifelse(is.na(passport$DMSlat), passport$DMSlat2, passport$DMSlat)
passport$DMSlon <- ifelse(is.na(passport$DMSlon), passport$DMSlon2, passport$DMSlon)

# convert from DMS to decimal degrees
passport$lat = measurements::conv_unit(passport$DMSlat, from = 'deg_dec_min', to = 'dec_deg')
passport$lon = measurements::conv_unit(passport$DMSlon, from = 'deg_dec_min', to = 'dec_deg')

passport$lat<-as.numeric(passport$lat)
passport$lon<-as.numeric(passport$lon)

#extract mapped features and put them in a SPDF
pasmap<-passport[which(!is.na(passport[,"lat"])),]
pass<-SpatialPointsDataFrame(data=pasmap, coords =pasmap[,c("lon", "lat")],
                             proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))  

#plot(pass, main="Map of Teff")

#extract extension and attributes
minmax<-extent(pass)
att<-pass@data 

#get elevation
alt<-getData('alt', country='ETH')

#get additional data
et1<-getData('GADM' , country="ETH", level=1)

#get bioclim at max resolution using tef extension to get the right tile
bioc <- getData("worldclim",var="bio",lon = minmax[1], lat=minmax[3],res=0.5)

#extract values
altitude<-extract(alt,pass)
values <- extract(bioc,pass)
biovar <- cbind.data.frame(ACCESSION=att[,1],altitude, values)

#bring together with passport
names(biovar)
names(passport)
passbio<-merge(passport, biovar, by = "ACCESSION", all=T)

#check number of accessions mapped
head(passbio)
length(which(!is.na(passbio$lat)))
length(which(!is.na(passbio$lon)))

#rearrange the dataset
names(passbio)
passbio<-passbio[,c(1:11,37,38,33,34,39,16:32,40:58)]

#fix the occurrence of comma in LOCALITY column
hit<-apply(passbio, 2, function(x) grep(",", x))
passbio[hit$LOCALITY,"LOCALITY"]<-sub(",", " ", passbio[hit$LOCALITY,"LOCALITY"])

length(which(!is.na(passbio$ID)))

#save object for 
save(passbio, file="../metadata/passport.data.qualitative.traits.bioclim_OK.Rdata")
write.csv(passbio, file="../metadata/passport.data.qualitative.traits.bioclim_OK.csv", row.names=F, quote=F)

