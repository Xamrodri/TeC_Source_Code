################################################################################
# T&C OUTPUT:
#  READ WEEKLY MAPS FROM RASTER-FILES & PLOT EARLY WILTING
#
# 
# TODO: -
# 
# NEW:  
#       
#
#
# 2023/11/07
#
# Pascal Buri | High Mountain Glaciers and Hydrology | 
#  Swiss Federal Institute for Forest, Snow and Landscape Research, WSL |
#  Zürcherstrasse 111, 8903 Birmensdorf | pascal.buri@wsl.ch
################################################################################
# clear entire workspace (excl. packages)
rm(list = ls())
gc()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')


#===============================================================================
# SETTINGS
#===============================================================================
## define output name (=main folder name)
# outnm<-'v20231006_Aaretal_2557days'
# outnm<-'v20231006_Aargauer_Rhein_2557days'
# outnm<-'v20231006_Berner_Oberland_2557days'
# outnm<-'v20231006_Bodensee_2557days'
# outnm<-'v20231006_Engadin_2557days'
# outnm<-'v20231006_Glatt_2557days'
# outnm<-'v20231006_Hinterrhein_2557days'
# outnm<-'v20231006_Jura_2557days'
# outnm<-'v20231006_Leman_2557days'
# outnm<-'v20231006_Limmat_2557days'
# outnm<-'v20231006_Neuchatel_2557days'
# outnm<-'v20231006_Oberwallis_2557days'
# outnm<-'v20231006_Reuss_2557days'
# outnm<-'v20231006_Rheintal_2557days'
# outnm<-'v20231006_Saane_2557days'
outnm<-'v20231006_Schaffhausen_2557days'
# outnm<-'v20231006_Ticino_2557days'
# outnm<-'v20231006_Thur_2557days'
# outnm<-'v20231006_Toess_2557days'
# outnm<-'v20231006_Unterwallis_369days'
# outnm<-'v20231006_Vorderrhein_2557days'
# outnm<-'v20231006_Zuercher_Rhein_2557days'
###
# outnm<-'v20231006_Test-Aletsch_2557days'
# outnm<-'v20231006_Test-Grossbach_2557days'
# outnm<-'v20231006_Test-Poschiavino_2557days'
# outnm<-'v20231006_Test-Rein_da_Sumvitg_2557days'
# outnm<-'v20231006_Test-Rietholzbach_2557days'
# outnm<-'v20231006_Test-Roggiasca_2557days'
# outnm<-'v20231006_Test-Rosegbach_2557days'


#===============================================================================
# PATHS
#===============================================================================
# detect current machine
if(Sys.getenv('computername') == 'PASCAL-THINK'){
  root<-'D:'
  path_onst<-'C:/Users/Pascal/switchdrive/WSL'  #Path to online storage
}
if(Sys.getenv('computername') == 'WSL26871'){
  root<-'E:'
  path_onst<-'C:/Users/Buri/switchdrive/WSL'  #Path to online storage
}

## Path to stored data
path_data<-root&'/T&C/OUTPUTSTORAGE/DISTRIBUTED/'&outnm

## Path to input data
path_input<-root&'/T&C/INPUTS'

## Path to store data
path_store<-root&'/T&C/POSTPROCESSED/DISTRIBUTED'

## Path to code scripts etc.
path_code<-path_onst&'/T&C/Code/POSTprocessing/CH2022'

## Path to functions file
path_func<-path_onst&'/T&C/Code/Functions/Functions.r'

## Path to GIS data
path_gis<-root&'/QGIS/TC_CH'

##define filename of hillshade
fn_hs<-path_gis&'/DEM/HS_AW3D250_EPSG-2056.tif'

##define filename of Brun2020-earlywilting map
fn_ewm<-path_gis&'/EarlyWilting/Early_wilting_32TMT.tif'


#===============================================================================
# LOAD PACKAGES
#===============================================================================
# define packages to load
library(raster)
library(sf)
library(rgdal)         
library(rgeos)  
library(reshape2)
library(ggplot2)
library(stringr)      #str_replace()
library(wesanderson)  #wes_palette()
library(erer)         #write.list()
library(lutz)         #tz_lookup()
library(viridis)      #scale_color_viridis()
library(ggsci)        #scale_color_jco()
library(R.matlab)     #readMat()
library(datetimeutils)#matlab2POS()
library(scanstatistics) #flipud()
library(matlab)       #reshape()
library(maptools)     #unionSpatialPolygons()
library(plotly)
library(doParallel)   #foreach()
library(doSNOW)
library(rasterVis)      #levelplot()

# define no. of digits printed in console
options('scipen'=100, 'digits'=4)

## define corresponding timezone
# TZ<-'Europe/Zurich'  #considers summer time (CET/CEST)
TZ<-'Africa/Algiers'  #does not have daylight saving

# set time zone & change local time from german to english for plotting
# (check in sessionInfo() ) 
Sys.setenv(TZ=TZ)
Sys.setlocale('LC_TIME', 'C')

# access functions file
source(path_func)


#===============================================================================
# READ SIMULATION FILE << AAA_TC_XXX.m >> OR << BBB_TC_XXX_reinit.m >>
#===============================================================================
setwd(path_data)
files<-list.files()
##read AAA-file
simfile<-files[startsWith(files,'AAA_')]
simfile<-readLines(simfile,warn=FALSE)
reinit=FALSE


#===============================================================================
# GET SIMULATION NAMES
#===============================================================================
VERSION<-grep(pattern="^VERSION",simfile)
VERSION<-simfile[VERSION]
VERSION<-str_split(VERSION,"'")[[1]][2]

DOMAIN<-grep(pattern="^DOMAIN",simfile)
DOMAIN<-simfile[DOMAIN]
DOMAIN<-str_split(DOMAIN,"'")[[1]][2]

simnm<-VERSION&'_'&DOMAIN
simnm1<-simnm

# ##check for BBB-file (re-initialization)
# simnm2<-files[startsWith(files,'BBB_')]
# if(length(simnm2) == 1){
#   reinit=TRUE
#   simfile2<-readLines(simnm2,warn=FALSE)
#   ##get secondary simulation name
#   ## VERSION
#   VERSION<-grep(pattern="VERSION",x=simfile2)
#   VERSION<-simfile2[VERSION[1]]
#   VERSION<-str_match(VERSION,"'\\s*(.*?)\\s*'")[,2]
#   ## DOMAIN
#   DOMAIN<-grep(pattern="DOMAIN",x=simfile2)
#   DOMAIN<-simfile2[DOMAIN[1]]
#   DOMAIN<-str_match(DOMAIN,"'\\s*(.*?)\\s*'")[,2]
#   simnm2<-VERSION&'_'&DOMAIN
#   simnm<-simnm2
# }
# rm(files)

if(reinit == F){simf<-simfile}
if(reinit == T){simf<-simfile2}


#===============================================================================
# CREATE NEW FOLDERS
#===============================================================================
new_dir<-path_store&'/'&outnm
dir.create(new_dir&'/PLOTS')
path_fig<-new_dir&'/PLOTS/WILTING'
dir.create(path_fig)


#===============================================================================
# READ MODELLED WEEKLY MAPS
#===============================================================================
setwd(new_dir&'/RASTERS/WEEKLY')

###### LAI_H ########
VOI<-'LAI_H'
files<-list.files()
br<-brick(files[which(startsWith(files,VOI))[1]])
projec<-projection(br)
##names(br)
## 07.Apr.2018-27.Oct.2018
idx<-132:161
br<-br[[idx]]
ext<-extent(br)
bb<-as(ext,'SpatialPolygons')
projection(bb)<-projec
##check
#plot(br[[1]])
#plot(bb,add=T)
dates<-str_split(names(br),'X')
dates<-do.call('cbind',dates)[2,]
dates<-str_replace(dates,'[.]',' ')
dates<-str_replace(dates,'[.]',' ')
DATES<-as.POSIXct(dates,format='%d %b %Y')
lai_h<-br
rm(VOI,br)

###### LAI_L ########
VOI<-'LAI_L'
br<-brick(files[which(startsWith(files,VOI))[1]])
br<-br[[idx]]
lai_l<-br
rm(VOI,br)

lai<-lai_h+lai_l
rm(lai_h,lai_l)


###### LAIdead_H ########
VOI<-'LAIdead_H'
br<-brick(files[which(startsWith(files,VOI))[1]])
br<-br[[idx]]
laidead_h<-br
rm(VOI,br)

###### LAIdead_L ########
VOI<-'LAIdead_L'
br<-brick(files[which(startsWith(files,VOI))[1]])
br<-br[[idx]]
laidead_l<-br
rm(VOI,br)

laidead<-laidead_h+laidead_l
rm(laidead_h,laidead_l)


#===============================================================================
# READ MODELLED DAILY MAPS
#===============================================================================
setwd(new_dir&'/RASTERS/DAILY')

###### O_H ########
VOI<-'OH'
files<-list.files()
br<-brick(files[which(startsWith(files,VOI))[1]])
##names(br)
## 01.Apr.2018-01.Oct.2018
idx<-549:763
br<-br[[idx]]
o_h<-br
rm(VOI,br)

###### O_L ########
VOI<-'OL'
files<-list.files()
br<-brick(files[which(startsWith(files,VOI))[1]])
br<-br[[idx]]
o_l<-br
rm(VOI,br)

o<-o_h+o_l
names(o)<-names(o_h)
rm(o_h,o_l)


#===============================================================================
# READ HILLSHADE
#===============================================================================
hs<-raster(fn_hs)
hs<-crop(hs,bb)


#===============================================================================
# READ INITIAL CONDITIONS MAPS
#===============================================================================
# %%% #1) Evergreen Tree Forest
# %%% #2) Mixed Evergreen (50%) Decidous (50%)
# %%% #3) Decidous
# %%% #4) Meadow - Grasses
# %%% #5) Rock
# %%% #6) Water
# %%% #7) Ice (bedrock)
pft<-raster(path_store&'/'&outnm&'/RASTERS/INIT-COND/VEG_CODE.tif')
osat<-raster(path_store&'/'&outnm&'/RASTERS/INIT-COND/osat.tif')
ohy<-raster(path_store&'/'&outnm&'/RASTERS/INIT-COND/ohy.tif')


#===============================================================================
# GET MODELLING DOMAIN SHAPEFILE
#===============================================================================
msk<-st_read(path_gis&'/ModellingRegions/ModellingRegion_'&DOMAIN&'.shp')
msk<-st_zm(msk,drop=TRUE,what="ZM")
msk<-as(msk,'Spatial')
projection(msk)<-projec


##############################################
# DERIVE EFFECTIVE SOIL SATURATION (PER PFT)
#
# Se=(O-Ohy)/(Osat-Ohy)
#   Osat = saturated water content [-]   
#   Ohy = residual or hygroscopic moisture content [-] 

se<-(o-ohy)/(osat-ohy)
se[se < 0]<-0
se[se > 0.8]<-NA ### hack!!!
names(se)<-names(o)
##############################################


#===============================================================================
# READ OBSERVED MAPS & AGGREGATE TO MODEL SPATIAL RESOLUTION
#===============================================================================
##Map from Brun et al., 2020 (GCB)
## Values: 0 = early-witling absence, 1 = early-wilting presence, NA = non-forested area
## Resolution: 10x10 m
## Projection: ETRS89-extended / LAEA Europe (EPSG:3035)
obs<-raster(fn_ewm)
bb2<-spTransform(bb,projection(obs))
obs<-crop(obs,bb2)
f<-res(lai)[1]/res(obs)[1]
obs<-raster::aggregate(obs,fact=f,fun='modal')
obs<-projectRaster(from=obs,crs=projec,method='ngb')
obs<-crop(obs,bb)

plot(obs)

# PLOT REFERENCE EARLY WILTING MAP
r<-obs
r[r < 1]<-NA
pal<-viridis(2,option='H')
tit<-"Presence of early wilting 2018 \n (Brun et al., 2020, GCB)"
fn<-path_fig&'/obs_earlywilting_2018___B2020.png'
png(file=fn,units='in',width=9,height=7.5,res=300)
plot(hs,col=grey(1:100/100),legend=F,main=tit,axes=FALSE)
plot(r,col=pal,add=T,alpha=0.8,legend=F)
plot(msk,add=T)
dev.off()
rm(pal,fn)


#===============================================================================
# PLOT LAI
#===============================================================================
##names(lai)
i1<-15
i2<-29

##########################################################
######### PLOT ABSOLUTE CHANGE (DOYs 195-293)
r<-lai[[i2]]-lai[[i1]]  ##difference in LAI between DOYs 195 and 293
r[is.na(pft)]<-NA                   ##domain mask

pal<-rev(viridis(100,option='H'))
tit<-"Change in LAI between 14 Jul and 20 Oct 2018 [m2/m2] \n (DOYs 195 and 293)"
fn<-path_fig&'/mod_LAI-change_doys195-293_2018.png'
png(file=fn,units='in',width=9,height=7.5,res=300)
plot(hs,col=grey(1:100/100),legend=F,main=tit,axes=FALSE)
plot(r,col=pal,add=T,alpha=0.8,legend=T)
plot(msk,add=T)
dev.off()
rm(r,pal,tit,fn)



##########################################################
######### PLOT RELATIVE CHANGE (DOYs 195-293)
r<-((lai[[i2]]-lai[[i1]])/lai[[i1]])*100  ##relative difference in LAI between DOYs 195 and 293
r[is.na(pft)]<-NA                   ##domain mask

pal<-rev(viridis(100,option='H'))
tit<-"Relative change in LAI between 14 Jul and 20 Oct 2018 [%] \n (DOYs 195 and 293)"
fn<-path_fig&'/mod_LAI-relativechange_doys195-293_2018.png'
png(file=fn,units='in',width=9,height=7.5,res=300)
plot(hs,col=grey(1:100/100),legend=F,main=tit,axes=FALSE)
plot(r,col=pal,add=T,alpha=0.8,legend=T)
plot(msk,add=T)
dev.off()
rm(r,pal,tit,fn)


#===============================================================================
# PLOT EFFECTIVE SOIL SOIL SATURATION
#===============================================================================
##names(se)
i1<-105
i2<-205

##########################################################
######### PLOT ABSOLUTE CHANGE (DOYs 195-295)
r<-se[[i2]]-se[[i1]]  ##difference in Se between DOYs 195 and 295
r[is.na(pft)]<-NA                   ##domain mask

pal<-rev(viridis(100,option='H'))
tit<-"Change in effective soil sat. between 14 Jul and 22 Oct 2018 [-] \n (DOYs 195 and 295)"
fn<-path_fig&'/mod_Se-change_doys195-295_2018.png'
png(file=fn,units='in',width=9,height=7.5,res=300)
plot(hs,col=grey(1:100/100),legend=F,main=tit,axes=FALSE)
plot(r,col=pal,add=T,alpha=0.8,legend=T)
plot(msk,add=T)
dev.off()
rm(r,pal,tit,fn)



##########################################################
######### PLOT RELATIVE CHANGE (DOYs 195-295)
r<-(se[[i2]]-se[[i1]])/se[[i1]]   ##relative difference in Se between DOYs 195 and 295
r[r > 1]<-1
r<-r*100  
r[is.na(pft)]<-NA                   ##domain mask

pal<-rev(viridis(100,option='H'))
tit<-"Relative change in effective soil sat. between 14 Jul and 27 Oct 2018 [%] \n (DOYs 195 and 295)"
fn<-path_fig&'/mod_Se-relativechange_doys195-295_2018.png'
png(file=fn,units='in',width=9,height=7.5,res=300)
plot(hs,col=grey(1:100/100),legend=F,main=tit,axes=FALSE)
plot(r,col=pal,add=T,alpha=0.8,legend=T)
plot(msk,add=T)
dev.off()
rm(r,pal,tit,fn)




#===============================================================================
# PLOT WEEKLY LAI MAPS (JULY-OCT)
#===============================================================================
## names(lai)
## 7 Jul - 27 Oct
r<-lai[[14:30]]
r[is.na(pft)]<-NA                   ##domain mask
mx<-max(values(max(r)),na.rm=T)
mx<-mroundu(mx,0.5)

lab1<-format(as.Date(names(r),format="X%d.%b.%Y"),format="%b %d")
lab2<-format(as.Date(names(r),format="X%d.%b.%Y"),format="%m.%d")
lab2 <- as.numeric(lab2)
r<-stackApply(r,lab2,fun=mean)
names(r) <- lab1
breaks<-seq(0,mx,by=0.5)
# cols<-colorRampPalette(c("red", "yellow", "lightgreen"))(length(breaks)-1)
cols<-viridis(100,option='D')

fn<-path_fig&'/mod_weekly-LAI_2018.png'
png(file=fn,units='in',width=10,height=7,res=300)
rasterVis::levelplot(r,at=breaks,col.regions=cols,main='LAI   [m2/m2]')
dev.off()
rm(r,cols,lab1,lab2,fn)


#===============================================================================
# PLOT WEEKLY Se MAPS (JULY-OCT)
#===============================================================================
## names(se)
## 7 Jul - 27 Oct
r<-se[[seq(98,210,by=7)]]
r[is.na(pft)]<-NA                   ##domain mask
mx<-max(values(max(r)),na.rm=T)
mx<-mroundu(mx,0.1)

lab1<-format(as.Date(names(r),format="X%d.%b.%Y"),format="%b %d")
lab2<-format(as.Date(names(r),format="X%d.%b.%Y"),format="%m.%d")
lab2 <- as.numeric(lab2)
r<-stackApply(r,lab2,fun=mean)
names(r) <- lab1
breaks<-seq(0,mx,by=0.1)
# cols<-colorRampPalette(c("red", "yellow", "lightgreen"))(length(breaks)-1)
cols<-rev(viridis(100,option='C'))

fn<-path_fig&'/mod_weekly-Se_2018.png'
png(file=fn,units='in',width=10,height=7,res=300)
rasterVis::levelplot(r,at=breaks,col.regions=cols,main='Effective soil saturation   [-]')
dev.off()
rm(r,cols,lab1,lab2,fn)




###
print('+++END SCRIPT "PLOT EARLY WILTING"+++   '&outnm)


