################################################################################
# T&C OUTPUT:
#  READ INITIAL CONDITIONS MAPS & WRITE TO TIFF-FILES
#
# 
# TODO: -
# 
# NEW:  
#       
#
#
# 2023/10/06
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
##define number of cores to be used
# n.cores<-parallel::detectCores(all.tests = FALSE, logical = TRUE) - 3
# n.cores<-parallel::detectCores(all.tests = FALSE, logical = TRUE)
n.cores<-1  

## define output name (=main folder name)
# outnm<-'v20231213_Aaretal_2557days'
# outnm<-'v20231213_Aargauer_Rhein_2557days'
# outnm<-'v20231213_Berner_Oberland_2557days'
# outnm<-'v20231213_Bodensee_2557days'
# outnm<-'v20231213_Engadin_2557days'
outnm<-'v20231213_Glatt_2557days'
outnm<-'v20231213_Hinterrhein_2557days'
# outnm<-'v20231213_Jura_2557days'
# outnm<-'v20231213_Leman_2557days'
# outnm<-'v20231213_Limmat_2557days'
# outnm<-'v20231213_Neuchatel_2557days'
# outnm<-'v20231213_Oberwallis_2557days'
# outnm<-'v20231213_Reuss_2557days'
# outnm<-'v20231213_Rheintal_2557days'
# outnm<-'v20231213_Saane_2557days'
# outnm<-'v20231213_Schaffhausen_2557days'
# outnm<-'v20231213_Ticino_2557days'
# outnm<-'v20231213_Thur_2557days'
# outnm<-'v20231213_Toess_2557days'
# outnm<-'v20231213_Unterwallis_2557days'
# outnm<-'v20231213_Vorderrhein_2557days'
# outnm<-'v20231213_Zuercher_Rhein_2557days'
###
# outnm<-'v20231213_Test-Aletsch_2557days'
# outnm<-'v20231211_Test-Grossbach_2557days'
# outnm<-'v20231213_Test-Poschiavino_2557days'
# outnm<-'v20231213_Test-Rein_da_Sumvitg_2557days'
# outnm<-'v20231213_Test-Rietholzbach_2557days'
# outnm<-'v20231213_Test-Roggiasca_2557days'
# outnm<-'v20231211_Test-Rosegbach_2557days'


#===============================================================================
# PATHS
#===============================================================================
#X1
if(Sys.getenv('computername') == 'PASCAL-THINK'){ 
  root<-'D:'                                    ## Root
  path_onst<-'C:/Users/Pascal/switchdrive/WSL'  ## Path to online storage
}
#T490
if(Sys.getenv('computername') == 'WSL26871'){ 
  root<-'E:'                                   ## Root
  path_onst<-'C:/Users/Buri/switchdrive/WSL'   ## Path to online storage
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
path_gis<-root&'/QGIS'


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
library(doParallel)   #foreach()
library(doSNOW)
library(rasterVis)
library(RColorBrewer)
library(viridisLite)
library(extrafont)
# font_import()
loadfonts(device = "win")

# define no. of digits printed in console
options('scipen'=100, 'digits'=4)

## define corresponding timezone
# TZ<-'Europe/Zurich'  #considers summer time (CET/CEST)
TZ<-'Africa/Algiers'  #does not have daylight saving

# set time zone & change local time from german to english for plotting
# (check in sessionInfo() ) 
Sys.setenv(TZ=TZ)
Sys.setlocale('LC_TIME', 'C')

# register cores for parallel functionality
registerDoParallel(cores=n.cores)   #put "stopImplicitCluster()" at the end

# access functions file
source(path_func)


#===============================================================================
# CREATE NEW FOLDERS
#===============================================================================
new_dir<-path_store&'/'&outnm
dir.create(new_dir)
dir.create(new_dir&'/RASTERS')
dir.create(new_dir&'/RASTERS/INIT-COND')
dir.create(new_dir&'/SHAPEFILES')
dir.create(new_dir&'/PLOTS')
dir.create(new_dir&'/PLOTS/INIT-COND')


#===============================================================================
# READ PARAMETER FILE << OUTPUT_MANAGER_PAR.m >>
#===============================================================================
setwd(path_data)
paramfile<-readLines('OUTPUT_MANAGER_PAR.m')
## if you want to match metacharacters (e.g. parentheses etc.) 
## you have to "scape" the metacharacter with two backslashes


#===============================================================================
# READ SIMULATION FILE << AAA_TC_XXX.m >> OR << BBB_TC_XXX_reinit.m >>
#===============================================================================
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
# GET INPUT FOLDER NAME
#===============================================================================
idx<-grep(pattern="^GEODATAINPUT = \\[",simf)
fo<-simf[idx]
fo<-str_split(fo,"___v")[[1]][2]
inputfo<-path_input&'/GeodataInput_'&DOMAIN&'___v'&as.character(gsub("\\D","",fo))
rm(idx,fo)


#===============================================================================
# LOAD GEODATA & GEOREF INPUT-FILES << Geodata_XXX.mat >> & << GEOREF_XXX.mat >>
#===============================================================================
setwd(inputfo)
files<-list.files()
geodatafile<-readMat(files[which(str_starts(files,'Geodata_'))])

### RASTER DIMENSIONS & coordinates
m_cell<-geodatafile[['m']][1,1]
n_cell<-geodatafile[['n']][1,1]
resol<-geodatafile[['cellsize']][1,1]
x<-geodatafile[['x']][,1]
y<-geodatafile[['y']][,1]

### PROJECTION
georef<-readLines(files[which(str_starts(files,'GEOREF_'))])
epsg<-str_split(georef[[1]][1],'-')[[1]]
epsg<-as.numeric(epsg[length(epsg)])
projec<-make_EPSG()
projec<-projec[which(projec$code == epsg),'prj4']

### EXTENT
ext<-c(min(x)-resol/2,max(x)+resol/2,min(y)-resol/2,max(y)+resol/2)

### MASK
msk<-geodatafile[['MASK']]
msk<-matlab::reshape(msk,m_cell,n_cell) 
msk<-raster::raster(apply(msk,2,rev))
raster::projection(msk)<-as.character(projec)
raster::extent(msk)<-ext
writeRaster(msk,new_dir&'/RASTERS/INIT-COND/msk.tif','GTiff',overwrite=TRUE)

##define land cover classes
lc_code<-c("Evergreen forest","Mixed forest","Decidous forest","Grass",'Rock',"Water","Glacier ice")


#==================================================================================
# GENERATE SHAPEFILE FROM CATCHMENT-RASTER & WRITE TO FILE
#==================================================================================
setwd(new_dir&'/SHAPEFILES')
fn<-'msk.shp'
if(file.exists(fn)){
  basin<-st_read(fn,quiet=T)
  basin2<-basin
  basin<-as(basin,'Spatial')
}
if(!file.exists(fn)){
  xy<-coordinates(msk)
  idx<-which(!is.na(as.vector(msk)))
  xy<-xy[idx,]
  pts<-SpatialPoints(cbind(xy[,'x'],xy[,'y']))
  projection(pts)<-projection(msk)
  buf<-gBuffer(pts,byid=FALSE,width=res(msk)[1])  #expensive!
  basin<-gBuffer(buf,byid=FALSE,width=-res(msk)[1]*0.5)
  basin<-unionSpatialPolygons(basin,rep(1,length(basin)))
  rm(xy,idx,buf,pts)
  basin2<-st_as_sf(basin)
  ##write to file
  st_write(basin2,fn,
           driver="ESRI Shapefile",delete_dsn=TRUE)
}


setwd(path_data&'/OUTPUTS')
#===============================================================================
# WRITE & PLOT HYGROSCOPIC & SATURATION WATER CONTENT
#===============================================================================
files<-list.files()

### OHY MAP
ohyfile<-readMat(files[which(str_ends(files,'Ohy_1.mat'))])
ohy<-ohyfile[['Ohy.OUT']]
ohy<-apply(ohy,2,rev)
ohy<-raster(nrow=m_cell,ncol=n_cell,vals=ohy)
projection(ohy)<-as.character(projec)
extent(ohy)<-ext
ohy[msk == 0]<-NA
writeRaster(ohy,new_dir&'/RASTERS/INIT-COND/ohy.tif','GTiff',overwrite=TRUE)
png(file=new_dir&'/PLOTS/INIT-COND/ohy.png',units='in',width=5,height=5,res=300,family='Arial')
plot(ohy,axes=FALSE,box=FALSE,main='Hygroscopic moisture content',cex.main=1.5,legend=F,
     axis.args=list(cex.axis=1.5))
plot(ohy,legend.only=TRUE,
     legend.width=1,legend.shrink=0.75,
     axis.args=list(
       cex.axis=1.5
     ),
     legend.args=list(text='[-]',side=3,font=1,line=0.5,cex=1.5))
dev.off()

### OSAT MAP
osatfile<-readMat(files[which(str_ends(files,'Osat_1.mat'))])
osat<-osatfile[['Osat.OUT']]
osat<-apply(osat,2,rev)
osat<-raster(nrow=m_cell,ncol=n_cell,vals=osat)
projection(osat)<-as.character(projec)
extent(osat)<-ext
osat[msk == 0]<-NA
writeRaster(osat,new_dir&'/RASTERS/INIT-COND/osat.tif','GTiff',overwrite=TRUE)
png(file=new_dir&'/PLOTS/INIT-COND/osat.png',units='in',width=5,height=5,res=300,family='Arial')
plot(osat,axes=FALSE,box=FALSE,main='Saturated water content',cex.main=1.5,legend=F,
     axis.args=list(cex.axis=1.5))
plot(osat,legend.only=TRUE,
     legend.width=1,legend.shrink=0.75,
     axis.args=list(
       cex.axis=1.5
     ),
     legend.args=list(text='[-]',side=3,font=1,line=0.5,cex=1.5))
dev.off()


#===============================================================================
# WRITE RASTERS FROM GEODATA-INPUT TO TIFF-FILES
#===============================================================================
# > names(geodatafile)
# [1] "Aacc"          "DEB.MAP"       "DTM"           "DTM.orig"      "GLA.MAP"      
# [6] "GLH"           "MASK"          "PCLA"          "MCOS"          ""             
# [11] "PORG"          "PSAN"          "SN"            "SNOWD"         "SOIL.TH"      
# [16] "T.flow"        "VEG.CODE"      "Xout"          "Xoutlet"       "Yout"         
# [21] "Youtlet"       "cellsize"      "m"             "n"             "region.ln"    
# [26] "x"             "x.all"         "xllcorner"     "xllcorner.all" "y"            
# [31] "y.all"         "yllcorner"     "yllcorner.all" "" 
vars<-c("Aacc","DEB.MAP","DTM","GLH","SN","PCLA","PORG","PSAN","SNOWD","SOIL.TH")
tit_v<-c("log(upslope Area)","Debris thickness","Elevation","Glacier ice thickness",
         "Stream network","Clay content (soil)","Organic content (soil)","Sand content (soil)",
         "Snow depth","Soil depth")
un_v<-c("[cells]","[mm]","[m a.s.l.]","[m]","[0,1]","[%]","[%]","[%]","[m]","[m]")


for(v in 1:length(vars)){
  VAR<-vars[v]
  VAR2<-str_replace(VAR,'\\.','\\_')
  mat<-geodatafile[[VAR]]
  mat<-matlab::reshape(mat,m_cell,n_cell) 
  r<-raster::raster(apply(mat,2,rev))
  rm(mat)
  raster::projection(r)<-as.character(projec)
  raster::extent(r)<-ext
  r[msk == 0]<-NA
  
  # if(VAR %in% c("DEB.MAP","GLH","PCLA","PORG","PSAN","SNOWD","SOIL.TH")){r[r == 0]<-NA}
  cols<-viridis(100)
  if(VAR == 'Aacc'){r<-log(r)}
  
  writeRaster(r,new_dir&'/RASTERS/INIT-COND/'&VAR2&'.tif','GTiff',overwrite=TRUE)
  
  png(file=new_dir&'/PLOTS/INIT-COND/'&VAR2&'.png',units='in',width=6,height=6,res=300,family='Arial')
  plot(r,axes=FALSE,box=FALSE,main=tit_v[v],cex.main=1.5,legend=F,
       axis.args=list(cex.axis=1.5),col=cols)
  plot(r,col=cols,legend.only=TRUE,
       legend.width=1,legend.shrink=0.75,
       axis.args=list(
         cex.axis=1.5
       ),
       legend.args=list(text=un_v[v],side=3,font=1,line=0.5,cex=1.5))
  dev.off()
  rm(VAR,VAR2,r)
}


## LAND COVER (NOT WORKING IN LOOP!)
VAR<-"VEG.CODE"
VAR2<-str_replace(VAR,'\\.','\\_')
mat<-geodatafile[[VAR]]
mat<-matlab::reshape(mat,m_cell,n_cell) 
r<-raster::raster(apply(mat,2,rev))
rm(mat)
raster::projection(r)<-as.character(projec)
raster::extent(r)<-ext
r[msk == 0]<-NA
### LAND COVER (CONVERT TO CATEGORIAL RASTER)
r<-as.factor(r)# Values as factor
rat<-levels(r)[[1]]# Extract attribute table
rat[["Class"]]<-lc_code[rat$ID]# Set custom breaks
levels(r)<-rat# Add back RAT

writeRaster(r,new_dir&'/RASTERS/INIT-COND/'&VAR2&'.tif','GTiff',overwrite=TRUE)

# cols<-colorRampPalette(brewer.pal(8,"Dark2"))(length(lc_code))
cols<-c('#008B00','#ADFF2F','#C8611F','#BBA90B','#778899','#3A5FCD','#87CEFA')
cols<-cols[levels(r)[[1]]$ID] ##pick colours
png(file=new_dir&'/PLOTS/INIT-COND/'&VAR2&'.png',units='in',width=6,height=6,res=300,family='Arial')
rasterVis::levelplot(r,main="Land cover",col.regions=cols,
                     par.settings = list(axis.line = list(col = "transparent"), 
                                         strip.background = list(col = 'transparent'), 
                                         strip.border = list(col = 'transparent')))
dev.off()
rm(VAR,VAR2,r)



###
print('+++END SCRIPT <<READ/PLOT INIT COND>>+++   '&outnm)


