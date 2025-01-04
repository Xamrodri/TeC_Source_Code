################################################################################
# T&C ENTIRE CH DOMAIN:
#  READ INITIAL CONDITIONS MAPS & WRITE TO TIFF-FILES
#
# 
# TODO: -
# 
# NEW:  
#       
#
#
# 2024/01/17
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
##define input folder name (=main folder name)
inputfo<-'GeodataInput_CH___v20240117'

##define CH-boundaries filename
fnBord = 'SwissTLM3D/CH_boundaries_EPSG-2056.shp';


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

# access functions file
source(path_func)


#===============================================================================
# CREATE NEW FOLDERS
#===============================================================================
outnm<-str_split(inputfo,'GeodataInput_')[[1]][2]
new_dir<-path_store&'/'&outnm
dir.create(new_dir)
dir.create(new_dir&'/RASTERS')
dir.create(new_dir&'/RASTERS/INIT-COND')
dir.create(new_dir&'/SHAPEFILES')
dir.create(new_dir&'/PLOTS')
dir.create(new_dir&'/PLOTS/INIT-COND')


#===============================================================================
# LOAD GEODATA & GEOREF INPUT-FILES << Geodata_XXX.mat >> & << GEOREF_XXX.mat >>
#===============================================================================
setwd(path_input&'/'&inputfo)
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
# READ CH-BORDERS & WRITE TO FILE
#==================================================================================
bord<-st_read(path_gis&'/'&fnBord,quiet=T)
setwd(new_dir&'/SHAPEFILES')
bord2<-bord
bord<-as(bord,'Spatial')

##write to file
st_write(bord2,'msk.shp',
         driver="ESRI Shapefile",delete_dsn=TRUE)


#===============================================================================
# WRITE RASTERS FROM GEODATA-INPUT TO TIFF-FILES
#===============================================================================
# > names(geodatafile)
# [1] "DEB.MAP"   "DTM"       "GLA.MAP2"  "GLH"       "MASK"      "PCLA"      "PORG"     
# [8] "PSAN"      "SN"        "SNOWD"     "SOIL.TH"   "SWEQ"      "VEG.CODE"  "cellsize" 
# [15] "m"         "n"         "x"         "xllcorner" "y"         "yllcorner"
vars<-c("DEB.MAP","DTM","GLH",
        "SN","PCLA","PORG","PSAN",
        "SNOWD","SOIL.TH")
tit_v<-c("Debris thickness","Elevation","Glacier ice thickness",
         "Stream network","Clay content (soil)","Organic content (soil)","Sand content (soil)",
         "Snow depth","Soil depth")
un_v<-c("[mm]","[m a.s.l.]","[m]","[0,1]","[%]","[%]","[%]","[m]","[m]")


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


