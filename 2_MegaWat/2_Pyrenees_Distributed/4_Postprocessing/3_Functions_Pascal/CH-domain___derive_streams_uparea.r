################################################################################
# T&C ENTIRE CH DOMAIN:
#  DERIVE STREAMNETWORK UPSTREAM AREAS
#
# 
# TODO: -
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
##define output folder name (=main folder name)
outnm<-'CH___v20240117'


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

## Path to postprocessed data
path_data<-root&'/T&C/POSTPROCESSED/DISTRIBUTED/'&outnm

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
library(Rsagacmd)
library(terra)
library(sf)
library(reshape2)
library(ggplot2)
library(stringr)      #str_replace()
library(wesanderson)  #wes_palette()
library(lutz)         #tz_lookup()
library(viridis)      #scale_color_viridis()
library(ggsci)        #scale_color_jco()
library(scales)       #date_format()
library(rasterVis)    #levelplot
library(gridExtra)    #grid.arrange()
library(RColorBrewer) #brewer.pal()
library(lubridate)    #floor_date(), hour(), %m+% 
library(extrafont)
library(grid)
library(gridtext)
library(xtable)
library(jpeg)
library(rworldxtra)
library(tcltk2)        #tkProgressBar()
library(doParallel)   #foreach()
library(doSNOW)

###font_import()
loadfonts(device = "win")
# define no. of digits printed in console
options('scipen'=100,'digits'=4,warn=1)

## define corresponding timezone
# TZ<-'Europe/Zurich'  #considers summer time (CET/CEST)
TZ<-'Africa/Algiers'  #does not have daylight saving

# set time zone & change local time from german to english for plotting
# (check in sessionInfo() ) 
Sys.setenv(TZ=TZ)
Sys.setlocale('LC_TIME', 'C')

# register cores for parallel functionality
# registerDoParallel(cores=n.cores)   #put "stopImplicitCluster()" at the end

## plot specs
source(path_func)
base_size<-15
base_family<-'Arial'


#===============================================================================
# EXTRACT DOMAIN NAME
#===============================================================================
DOMAIN<-str_split(outnm,'_')[[1]][1]


#===============================================================================
# READ INIT. COND. RASTERS
#===============================================================================
dem<-terra::rast(path_data&'/RASTERS/INIT-COND/DTM.tif')   ## DEM
dem_m<-as.matrix(dem)
slo<-terra::terrain(dem,v='slope',unit='degrees')  ## slope
asp<-terra::terrain(dem,v='aspect',unit='degrees') ## aspect
fd<-terra::terrain(dem,v='flowdir') ## flow direction
sn<-terra::rast(path_data&'/RASTERS/INIT-COND/SN.tif')     ## streamnetwork
lc<-terra::rast(path_data&'/RASTERS/INIT-COND/VEG_CODE.tif')   ## LC
sn[sn == 0]<-NA


#===============================================================================
# MASK STREAMNETWORK WITH LAKES/WATER
#===============================================================================
#   #1) Evergreen Tree Forest
#   #2) Mixed Evergreen (50%) Decidous (50%)
#   #3) Decidous
#   #4) Meadow - Grasses
#   #5) Rock
#   #6) Water
#   #7) Ice (bedrock)
sn[lc == 6]<-NA


#==================================================================================
# DEFINE POURPOINTS
#==================================================================================
ppts<-as.data.frame(crds(sn,na.rm=T))
ppts<-vect(ppts,geom=c("x","y"))
##### select subset ######
ppts<-ppts[seq(1,length(ppts),by=50)]
##########################
plot(dem)
plot(ppts, add=T)


#===================================================
# CONNECT TO SAGA-GIS APPLICATION
# https://github.com/stevenpawley/Rsagacmd
#===================================================
if(!exists('saga')){
  saga<-saga_gis(saga_bin='C:/Program Files/SAGA/saga_cmd.exe')
}


#===================================================
# REMOVE SINKS IN DEM
# https://saga-gis.sourceforge.io/saga_tool_doc/9.2.0/ta_preprocessor_3.html
#===================================================
if(!exists('dem_filled')){
  dem_filled<-saga$ta_preprocessor$fill_sinks_planchon_darboux_2001(dem=dem,
                                                                    result='filled_DEM.sgrd')
}


new_dir<-path_data&'/RASTERS/UPAREA'
dir.create(new_dir)
setwd(new_dir)
#==================================================================================
# DERIVE UPSTREAM AREA PER POURPOINT (FOR LOOP)
#==================================================================================
## define method to determine flow directions
meth<-0 # Deterministic 8
# meth<-1 # Deterministic Infinity
# meth<-2 # Multiple Flow Direction
# meth<-3 # Multiple Triangular Flow Directon
# meth<-4 # Multiple Maximum Downslope Gradient Based Flow Directon

df<-terra::extract(dem,ppts,ID=T,xy=T,cells=T)
df<-round(df,0)
df$ncell<-NA
rm(ppts)
n<-nrow(df)
pb<-txtProgressBar(max=n,style=3)
r_ls<-list()
print(paste0('*** start: ',Sys.time(),' ***********************************'))
print('*** '&DOMAIN&': Processing upstream area for '&n&' pourpoints ***')
for(i in 1:n){
# for(i in c(111,222,333,444,555,666)){    ##testing
  setTxtProgressBar(pb,i)
  r<-UPAREA(df[i,'x'],df[i,'y'],dem_filled,meth)
  df$ncell[i]<-terra::freq(r,value=1)$count
  r_ls[[i]]<-r
  rm(r)
}
close(pb)
file.remove(dir(new_dir,pattern="^ua.",full.names=TRUE)) ##delete temp. files
br<-rast(r_ls)
rm(r_ls)
print(paste0('************************************ done: ',Sys.time(),' ***'))
print('***')

names(br)<-'id'&df$ID&'_'&df$DTM&'m'
write.table(df,new_dir&'/uparea_id___n'&n&'.check',dec='.',sep='\t',quote=FALSE,
            col.names=TRUE,row.names=FALSE)
writeRaster(br,file=new_dir&'/uparea_br___n'&n&'.grd',overwrite=T)
rm(df,br)


##alternative: 
# WriteNCDF {greenbrown} ?
##OR:
# test <- brick('/TavgM_1981.gri')
# writeRaster(test, "rstack.nc", overwrite=TRUE, format="CDF",     varname="Temperature", varunit="degC", 
#             longname="Temperature -- raster stack to netCDF", xname="Longitude",   yname="Latitude", zname="Time (Month)")

### r<-brick(new_dir&'/uparea_br___n'&n&'.grd')


# #==================================================================================
# # DERIVE UPSTREAM AREA PER POURPOINT (PARALLEL LOOP)
# #==================================================================================
# ## define method to determine flow directions
# meth<-0 # Deterministic 8
# # meth<-1 # Deterministic Infinity
# # meth<-2 # Multiple Flow Direction
# # meth<-3 # Multiple Triangular Flow Directon
# # meth<-4 # Multiple Maximum Downslope Gradient Based Flow Directon
# 
# ## prepare parallel loop
# xy<-crds(ppts)
# n<-length(ppts)
# pb<-txtProgressBar(max=n,style=3)
# progr<-function(n) setTxtProgressBar(pb,n)
# opts<-list(progress=progr)
# pkgs<-c('Rsagacmd','terra')   #packages needed in UPAREA()
# cl<-makeCluster(2);  registerDoSNOW(cl)
# 
# ## parallel loop
# print(paste0('*** start: ',Sys.time(),' ***********************************'))
# print('*** '&DOMAIN&': Processing upstream area for '&n&' pourpoints ***')
# LS<-foreach(i = icount(n),.combine=c,.options.snow=opts,.packages=pkgs) %dopar% {
#   UPAREA(xy[i,],dem_filled,meth)
# }
# close(pb); stopCluster(cl) 
# br<-raster::brick(LS)
# rm(LS)
# print(paste0('************************************ done: ',Sys.time(),' ***'))
# print('***')
# names(br)<-iDate
# writeRaster(br,file=new_dir&'/RASTERS/uparea_br.grd',overwrite=T)
# rm(br)
# stopImplicitCluster()



###
print('+++END SCRIPT <<DERIVE STREAM UPAREAS>>+++   '&outnm)



