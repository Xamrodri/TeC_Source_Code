################################################################################
# T&C OUTPUT STEP 2:
#  DERIVE STREAMNETWORK UPSTREAM AREAS
#
# 
# TODO: -
#
#
# 2023/12/14
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
# #define number of cores to be used
# n.cores<-parallel::detectCores(all.tests = FALSE, logical = TRUE) - 3
# n.cores<-parallel::detectCores(all.tests = FALSE, logical = TRUE)
# n.cores<-2

## define output name (=main folder name)
# outnm<-'v20231213_Aaretal_2557days'
# outnm<-'v20231213_Aargauer_Rhein_2557days'
# outnm<-'v20231213_Berner_Oberland_2557days'
# outnm<-'v20231213_Bodensee_2557days'
# outnm<-'v20231213_Engadin_2557days'
outnm<-'v20231213_Glatt_2557days'
# outnm<-'v20231213_Hinterrhein_2557days'
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
outnm<-'v20231213_Toess_2557days'
# outnm<-'v20231213_Unterwallis_2557days'
outnm<-'v20231213_Vorderrhein_2557days'
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
DOMAIN<-str_split(outnm,'_')[[1]]
DOMAIN<-DOMAIN[-c(1,length(DOMAIN))]


#===============================================================================
# READ INIT. COND. RASTERS
#===============================================================================
dem<-terra::rast(path_data&'/RASTERS/INIT-COND/DTM.tif')   ## DEM
dem_m<-as.matrix(dem)
slo<-terra::terrain(dem,v='slope',unit='degrees')  ## slope
asp<-terra::terrain(dem,v='aspect',unit='degrees') ## aspect
fd<-terra::terrain(dem,v='flowdir') ## flow direction
acc<-terra::rast(path_data&'/RASTERS/INIT-COND/Aacc.tif')  ## upstream accumulation area
sn<-terra::rast(path_data&'/RASTERS/INIT-COND/SN.tif')     ## streamnetwork
sn[sn == 0]<-NA


#==================================================================================
# DEFINE POURPOINTS
#==================================================================================
ppts<-as.data.frame(crds(sn,na.rm=T))
ppts<-vect(ppts,geom=c("x","y"))
##### select subset ######
ppts<-ppts[seq(1,length(ppts),by=3)]
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



