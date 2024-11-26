################################################################################
# T&C OUTPUT:
#  READ WEEKLY MAPS & WRITE TO RASTER-BRICKS (.GRD)
#
# 
# TODO: -
# 
# NEW:  
#       
#
#
# 2023/08/29
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
# outnm<-'v20231006_Aaretal_404days'
# outnm<-'v20231006_Aargauer_Rhein_2557days'
# outnm<-'v20231006_Berner_Oberland_779days'
# outnm<-'v20231006_Bodensee_2557days'
outnm<-'v20231006_Engadin_1377days'
# outnm<-'v20231006_Glatt_2557days'
outnm<-'v20231006_Hinterrhein_1807days'
outnm<-'v20231006_Jura_1135days'
outnm<-'v20231006_Leman_1183days'
outnm<-'v20231006_Limmat_1305days'
# outnm<-'v20231006_Neuchatel_813days'
outnm<-'v20231006_Oberwallis_1498days'
outnm<-'v20231006_Reuss_1189days'
# outnm<-'v20231006_Rheintal_402days'
# outnm<-'v20231006_Saane_566days'
# outnm<-'v20231006_Schaffhausen_2557days'
outnm<-'v20231006_Ticino_625days'
# outnm<-'v20231006_Thur_841days'
# outnm<-'v20231006_Toess_1130days'
# outnm<-'v20231006_Unterwallis_376days'
outnm<-'v20231006_Vorderrhein_1122days'
# outnm<-'v20231006_Zuercher_Rhein_2557days'
###
# outnm<-'v20231006_Test-Aletsch_2557days'
# outnm<-'v20231006_Test-Grossbach_2557days'
# outnm<-'v20231006_Test-Poschiavino_2557days'
# outnm<-'v20231106_Test-Rein_da_Sumvitg_2557days'
# outnm<-'v20231006_Test-Rietholzbach_2557days'
# outnm<-'v20231006_Test-Roggiasca_2557days'
# outnm<-'v20231106_Test-Rosegbach_2557days'


#===============================================================================
# PATHS
#===============================================================================
### Root
# root<-'N:/gebhyd/8_Him/Personal_folders/Pascal'
root<-'E:'

## Path to online storage
path_onst<-'C:/Users/Buri/switchdrive/WSL' 

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
library(plotly)
library(doParallel)   #foreach()
library(doSNOW)

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
dir.create(new_dir&'/RASTERS/WEEKLY')


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
# GET MODEL PERIOD DATES
#===============================================================================
# initial timestamp
ts<-grep(pattern="^mod_ini_Date",simf)
ts<-simf[ts]
ts<-str_split(ts,"'")[[1]][2]
date_ini<-as.POSIXct(ts,format='%d-%b-%Y %H:%M:%S')
rm(ts)
# final timestamp
ts<-grep(pattern="^mod_fin_Date",simf)
ts<-simf[ts]
ts<-str_split(ts,"'")[[1]][2]
date_fin<-as.POSIXct(ts,format='%d-%b-%Y %H:%M:%S')
rm(ts)
Date<-seq(date_ini,date_fin,by=3600)


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


#===============================================================================
# ADD VARIABLE NAMES & CORRESPONDING UNITS
#===============================================================================
setwd(path_code)
setwd('..')
var_df<-read.table('TC_Variables.txt',header=T,sep='\t')
var_df[]<-lapply(var_df,as.character)


#===============================================================================
# READ HEADERS OF SPATIAL OUTPUT MAPS
#===============================================================================
rawfile<-paramfile
# find lines where output variables are defined
idx1<-grep(pattern="START <<OUTPUT_WE_ZZZ_SPATIAL_YYY.dat>>",x=rawfile)
idx2<-grep(pattern="END <<OUTPUT_WE_ZZZ_SPATIAL_YYY.dat>>",x=rawfile)
if(length(idx1) == 0 | length(idx2) == 0){stop('Lines in "paramfile" not found!')}
rawfile<-rawfile[idx1:idx2]
# find all variable-lines
vars<-grepl(pattern="'",x=rawfile)
rawfile<-rawfile[vars]
vars<-lapply(strsplit(rawfile,','), function(x) x[1])
vars<-do.call('rbind',vars)
vars<-lapply(strsplit(vars,"'"), function(x) x[2])
vars<-do.call('c',vars)
vars<-vars[endsWith(vars,'_we')] ##drop non-relevant variables
vars<-str_remove(vars,'_spatial_we')
rm(idx1,idx2,rawfile)


#===============================================================================
# READ SPATIAL OUTPUT MAPS
#===============================================================================
# get all output files, put in temporal order
setwd(path_data&'/OUTPUTS')
fn_v<-list.files(pattern='OUTPUT_WE_'&simnm&'_SPATIAL_')
if(reinit){
  fn_v<-c(list.files(pattern='OUTPUT_WE_'&simnm1&'_SPATIAL_'),fn_v)
}
iTS<-stringr::str_remove(fn_v,'.mat$')
iTS<-as.numeric(sub('.*_SPATIAL_','',iTS))#extract intermediate timesteps
idx<-base::order(iTS)  #ordering index
iTS<-iTS[idx]
fn_v<-fn_v[idx]
iDate<-format(Date[iTS],'%d %b %Y')
PERIOD_LO<-iDate[1]&' - '&iDate[length(iDate)]
PERIOD_LO<-paste0(str_split(PERIOD_LO,' ')[[1]][-c(1,5)],collapse=' ')
PERIOD_SH<-paste0(str_split(PERIOD_LO,' ')[[1]],collapse='')

rmidx<-which(is.na(iDate)) ##remove "empty" months (e.g. end of simulations)
if(length(rmidx != 0)){
  iTS<-iTS[-rmidx]
  fn_v<-fn_v[-rmidx] 
  iDate<-iDate[-rmidx]
}


###### LOOP OVER VARIABLES #######
n<-length(fn_v)
pb<-txtProgressBar(max=n,style=3)
progr<-function(n) setTxtProgressBar(pb,n)
opts<-list(progress=progr)
pkgs<-c('R.matlab','matlab','raster')   #packages needed in MAT2RAST()

for(v in 1:length(vars)){
  VAR<-vars[v]
  VAR2<-str_replace(VAR,'_','.')
  cl<-makeCluster(2);  registerDoSNOW(cl)
  print(paste0('*** start: ',Sys.time(),' ***********************************'))
  print('*** '&DOMAIN&': Processing weekly << '&VAR&' >> outputs ('&n&' weeks)  |||  Variable '&v&' of '&length(vars)&' ***')
  LS<-foreach(i = icount(n),.combine=c,.options.snow=opts,.packages=pkgs) %dopar% {
    MAT2RAST(fn_v[i],paste0(VAR2,'.spatial.we'),m_cell,n_cell)
  }
  close(pb); stopCluster(cl) 
  br<-raster::brick(LS)
  rm(LS)
  print(paste0('************************************ done: ',Sys.time(),' ***'))
  print('***')
  names(br)<-iDate
  writeRaster(br,file=new_dir&'/RASTERS/WEEKLY/'&VAR&'___'&PERIOD_SH&'___br.grd',overwrite=T)
  rm(VAR,VAR2,br)
}



###
print('+++END SCRIPT "READ WEEKLY MAPS"+++   '&outnm)









