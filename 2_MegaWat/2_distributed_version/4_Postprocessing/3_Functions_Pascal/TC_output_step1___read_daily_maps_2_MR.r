################################################################################
# T&C OUTPUT:
#  READ DAILY MAPS & WRITE TO RASTER-BRICKS (.GRD)
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
# outnm<-'v20231213_Aaretal_2557days'
# outnm<-'v20231213_Aargauer_Rhein_2557days'
# outnm<-'v20231213_Berner_Oberland_2557days'
# outnm<-'v20231213_Bodensee_2557days'
# outnm<-'v20231213_Engadin_2557days'
outnm<-'v20241115_Cinca_Mid_XXXXdays'
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
### Root
# root<-'N:/gebhyd/8_Him/Personal_folders/Pascal'
root<-'M:/19_ISTA/1_TC/3_Model_Source/'

## Path to online storage
path_onst<-'C:/Users/Buri/switchdrive/WSL' 

## Path to stored data
path_data<-root&'3_Output_files/'&outnm

## Path to input data
path_input<-root&'2_TeC_Source_Code/European_catchments/2_Catchments/Pyrenees_distributed/Inputs/'

## Path to store data
path_store<-root&'3_Output_files/2_Postprocessing_output/'

## Path to code scripts etc.
path_code<-path_onst&'/T&C/Code/POSTprocessing/CH2022'

## Path to functions file
path_func<-root&'2_TeC_Source_Code/European_catchments/1_Functions/2_Functions/1_Helper_functions/Functions.r'

## Path to GIS data
path_gis<-root&'/QGIS'


#===============================================================================
# LOAD PACKAGES
#===============================================================================
# define packages to load
library(raster)
library(sf)
#library(rgdal) --> Old function. Seek a replacement in packages sp or sf         
#library(rgeos) --> Old function. Seek a replacement in packages sp or sf
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
#library(maptools)     #unionSpatialPolygons() --> Old function. Seek a replacement in packages sp or sf
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
dir.create(new_dir&'/RASTERS/DAILY')


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

# ### BASIN MASK
msk<-geodatafile[['MASK']]
msk[msk < 1]<-NA
msk<-apply(msk,2,rev)   ## Matlab: flipud(DTM)
#check
# plot(raster(msk))
msk<-raster(nrow=m_cell,ncol=n_cell,vals=msk)
projection(msk)<-as.character(projec)


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
idx1<-grep(pattern="START <<OUTPUT_DA_ZZZ_SPATIAL_YYY.dat>>",x=rawfile)
idx2<-grep(pattern="END <<OUTPUT_DA_ZZZ_SPATIAL_YYY.dat>>",x=rawfile)
if(length(idx1) == 0 | length(idx2) == 0){stop('Lines in "paramfile" not found!')}
rawfile<-rawfile[idx1:idx2]
# find all variable-lines
vars<-grepl(pattern="'",x=rawfile)
rawfile<-rawfile[vars]
vars<-lapply(strsplit(rawfile,','), function(x) x[1])
vars<-do.call('rbind',vars)
vars<-lapply(strsplit(vars,"'"), function(x) x[2])
vars<-do.call('c',vars)
vars<-vars[endsWith(vars,'_da')] ##drop non-relevant variables
vars<-str_remove(vars,'_spatial_da')
rm(idx1,idx2,rawfile)


#===============================================================================
# READ SPATIAL OUTPUT MAPS
#===============================================================================
# get all output files, put in temporal order
setwd(path_data&'/OUTPUTS')
fn_v<-list.files(pattern='OUTPUT_DA_'&simnm&'_SPATIAL_')
if(reinit){
  fn_v<-c(list.files(pattern='OUTPUT_DA_'&simnm1&'_SPATIAL_'),fn_v)
}
iTS<-stringr::str_remove(fn_v,'.mat')
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
  print('*** '&DOMAIN&': Processing daily << '&VAR&' >> outputs ('&n&' days)  |||  Variable '&v&' of '&length(vars)&' ***')
  LS<-foreach(i = icount(n),.combine=c,.options.snow=opts,.packages=pkgs) %dopar% {
    MAT2RAST(fn_v[i],paste0(VAR2,'.spatial.da'),m_cell,n_cell)
  }
  close(pb); stopCluster(cl) 
  br<-raster::brick(LS)
  rm(LS)
  print(paste0('************************************ done: ',Sys.time(),' ***'))
  print('***')
  names(br)<-iDate
  writeRaster(br,file=new_dir&'/RASTERS/DAILY/'&VAR&'___'&PERIOD_SH&'___br.grd',overwrite=T)
  rm(VAR,VAR2,br)
}




###
print('+++END SCRIPT "READ DAILY MAPS"+++   '&outnm)










