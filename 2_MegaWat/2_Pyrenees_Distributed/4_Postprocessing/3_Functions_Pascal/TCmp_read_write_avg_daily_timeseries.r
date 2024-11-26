################################################################################
# T&C-MULTIPOINT OUTPUT:
#  READ HOURLY & DAILY TIMESERIES, AVERAGE TO DAILY AND WRITE TO FILE
#
# 
# TODO: -
# 
# NEW:  -drop soil layers below soil thickness or entirely if no soil present
#       
#
#
# 2023/08/17
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
outnm<-'v20230727_SP'
# outnm<-'v20230727_FLUX'


#===============================================================================
# PATHS
#===============================================================================
### Root
# root<-'N:/gebhyd/8_Him/Personal_folders/Pascal'
root<-'E:'

## Path to online storage
path_onst<-'C:/Users/Buri/switchdrive/WSL' 

## Path to stored data
path_data<-root&'/T&C/OUTPUTSTORAGE/MULTIPOINT/'&outnm

## Path to input data
path_input<-root&'/T&C/INPUTS'

## Path to store data
path_store<-root&'/T&C/POSTPROCESSED/MULTIPOINT'

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
# ADD VARIABLE NAMES & CORRESPONDING UNITS
#===============================================================================
setwd(path_code)
setwd('..')
var_df<-read.table('TCmp_Variables.txt',header=T,sep='\t')
var_df[]<-lapply(var_df,as.character)


#===============================================================================
# READ SIMULATION FILE << AAA_TCmp_XXX.m >>
#===============================================================================
setwd(path_data)
files<-list.files()
##read AAA-file
simfile<-files[startsWith(files,'AAA_')]
simfile<-readLines(simfile,warn=FALSE)


#===============================================================================
# READ PARAMETERS FILE << PARAMETERS_SOIL_XXX.m >>
#===============================================================================
paramfile<-grep(pattern="^fnParameter",simfile)
paramfile<-simfile[paramfile]
paramfile<-str_split(paramfile,"'")[[1]][2]
paramfile<-readLines(path_input&'/'&paramfile,warn=FALSE)

# find line where Zs (soil layer depths) is defined
rawfile<-paramfile
idx1<-grep(pattern="^Zs",x=rawfile)
rawfile<-rawfile[idx1[1]]
rawfile<-str_split(rawfile,'%')[[1]][1]
Zs<-as.numeric(str_extract_all(rawfile,"[0-9]+")[[1]])
rm(idx1,rawfile)


#===============================================================================
# GET POI-INFORMATION
#===============================================================================
fn<-grep(pattern="^fnPOIs",simfile)
fn<-simfile[fn]
fn<-str_split(fn,"'")[[1]][2]

pois_df<-read.table(path_input&'/'&fn,header=T)
rm(fn)


#===============================================================================
# CREATE NEW FOLDERS
#===============================================================================
new_dir<-path_store&'/'&outnm
dir.create(new_dir)
dir.create(new_dir&'/TABLES')


#===============================================================================
# READ MODEL OUTPUT FILE
#===============================================================================
setwd(path_data&'/OUTPUTS')
fn_ls<-list.files(pattern='_results.mat')
n<-length(fn_ls)

for(i in 1:n){
  print('+++ processing <<'&fn_ls[i]&'>>   ('&i&'/'&n&')')
  matfile<-readMat(fn_ls[i])  ###expensive!
  nms_v<-names(matfile)
  
  if(i == 1){
    Date<-matfile[['Datam']]
    Date<-apply(Date,1,paste,collapse="-")
    Date<-as.POSIXct(Date,format='%Y-%m-%d-%H')
  }
  
  idx<-which(nms_v %in% var_df$var)
  NMS<-nms_v[idx]
  rm(idx)
  
  ### read hourly & daily timeseries
  df_h<-data.frame(Date=Date)
  df_d<-data.frame(Date=unique(cut(df_h$Date,'day'))) #derive daily mean
  df_d$Date<-as.POSIXct(df_d$Date,format='%Y-%m-%d')
  
  for(j in 1:length(NMS)){
    NM<-NMS[j]
    idx<-which(var_df$var == NM)
    dat<-data.frame(matfile[[NM]])
    
    if(var_df$per_soil_layer[idx] == FALSE){colnames(dat)<-NM}
    if(var_df$per_soil_layer[idx] == TRUE){colnames(dat)<-NM&'_ly'&1:ncol(dat)}
    if(NM == 'B.H' || NM == 'B.L'){colnames(dat)<-NM&'_'&c('Fol','LivSap','FinRoo','CHR','FrFl','Hea','SDFol','Aux')}
    
    if(var_df$temporal_resolution[idx] == 'hourly')
    {df_h<-cbind.data.frame(df_h,dat)}
    if(var_df$temporal_resolution[idx] == 'daily')
    {df_d<-cbind.data.frame(df_d,dat)}
    
    rm(dat,idx,NM)
  }
  
  ### AVERAGE hourly to daily values
  df<-aggregate(df_h[,-1],list(Date=cut(df_h$Date,'day')),mean)[,-1]
  #################################
  
  ### subset soil layers to actual soil depth
  idx<-fn_ls[i]
  ID<-str_split(idx,'-')[[1]][1]
  if(startsWith(idx,'CH-')){ID<-str_split(idx,'-FLU_')[[1]][1]}
  idx<-which(pois_df$ID == ID)
  SOILTH<-pois_df$SOILTH[idx] #[m]
  if(SOILTH == 0){
    rmidx<-which(grepl('_ly',names(df)))
    df<-df[,-rmidx]
    rm(rmidx)
  }
  if(SOILTH > 0){
    ly<-which.min(abs(Zs-(SOILTH*1000))) #layer no. that corresponds to soil thickness (from top to bottom)
    if(ly < length(Zs)){
      rmly<-(ly+1):(length(Zs)-1)                   #layers to remove
      rmidx<-which(endsWith(names(df),'_ly'&rmly))
      df<-df[,-rmidx]                               #remove layers below soil thickness
      rm(rmly,rmidx)
    }
    rm(ly)
  }
  
  ##prepare data
  df<-cbind.data.frame(df_d[,-1],df)
  df<-df[,order(names(df))]
  df<-cbind.data.frame(Date=df_d$Date,df)
  df[,-1]<-round(df[,-1],5)
  df$Date<-as.character(df$Date)
  
  ##write to file
  fn<-pois_df$ID[idx]&'-'&pois_df$ID2[idx]&'_'&pois_df$elevation[idx]&'m'
  fn<-new_dir&'/TABLES/'&fn
  write.table(df,fn&'___daily_avg.txt',dec='.',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
  rm(matfile,nms_v,NMS,idx,ID,fn,df,df_h,df_d)
}









