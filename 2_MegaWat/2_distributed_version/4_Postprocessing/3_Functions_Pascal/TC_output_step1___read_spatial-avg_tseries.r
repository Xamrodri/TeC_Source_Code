################################################################################
# T&C OUTPUT:
#  READ TIMESERIES OF SPATAL AVERAGES & WRITE TO TEXT FILES
#
# 
# TODO: -
# 
# NEW:  
#       
#
#
# 2023/09/30
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
new_dir<-path_store&'/'&outnm
dir.create(new_dir)
dir.create(new_dir&'/TABLES')


#===============================================================================
# DETECT IF REINITIATED
#===============================================================================
setwd(path_data&'/OUTPUTS')
files<-list.files()
files<-files[startsWith(files,'reinit_')]
if(length(files) > 1){reinit=TRUE}
if(length(files) < 1){reinit=FALSE}
rm(files)


#===============================================================================
# READ PARAMETER FILE << OUTPUT_MANAGER_PAR.m >>
#===============================================================================
setwd(path_data)
paramfile<-readLines('OUTPUT_MANAGER_PAR.m')
## if you want to match metacharacters (e.g. parentheses etc.) 
## you have to "scape" the metacharacter with two backslashes


#===============================================================================
# READ SIMULATION FILE << AAA_TC_XXX.m >> OR << BBB_TC_XXX_reinit.m >>
#  & GET SIMULATION NAMES
#===============================================================================
files<-list.files()
##read AAA-file
simfile<-files[startsWith(files,'AAA_')]
simfile<-readLines(simfile,warn=FALSE)

VERSION<-grep(pattern="^VERSION",simfile)
VERSION<-simfile[VERSION]
VERSION<-str_split(VERSION,"'")[[1]][2]

DOMAIN<-grep(pattern="^DOMAIN",simfile)
DOMAIN<-simfile[DOMAIN]
DOMAIN<-str_split(DOMAIN,"'")[[1]][2]

simnm<-VERSION&'_'&DOMAIN
simnm1<-simnm

##check for BBB-file (re-initialization)
if(reinit){
  simnm2<-files[startsWith(files,'BBB_')]
  simfile2<-readLines(simnm2,warn=FALSE)
  ##get secondary simulation name
  ## VERSION
  VERSION<-grep(pattern="^VERSION",x=simfile2)
  VERSION<-simfile2[VERSION[1]]
  VERSION<-str_match(VERSION,"'\\s*(.*?)\\s*'")[,2]
  ## DOMAIN
  DOMAIN<-grep(pattern="^DOMAIN",x=simfile2)
  DOMAIN<-simfile2[DOMAIN[1]]
  DOMAIN<-str_match(DOMAIN,"'\\s*(.*?)\\s*'")[,2]
  simnm2<-VERSION&'_'&DOMAIN
  simnm<-simnm2
}
rm(files)


#===============================================================================
# GET MODEL PERIOD DATES
#===============================================================================
# initial timestamp
ts<-grep(pattern="^mod_ini_Date",simfile)
ts<-simfile[ts]
ts<-str_split(ts,"'")[[1]][2]
date_ini<-as.POSIXct(ts,format='%d-%b-%Y %H:%M:%S')
rm(ts)
# final timestamp
ts<-grep(pattern="^mod_fin_Date",simfile)
ts<-simfile[ts]
ts<-str_split(ts,"'")[[1]][2]
date_fin<-as.POSIXct(ts,format='%d-%b-%Y %H:%M:%S')
rm(ts)
Date<-seq(date_ini,date_fin,by=3600)


#===============================================================================
# READ LAND COVER CLASSES
#===============================================================================
rawfile<-paramfile
# find lines where classes are defined
idx1<-grep(pattern="START LAND COVER CLASSES",x=rawfile)
idx2<-grep(pattern="END LAND COVER CLASSES",x=rawfile)
if(length(idx1) == 0 | length(idx2) == 0){stop('Lines in "paramfile" not found!')} 
rawfile<-rawfile[idx1:idx2]
# find all variable-lines
vars<-grepl(pattern='[0-9]',x=rawfile)
rawfile<-rawfile[vars]
vars<-lapply(strsplit(rawfile,'%%%'), function(x) x[2])
vars<-do.call('rbind',vars)
vars<-lapply(strsplit(vars,' = '), function(x) c(x[1],x[2]))
vars<-do.call('rbind.data.frame',vars)
colnames(vars)<-c('Code','Class')
vars$Code<-as.numeric(as.character(vars$Code))
vars$Class<-as.character(vars$Class)
LandcoverClasses_df<-vars
##add type
LandcoverClasses_df$Type<-NA
LandcoverClasses_df[which(LandcoverClasses_df$Class %in% c('Evergreen','Mixed-Evergreen-Decidous','Decidous','Grass')),'Type']<-'vegetated'
LandcoverClasses_df[which(LandcoverClasses_df$Class %in% c('Rock')),'Type']<-'rocky'
LandcoverClasses_df[which(LandcoverClasses_df$Class %in% c('Water')),'Type']<-'water'
LandcoverClasses_df[which(LandcoverClasses_df$Class %in% c('Ice','Clean-ice','Debris-covered-ice')),'Type']<-'glacierized'
rm(idx1,idx2,rawfile,vars)


#===============================================================================
# READ HEADERS FOR "SPATIAL AVERAGE OVER THE WATERSHED" 
#===============================================================================
rawfile<-paramfile
# find lines where output variables are defined
idx1<-grep(pattern="START <<OUTPUT_ZZZ_AVG.dat>>",x=rawfile)
idx2<-grep(pattern="END <<OUTPUT_ZZZ_AVG.dat>>",x=rawfile)
if(length(idx1) == 0 | length(idx2) == 0){stop('Lines in "paramfile" not found!')} 
rawfile<-rawfile[idx1:idx2]
# find all variable-lines
vars<-grepl(pattern="fprintf",x=rawfile)
rawfile<-rawfile[vars]
vars<-lapply(strsplit(rawfile,','), function(x) x[3])
vars<-do.call('rbind',vars)
vars<-lapply(strsplit(vars,');'), function(x) x[1])
vars<-do.call('c',vars)
vars<-str_remove(vars,'_tg')
headers_SPAVG<-vars
rm(idx1,idx2,rawfile,vars)


#===============================================================================
# READ HEADERS FOR "SPATIAL AVERAGE OVER THE VEGETATION TYPE"
#===============================================================================
rawfile<-paramfile
# find lines where output variables are defined
idx1<-grep(pattern="START <<OUTPUT_ZZZ_AVG_PFT_YYY_code_XXX.dat>>",x=rawfile)
idx2<-grep(pattern="END <<OUTPUT_ZZZ_AVG_PFT_YYY_code_XXX.dat>>",x=rawfile)
if(length(idx1) == 0 | length(idx2) == 0){stop('Lines in "paramfile" not found!')}
rawfile<-rawfile[idx1:idx2]
# find all variable-lines
vars<-grepl(pattern="fprintf",x=rawfile)
rawfile<-rawfile[vars]
vars<-lapply(strsplit(rawfile,','), function(x) x[4])
vars<-do.call('rbind',vars)
vars<-lapply(strsplit(vars,'\\('), function(x) x[1])
vars<-do.call('c',vars)
vars<-str_remove(vars,'_tg')
headers_SPAVG_PFT<-vars
rm(idx1,idx2,rawfile,vars)


#===============================================================================
# READ HEADERS FOR "SPATIAL AVERAGE OVER EACH LAND COVER CLASS"
#===============================================================================
rawfile<-paramfile
# find lines where output variables are defined
idx1<-grep(pattern="START <<OUTPUT_ZZZ_AVG_LC_code_YYY.dat>>",x=rawfile)
idx2<-grep(pattern="END <<OUTPUT_ZZZ_AVG_LC_code_YYY.dat>>",x=rawfile)
if(length(idx1) == 0 | length(idx2) == 0){stop('Lines in "paramfile" not found!')}
rawfile<-rawfile[idx1:idx2]
# find all variable-lines
vars<-grepl(pattern="fprintf",x=rawfile)
rawfile<-rawfile[vars]
vars<-lapply(strsplit(rawfile,','), function(x) x[3])
vars<-do.call('rbind',vars)
vars<-lapply(strsplit(vars,'\\('), function(x) x[1])
vars<-do.call('c',vars)
vars<-str_remove(vars,'_tg_LC')
headers_SPAVG_LC<-vars
rm(idx1,idx2,rawfile,vars)


#===============================================================================
# ADJUST DATES-VECTOR IN CASE MODEL WAS NOT REACHING LAST PREDEFINED TIMESTEP
#===============================================================================
setwd(path_data&'/OUTPUTS')
fn<-'OUTPUT_'&simnm&'_AVG.dat'
df<-read.table(fn,header=FALSE,sep='\t',dec='.')
colnames(df)<-headers_SPAVG
rm(fn)
if(reinit){
  fn2<-'reinit___OUTPUT_'&simnm&'_AVG.dat'
  df2<-read.table(fn2,header=FALSE,sep='\t',dec='.')
  colnames(df2)<-headers_SPAVG
  n<-nrow(df); n2<-nrow(df2)
  ##check for overlap after restart
  idx1<-which(outer(df[(n-24):n,'Ws'],df2[1:25,'Ws'],'=='),arr.ind=T)
  idx2<-which(outer(df[(n-24):n,'Pre'],df2[1:25,'Pre'],'=='),arr.ind=T)
  if(all(idx1[,2] == idx2[,2])){
    double_ts<-sort(n+1-idx1[,2])
  }
  rm(n,idx1,idx2)
  
  df<-rbind.data.frame(df,df2)
  if(exists('double_ts')){df<-df[-double_ts,]}
  rm(fn2,df2)
}
nts<-nrow(df)
rm(df)
if(length(Date) > nts){
  Date<-Date[1:nts]
}


#===============================================================================
# READ OUTPUT OF "SPATIAL AVERAGE OVER THE WATERSHED"
#===============================================================================
fn<-path_data&'/OUTPUTS/OUTPUT_'&simnm&'_AVG.dat'
spavg_df<-read.table(fn,header=FALSE,sep='\t',dec='.')
rm(fn)
if(reinit){
  fn2<-path_data&'/OUTPUTS/reinit___OUTPUT_'&simnm&'_AVG.dat'
  spavg_df2<-read.table(fn2,header=FALSE,sep='\t',dec='.')
  
  spavg_df<-rbind.data.frame(spavg_df,spavg_df2)
  if(exists('double_ts')){spavg_df<-spavg_df[-double_ts,]}
  rm(fn2,spavg_df2)
}
colnames(spavg_df)<-headers_SPAVG
spavg_df<-as.data.frame(data.matrix(spavg_df,rownames.force=NA)) #all cols to numeric
if(nts < nrow(spavg_df)){spavg_df<-spavg_df[-nts+1,]} ##hack!
spavg_df<-cbind.data.frame(Date=Date,spavg_df)
##write to file
write.table(spavg_df,new_dir&'/TABLES/SPAVG.txt',dec='.',sep='\t',quote=FALSE,
            col.names=TRUE,row.names=FALSE)



#===============================================================================
# READ OUTPUT OF "SPATIAL AVERAGE OVER THE VEGETATION TYPE"
#===============================================================================
VegetationTypes_df<-LandcoverClasses_df[which(LandcoverClasses_df$Type == 'vegetated'),]
n<-nrow(VegetationTypes_df)
for(i in 1:n){
  code<-VegetationTypes_df$Code[i]
  
  print(DOMAIN&': reading spatial average over PFT '&code&'  ('&i&'/'&n&')')
  
  ###PFT 1
  fn<-'OUTPUT_'&simnm&'_AVG_PFT_1_code_'&code&'.dat'
  if(file.exists(fn)){
    df<-read.table(fn,header=FALSE,sep='\t',dec='.')
    
    if(reinit){
      fn2<-'reinit___OUTPUT_'&simnm&'_AVG_PFT_1_code_'&code&'.dat'
      df2<-read.table(fn2,header=FALSE,sep='\t',dec='.')
      
      df<-rbind.data.frame(df,df2)
      if(exists('double_ts')){df<-df[-double_ts,]}
      rm(fn2,df2)
    }
    colnames(df)<-headers_SPAVG_PFT
    df<-as.data.frame(data.matrix(df,rownames.force=NA)) #all cols to numeric
    if(nts < nrow(df)){df<-df[-nts+1,]} ##hack!
    ifelse(nts > nrow(df),
           df<-cbind.data.frame(Date=Date[1:nrow(df)],df), ##hack!
           df<-cbind.data.frame(Date=Date,df[1:nts,])
    )
    ##write to file
    write.table(df,new_dir&'/TABLES/SPAVG-PFT1___'&
                  VegetationTypes_df$Class[i]&'.txt',dec='.',
                sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
    rm(df,fn)
  }
  ###PFT 2
  fn<-'OUTPUT_'&simnm&'_AVG_PFT_2_code_'&code&'.dat'
  if(file.exists(fn)){
    df<-read.table(fn,header=FALSE,sep='\t',dec='.')
    
    if(reinit){
      fn2<-'reinit___OUTPUT_'&simnm&'_AVG_PFT_2_code_'&code&'.dat'
      df2<-read.table(fn2,header=FALSE,sep='\t',dec='.')
      
      df<-rbind.data.frame(df,df2)
      if(exists('double_ts')){df<-df[-double_ts,]}
      rm(fn2,df2)
    }
    colnames(df)<-headers_SPAVG_PFT
    df<-as.data.frame(data.matrix(df,rownames.force=NA)) #all cols to numeric
    if(nts < nrow(df)){df<-df[-nts+1,]} ##hack!
    ifelse(nts > nrow(df),
           df<-cbind.data.frame(Date=Date[1:nrow(df)],df), ##hack!
           df<-cbind.data.frame(Date=Date,df[1:nts,])
    )
    ##write to file
    write.table(df,new_dir&'/TABLES/SPAVG-PFT2___'&
                  VegetationTypes_df$Class[i]&'.txt',dec='.',
                sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
    rm(df,fn)
  }
  rm(code)
}
rm(n)


#===============================================================================
# READ OUTPUT OF "SPATIAL AVERAGE OVER EACH LAND COVER CLASS"
#===============================================================================
n<-nrow(LandcoverClasses_df)
for(i in 1:n){
  code<-LandcoverClasses_df$Code[i]
  
  print(DOMAIN&': reading spatial average over LC class '&code&'  ('&i&'/'&n&')')
  
  fn<-'OUTPUT_'&simnm&'_AVG_LC_code_'&code&'.dat'
  if(file.exists(fn)){
    df<-read.table(fn,header=FALSE,sep='\t',dec='.')
    
    if(reinit){
      fn2<-'reinit___OUTPUT_'&simnm&'_AVG_LC_code_'&code&'.dat'
      df2<-read.table(fn2,header=FALSE,sep='\t',dec='.')
      
      df<-rbind.data.frame(df,df2)
      if(exists('double_ts')){df<-df[-double_ts,]}
      rm(fn2,df2)
    }
    colnames(df)<-headers_SPAVG_LC
    df<-as.data.frame(data.matrix(df,rownames.force=NA)) #all cols to numeric
    if(nts < nrow(df)){df<-df[-nts+1,]} ##hack!
    ifelse(nts > nrow(df),
           df<-cbind.data.frame(Date=Date[1:nrow(df)],df), ##hack!
           df<-cbind.data.frame(Date=Date,df[1:nts,])
    )
    ##write to file
    write.table(df,new_dir&'/TABLES/SPAVG-LC___'&
                  LandcoverClasses_df$Class[i]&'.txt',dec='.',
                sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
    rm(df,fn)
  }
  rm(code)
}
rm(n)

###
print('+++END SCRIPT <<READ SPATIAL-AVG TSERIES>>+++   '&outnm)
