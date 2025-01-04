################################################################################
# T&C OUTPUT:
#  READ TIMESERIES OF TRACKED PIXELS & WRITE TO TEXT FILES
#
# 
# TODO: -
# 
# NEW:  
#       
#
#
# 2023/11/13
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
library(data.table)   #%like%

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
dir.create(new_dir&'/SHAPEFILES')
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
date_ini<-as.POSIXct(ts,format='%d-%b-%Y %H:%M:%S',tz=TZ)
rm(ts)
# final timestamp
ts<-grep(pattern="^mod_fin_Date",simfile)
ts<-simfile[ts]
ts<-str_split(ts,"'")[[1]][2]
date_fin<-as.POSIXct(ts,format='%d-%b-%Y %H:%M:%S',tz=TZ)
rm(ts)
Date<-seq(date_ini,date_fin,by=3600)


#===============================================================================
# GET INPUT FOLDER NAME
#===============================================================================
idx<-grep(pattern="^GEODATAINPUT = \\[",simfile)
fo<-simfile[idx]
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
# READ INITIAL CONDITIONS RASTERS
#===============================================================================
### DEBRIS THICKNESS
debth<-geodatafile[['DEB.MAP']]
debth[debth < 0]<-NA
debth<-apply(debth,2,rev) 
debth<-raster(nrow=m_cell,ncol=n_cell,vals=debth)
projection(debth)<-as.character(projec)
extent(debth)<-ext

### DEM
dem<-geodatafile[['DTM']]
dem[dem < 0]<-NA
dem<-apply(dem,2,rev)   ## Matlab: flipud(DTM)
dem<-raster(nrow=m_cell,ncol=n_cell,vals=dem)
projection(dem)<-as.character(projec)
extent(dem)<-ext

### GLACIER THICKNESS
glath<-geodatafile[['GLH']]
glath[glath < 0]<-NA
glath<-apply(glath,2,rev)
glath<-raster(nrow=m_cell,ncol=n_cell,vals=glath)
projection(glath)<-as.character(projec)
extent(glath)<-ext

### PCLA MAP
pcla<-geodatafile[['PCLA']]
pcla[pcla < 0]<-NA
pcla<-apply(pcla,2,rev) 
pcla<-raster(nrow=m_cell,ncol=n_cell,vals=pcla)
projection(pcla)<-as.character(projec)
extent(pcla)<-ext

### PORG MAP
porg<-geodatafile[['PORG']]
porg[porg < 0]<-NA
porg<-apply(porg,2,rev) 
porg<-raster(nrow=m_cell,ncol=n_cell,vals=porg)
projection(porg)<-as.character(projec)
extent(porg)<-ext

### PSAN MAP
psan<-geodatafile[['PSAN']]
psan[psan < 0]<-NA
psan<-apply(psan,2,rev)
psan<-raster(nrow=m_cell,ncol=n_cell,vals=psan)
projection(psan)<-as.character(projec)
extent(psan)<-ext

### SOIL THICKNESS
soilth<-geodatafile[['SOIL.TH']]
soilth[soilth < 0]<-NA
soilth<-apply(soilth,2,rev) 
soilth<-raster(nrow=m_cell,ncol=n_cell,vals=soilth)
projection(soilth)<-as.character(projec)
extent(soilth)<-ext

### SNOW DEPTH
snod<-geodatafile[['SNOWD']]
snod[snod < 0]<-NA
snod<-apply(snod,2,rev)  
snod<-raster(nrow=m_cell,ncol=n_cell,vals=snod)
projection(snod)<-as.character(projec)
extent(snod)<-ext

### LAND COVER CODE
lcc<-geodatafile[['VEG.CODE']]
lcc[lcc < 0]<-NA
lcc<-apply(lcc,2,rev) 
lcc<-raster(nrow=m_cell,ncol=n_cell,vals=lcc)
projection(lcc)<-as.character(projec)
extent(lcc)<-ext

### CALCULATE SLOPE & ASPECT
slo<-terrain(dem,'slope','degrees')
asp<-terrain(dem,'aspect','degrees')
extent(slo)<-extent(dem)
extent(asp)<-extent(dem)


#===============================================================================
# READ POIs COORDINATES & CONVERT TO SPATIAL POINTS
#===============================================================================
setwd(inputfo)
files<-list.files()
pois<-read.csv(files[which(str_starts(files,'POI_'))],header=TRUE, sep=',')
pois<-st_as_sf(pois,coords=c("chx","chy"),remove=FALSE)
st_crs(pois)<-projec
pois<-as(pois,'Spatial')


#===============================================================================
# EXTRACT POINT-SPECIFIC MODEL INPUT VARIABLES
#===============================================================================
pois$ASP<-raster::extract(asp,pois)
pois$DEBTH<-raster::extract(debth,pois)
pois$ELE_DEM<-raster::extract(dem,pois)
pois$GLATH<-raster::extract(glath,pois)
pois$LCC<-raster::extract(lcc,pois)
pois$PCLA<-raster::extract(pcla,pois)
pois$PORG<-raster::extract(porg,pois)
pois$PSAN<-raster::extract(psan,pois)
pois$SLO<-raster::extract(slo,pois)
pois$SNOD<-raster::extract(snod,pois)
pois$SOILTH<-raster::extract(soilth,pois)


#===============================================================================
# WRITE POI-DATA TO .TXT-FILE
#===============================================================================
df<-as.data.frame(pois)[,1:ncol(pois)]

## 1 digit
idx<-which(names(df) %in% c('ASP','PCLA','PORG','PSAN','SLO'))
df[,idx]<-round(df[,idx],1)
## 3 digits
idx<-which(names(df) %in% c('SNOD','SOILTH'))
df[,idx]<-round(df[,idx],3)
rm(idx)

## sort by name
poi_df<-df[order(df$name),]
rm(df)


################################
##check & correct pixel index
setwd(path_data&'/OUTPUTS')
idxs<-list.files(pattern='_PIXEL_')
idxs<-str_split(idxs,'_',simplify=T)
idxs<-idxs[idxs %like% '.dat']
idxs<-str_split(idxs,'.dat',simplify=T)
idxs<-sort(unique(as.numeric(idxs[,1])))

for(i in 1:nrow(poi_df)){
  if(file_test("-f",'OUTPUT_'&simnm&'_PIXEL_'&poi_df$idx[i]&'.dat') == FALSE){
    idx<-which.min(abs(idxs-poi_df$idx[i])) #take ~next (index-wise) pixel & overwrite
    idx<-idxs[idx]
    poi_df$idx[i]<-idx
    rm(idx)
  }
}
rm(idxs)
################################

## write (updated) table
write.table(poi_df,new_dir&'/TABLES/POIs.txt',dec='.',sep='\t',quote=FALSE,
            col.names=TRUE,row.names=FALSE)


#===============================================================================
# WRITE POIs TO .SHP-FILE
#===============================================================================
st_write(st_as_sf(pois),new_dir&'/SHAPEFILES/POIs.shp',
         driver="ESRI Shapefile",delete_dsn=TRUE)


#===============================================================================
# READ HEADERS FOR "TRACKED PIXELS TIME SERIES" (NORMAL OUTPUT)
#===============================================================================
rawfile<-paramfile
# find lines where output variables are defined
idx1<-grep(pattern="START <<OUTPUT_ZZZ_PIXEL_YYY.dat>>",x=rawfile)
idx2<-grep(pattern="END <<OUTPUT_ZZZ_PIXEL_YYY.dat>>",x=rawfile)
if(length(idx1) == 0 | length(idx2) == 0){stop('Lines in "paramfile" not found!')}
rawfile<-rawfile[idx1:idx2]
# find all variable-lines
vars<-grepl(pattern="fprintf",x=rawfile)
rawfile<-rawfile[vars]
vars<-lapply(strsplit(rawfile,','), function(x) x[4])
vars<-do.call('rbind',vars)
vars<-lapply(strsplit(vars,'\\(i'), function(x) x[1])
vars<-do.call('c',vars)
headers_PIXEL<-vars
rm(idx1,idx2,rawfile,vars)


#===============================================================================
# READ HEADERS FOR "TRACKED PIXELS TIME SERIES" (SOIL OUTPUT)
#===============================================================================
rawfile<-paramfile
# find lines where output variables are defined
idx1<-grep(pattern="START <<OUTPUT_ZZZ_PIXEL_SOIL_YYY.dat>>",x=rawfile)
idx2<-grep(pattern="END <<OUTPUT_ZZZ_PIXEL_SOIL_YYY.dat>>",x=rawfile)
if(length(idx1) == 0 | length(idx2) == 0){stop('Lines in "paramfile" not found!')}
rawfile<-rawfile[idx1:idx2]
# find all variable-lines
vars<-grepl(pattern="fprintf",x=rawfile)
rawfile<-rawfile[vars]
vars<-lapply(strsplit(rawfile,','), function(x) x[3])
vars<-do.call('rbind',vars)
vars<-lapply(strsplit(vars,'\\(i'), function(x) x[1])
vars<-do.call('c',vars)
headers_PIXEL_SOIL<-vars
rm(idx1,idx2,rawfile,vars)

# find line where ms_max is defined
rawfile<-simfile
idx1<-grep(pattern="ms_max",x=rawfile)
rawfile<-rawfile[idx1[1]]
ms_max<-as.numeric(str_extract_all(rawfile,"[0-9]+")[[1]])
last<-length(headers_PIXEL_SOIL) ##exeption "CK1"
rm(idx1,rawfile)
headers_PIXEL_SOIL<-c(rep(headers_PIXEL_SOIL[-last],ms_max),headers_PIXEL_SOIL[last])
headers_PIXEL_SOIL<-headers_PIXEL_SOIL&'_lr'&
  rep(1:ms_max,each=length(headers_PIXEL_SOIL)/ms_max)


#===============================================================================
### TO DELETE IN NEWER VERSIONS !!!!
#    !!!!
# FIND DUPLICATED TRACK PIXELS
# (e.g. when station & outlet are at the same pixel)
#===============================================================================
## idx<-which(duplicated(poi_df$idx))
## idx<-which(poi_df$idx == poi_df$idx[idx])
## which(poi_df$ID2[idx] == 'HYD')

# fn<-path_data&'/OUTPUTS/OUTPUT_'&simnm&'_PIXEL_'&poi_df$idx&'.dat'
# if(any(duplicated(fn))){
#   dupidx<-which(duplicated(fn))
#   fn<-fn[-dupidx]
#   rm(dupidx)
# }
# if(reinit){
  # fn2<-path_data&'/OUTPUTS/reinit___OUTPUT_'&simnm&'_PIXEL_'&poi_df$idx&'.dat'
  # if(any(duplicated(fn2))){
  #   dupidx<-which(duplicated(fn2))
  #   fn2<-fn2[-dupidx]
  #   rm(dupidx)
  # }
# }


fn<-path_data&'/OUTPUTS/OUTPUT_'&simnm&'_PIXEL_'&poi_df$idx&'.dat'
if(reinit){
  fn2<-path_data&'/OUTPUTS/reinit___OUTPUT_'&simnm&'_PIXEL_'&poi_df$idx&'.dat'
}


#===============================================================================
# ADJUST DATES-VECTOR IN CASE MODEL WAS NOT REACHING LAST PREDEFINED TIMESTEP
#===============================================================================
df<-read.table(path_data&'/OUTPUTS/OUTPUT_'&simnm&'_PIXEL_'&poi_df[1,'idx']&'.dat',
               header=FALSE,sep='\t',dec='.')
colnames(df)<-headers_PIXEL
if(reinit){
  df2<-read.table(path_data&'/OUTPUTS/reinit___OUTPUT_'&simnm&'_PIXEL_'&poi_df[1,'idx']&'.dat',
                  header=FALSE,sep='\t',dec='.')
  colnames(df2)<-headers_PIXEL
  n<-nrow(df); n2<-nrow(df2)
  ##check for overlap after restart
  idx1<-which(outer(df[(n-24):n,'Ws_S'],df2[1:25,'Ws_S'],'=='),arr.ind=T)
  idx2<-which(outer(df[(n-24):n,'Pre'],df2[1:25,'Pre'],'=='),arr.ind=T)
  if(all(idx1[,2] == idx2[,2])){
    double_ts<-sort(n+1-idx1[,2])
  }
  rm(n,idx1,idx2)
  
  df<-rbind.data.frame(df,df2)
  if(exists('double_ts')){df<-df[-double_ts,]}
  rm(df2)
}

# ###hack, to delete !!
# if(exists('dupidx')){
#   if(sum(df[1,]-df[2,],na.rm=T)+
#      sum(df[3,]-df[4,],na.rm=T)+
#      sum(df[5,]-df[6,],na.rm=T) == 0){
#     df<-df[-(seq(2,to=nrow(df),by=2)),]
#   }
# }
# ######

nts<-nrow(df)
rm(df)
if(length(Date) > nts){
  Date<-Date[1:nts]
}


#===============================================================================
# READ OUTPUT OF "TRACKED PIXELS TIME SERIES" (NORMAL OUTPUT)
#===============================================================================
n<-length(fn)
for(i in 1:n){
  print(DOMAIN&': reading POI <<'&poi_df$name[i]&' ('&poi_df$ID2[i]&') >>  ('&i&'/'&n&')')
  
  df<-read.table(fn[i],header=FALSE,sep='\t',dec='.')
  
  if(!is.na(df[2,10])){ ###hack!!
    if(df[1,10]-df[2,10] == 0 && df[1,13]-df[2,13] == 0 && df[1,24]-df[2,24] == 0){ ##hack!!
      df<-Nth.delete(df,2)
    }
    
    if(reinit){
      df2<-read.table(fn2[i],header=FALSE,sep='\t',dec='.')
      df<-rbind.data.frame(df,df2)
      if(exists('double_ts')){df<-df[-double_ts,]}
      rm(df2)
    }
    colnames(df)<-headers_PIXEL
    df<-as.data.frame(data.matrix(df,rownames.force=NA)) #all cols to numeric
    df<-df[,-which(is.na(colnames(df)))]
    # if(exists('dupidx')){     ###hack, to delete !!
    #   if(i %in% dupidx){
    #     df<-df[seq(1,nrow(df),by=2),]
    #   }
    # }
    if(nts < nrow(df)){df<-df[-nts+1,]}   ##hack!!
    if(nts < nrow(df)){df<-df[-nts+2,]}   ##hack!!
    if(nts < nrow(df)){df<-df[-nts+3,]}   ##hack!!
    if(nts < nrow(df)){df<-df[1:nts,]}   ##hack!!
    Date2<-Date
    if(nts > nrow(df)){Date2<-Date[1:nrow(df)]}   ##hack!!
    df<-cbind.data.frame(Date=Date2,df)
    ##write to file
    write.table(df,new_dir&'/TABLES/TrackedPixel___'&poi_df[i,'ID']&'-'&
                  poi_df[i,'ID2']&'___'&poi_df[i,'elevation']&'m.txt',dec='.',
                sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
  }
  rm(df)
}

# (QpointC/dth)*(cellsize^2)/(3600000); %% [m^3/s]


#===============================================================================
# READ OUTPUT OF "TRACKED PIXELS TIME SERIES" (SOIL OUTPUT)
#===============================================================================
for(i in 1:length(fn)){
  print(DOMAIN&': reading POI (soil output) <<'&poi_df$name[i]&' ('&poi_df$ID2[i]&') >>  ('&i&'/'&n&')')
  
  df<-read.table(fn[i],header=FALSE,sep='\t',dec='.')
  
  if(!is.na(df[2,10])){ ###hack!!
    if(df[1,10]-df[2,10] == 0 && df[1,13]-df[2,13] == 0 && df[1,24]-df[2,24] == 0){ ##hack!!
      df<-Nth.delete(df,2)
    }
    
    if(reinit){
      df2<-read.table(fn2[i],header=FALSE,sep='\t',dec='.')
      df<-rbind.data.frame(df,df2)
      if(exists('double_ts')){df<-df[-double_ts,]}
      rm(df2)
    }
    colnames(df)<-headers_PIXEL_SOIL
    df<-as.data.frame(data.matrix(df,rownames.force=NA)) #all cols to numeric
    # if(exists('dupidx')){  ###hack, to delete !!
    #   if(i %in% dupidx){
    #     df<-df[seq(1,nrow(df),by=2),]
    #   }
    # }
    if(nts < nrow(df)){df<-df[-nts+1,]}   ##hack!!
    if(nts < nrow(df)){df<-df[-nts+2,]}   ##hack!!
    if(nts < nrow(df)){df<-df[-nts+3,]}   ##hack!!
    if(nts < nrow(df)){df<-df[1:nts,]}   ##hack!!
    Date2<-Date
    if(nts > nrow(df)){Date2<-Date[1:nrow(df)]}   ##hack!!
    df<-cbind.data.frame(Date=Date2,df)
    ##write to file
    write.table(df,new_dir&'/TABLES/TrackedPixel-Soil___'&poi_df[i,'ID']&'-'&
                  poi_df[i,'ID2']&'___'&poi_df[i,'elevation']&'m.txt',dec='.',
                sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
  }
  rm(df)
}



###
print('+++END SCRIPT <<READ TRACKED PIXELS>>+++   '&outnm)

