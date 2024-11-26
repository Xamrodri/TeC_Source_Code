################################################################################
# T&C OUTPUT STEP 2:
#  PLOT ANNUAL AVG/SUMS OVER ELEVATION (FROM MONTHLY RASTERS)
#
# 
# TODO: -
#
#
# 2023/11/23
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

## define timespan & months to analyze
ts_init<-'2016-01-01 02:00:00'
# ts_final<-'2020-11-30 23:00:00'      #v20231006_Engadin_1906days
# ts_final<-'2019-12-31 23:00:00'      #v20231006_Jura_1562days
# ts_final<-'2020-09-30 23:00:00'      #v20231006_Hinterrhein_1853days
# ts_final<-'2019-10-31 23:00:00'      #v20231006_Oberwallis_1498days
# ts_final<-'2018-12-31 23:00:00'      #v20231006_Reuss_1257days
ts_final<-'2022-09-30 23:00:00'      #v20231006_Schaffhausen_2557days



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
library(scales)       #date_format()
library(networkD3)    #Sankey diagram
library(dplyr)        #Sankey diagram
library(webshot)      #save HTML screenshot
library(rasterVis)    #levelplot
library(gridExtra)    #grid.arrange()
library(metR)         #geom_text_contour()
library(RColorBrewer) #brewer.pal()
library(lubridate)    #floor_date(), hour(), %m+% 
library(extrafont)
library(grid)
library(gridExtra)
library(gridtext)
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
# ADD VARIABLE NAMES & CORRESPONDING UNITS
#===============================================================================
setwd(path_code)
setwd('..')
var_df<-read.table('TC_Variables.txt',header=T,sep='\t')
var_df[]<-lapply(var_df,as.character)


#===============================================================================
# READ DEM
#===============================================================================
## DEM
dem<-raster(path_data&'/RASTERS/INIT-COND/DTM.tif')
resol<-res(dem)[1]
slo<-raster::terrain(dem,opt='slope',unit='degrees')
asp<-raster::terrain(dem,opt='aspect',unit='degrees')


#===============================================================================
# RETRIEVE DATES FROM SIMULATION
#===============================================================================
#read catchment tables
setwd(path_data&'/TABLES')
spavg_df<-read.table('SPAVG.txt',header=T,sep='\t')
spavg_df$Date<-as.POSIXct(spavg_df$Date,format='%Y-%m-%d %H:%M:%S')
# 
# #read LC tables
# LC_A_df<-read.table('LC_Areas.txt',header=T,sep='\t')
# spavg_LC_ls<-list()
# for(i in 1:nrow(LC_A_df)){
#   nm<-as.character(LC_A_df$Class[i])
#   df<-read.table('SPAVG-LC___'&nm&'.txt',header=T,sep='\t')
#   df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S')
#   spavg_LC_ls[[i]]<-df
#   rm(nm,df)
# }
# names(spavg_LC_ls)<-LC_A_df$Class
# 
# #read Vegetation Type tables
# VegType_A_df<-read.table('VegType_Areas.txt',header=T,sep='\t')
# VegType_A_df$Type<-as.character(VegType_A_df$Type)
# VegType_A_df<-VegType_A_df[-(which(startsWith(VegType_A_df$Type,'Rock'))),]
# spavg_VegType_ls<-list()
# for(i in 1:nrow(VegType_A_df)){
#   nm<-as.character(VegType_A_df$Type[i])
#   df<-read.table('SPAVG-VegType___'&nm&'.txt',header=T,sep='\t')
#   df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S')
#   spavg_VegType_ls[[i]]<-df
#   rm(nm,df)
# }
# names(spavg_VegType_ls)<-VegType_A_df$Type


#===============================================================================
# TIME PERIOD OF INTEREST
#===============================================================================
ts_init<-as.POSIXct(ts_init,format='%Y-%m-%d %H:%M:%S')
ts_final<-as.POSIXct(ts_final,format='%Y-%m-%d %H:%M:%S')


PERIOD_SH<-format(ts_init,"%b%Y")&'-'&format(ts_final-86400,"%b%Y")
PERIOD_LO<-format(ts_init,"%d %b %Y")&' - '&format(ts_final-86400,"%d %b %Y")

t1<-spavg_df$Date[1]
tn<-spavg_df$Date[nrow(spavg_df)]
ini<-round(as.numeric(ts_init - t1)/30.42)+1
fin<-round(as.numeric(ts_final - t1)/30.42)
NMONTHS<-ini:fin
rm(t1,tn,ini,fin)

#index for selected period
idx_ts<-which(spavg_df$Date >= ts_init & spavg_df$Date <= ts_final)
dates_all<-spavg_df$Date
dates_period<-spavg_df$Date[idx_ts]
dates_period<-dates_period[1:(length(dates_period)-3)]


#===============================================================================
# DERIVE MONTHLY VALUES AVERAGED/WEIGHTED PER ELEVATION BAND
#===============================================================================
month_v<-unique(month(dates_period))
year_month_v<-unique(paste0(year(dates_period),'-',month(dates_period)))
idx<-c(which(diff(month(dates_period)) != 0),length(dates_period))
month_year_v<-unique(paste0(month.abb[month_v],' ',year(dates_period[idx])))
mon_v<-month.abb[month_v]
Dates<-seq(ts_init,ts_final,'1 month')


### LOOP: VARIABLES PER MONTH
eb<-200 #define elevation band vertical distance to aggreagte
df_peb_ls<-list()
df_cumpeb_ls<-list()
for(k in 1:length(year_month_v)){
  mon<-month_v[k]
  yearmon<-year_month_v[k]
  monyear<-month_year_v[k]
  
  print(DOMAIN&': processing rasters for << '&monyear&' >>')
  hrs<-which(paste0(year(dates_period),'-',month(dates_period)) == yearmon)
  nHours<-length(hrs) 
  nDays<-length(hrs)/24 
  rm(hrs)
  
  setwd(path_data&'/RASTERS/MONTHLY')
  fls<-list.files()
  
  ## 1: TOTAL EVAPOTRANSPIRATION 
  nms<-c('EG','EICE','EIn_H','EIn_L','EIn_rock','ESN','EWAT','T_H','T_L')
  r_ls<-list()
  for(i in 1:length(nms)){
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nHours    
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  asp_v<-extract(asp,pts[,1:2])
  slo_v<-extract(slo,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,asp=asp_v,slo=slo_v,vals=vals)
  df<-df[order(df$ele),]
  df1<-df                     
  catchmentA<-sum(as.vector(dem*0+1),na.rm=T)*(250^2)/10^6 #km2
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  
  ## 2: PRECIPITATION 
  nms<-'Pr'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]]
  r[is.na(r)]<-0   ###NAs from RhiresD
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df2<-df                 
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 3: ICE MELT
  nms<-'Imelt'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df3<-df                 
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 4: SNOW MELT
  nms<-'Smelt'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours     
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df4<-df                 
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 5: SNOW FALL
  nms<-'Pr_sno'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df5<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 6: RAIN
  nms<-'Pr_liq'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df6<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 7: SOIL WATER
  nms<-'V'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df7<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 8: FROZEN SOIL WATER
  nms<-'Vice'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df8<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 9: EVAPORATION FROM SNOW AND ICE
  nms<-c('EICE','ESN')
  r_ls<-list()
  for(i in 1:length(nms)){
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nHours         
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df9<-df
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  
  ## 10: TRANSPIRATION
  nms<-c('T_L','T_H')
  r_ls<-list()
  for(i in 1:length(nms)){   
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nHours      
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df10<-df
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  
  ## 11: GROUND EVAPORATION (SOIL, SNOW, WATER)
  nms<-c('EG','ESN','EWAT')
  r_ls<-list()
  for(i in 1:length(nms)){ 
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nHours       
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df11<-df
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  
  ## 12: EVAPORATION FROM INTERCEPTION (VEGETATION, ROCKS)
  nms<-c('EIn_H','EIn_L','EIn_rock')
  r_ls<-list()
  for(i in 1:length(nms)){     
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nHours       
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df12<-df
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  
  ## 13: SPLASH EROSION
  nms<-'er'  
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours        
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df13<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 14: INTERCEPTED WATER (STORAGE)
  nms<-'In'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df14<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 15: STORAGE IN FRACTURED ROCK (STORAGE)
  nms<-'FROCK'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df15<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 16: INCOMING LATERAL SUBSURFACE FLOW
  nms<-'Qlat_in'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours       
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df16<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 17: OUTGOING LATERAL SUBSURFACE FLOW
  nms<-'Qlat_out'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours    
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df17<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 18: INTERCEPTED VEGETATION WATER (STORAGE)
  nms<-'Inveg'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df18<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 19: WATER FLUX INCOMING TO THE SOIL
  nms<-'WIS'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df19<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 20: INFILTRATION EXCESS RUNOFF
  nms<-'Rh'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df20<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 21: INFILTRATION
  nms<-'f'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r<-r*nDays  
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df21<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 22: LEAKAGE
  nms<-'Lk'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours    
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df22<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 23: RUNON
  nms<-'q_runon'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours        
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df23<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 24: SATURATION EXCESS RUNOFF
  nms<-'Rd'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df24<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 25: SOIL EVAPORATION 
  nms<-'EG'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours      
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df25<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 26: SNOW DENSITY 
  nms<-'ros'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df26<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 27: SWE 
  nms<-'SWE'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df27<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 28: NET RADIATION 
  nms<-'Rn'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df28<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 29: INCOMING SHORTWAVE RADIATION 
  nms<-'Rsw'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]  
  r<-r[[k]] 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  r[is.na(r)]<-0  ####hack!
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df29<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 30: GROSS PRIMARY PRODUCTION 
  nms<-c('GPP_H','GPP_L')
  r_ls<-list()
  for(i in 1:length(nms)){    
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nDays
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df30<-df
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  
  ## 31: ABOVE GROUND NET PRIMARY PRODUCTION 
  nms<-c('ANPP_H','ANPP_L')
  r_ls<-list()
  for(i in 1:length(nms)){   
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nDays  
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df31<-df
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  ## 32: NET ASSIMILATION 
  nms<-c('An_H','An_L')
  r_ls<-list()
  for(i in 1:length(nms)){ 
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nHours*3600   
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df32<-df
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  
  ## 33: EVAPORATION FROM INTERCEPTION (VEGETATION)
  nms<-c('EIn_H','EIn_L')
  r_ls<-list()
  for(i in 1:length(nms)){       
    fn<-fls[str_detect(fls,paste0('^',nms[i],'___.*','.grd','$'))]
    st<-stack(fn)[[NMONTHS]]
    r_ls[[i]]<-st[[k]]          
    rm(fn,st)
  }
  r<-stack(r_ls)
  r<-calc(r,sum)              
  r<-r*nHours   
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df33<-df
  rm(nms,r_ls,r,pts,ele_v,vals,df)
  
  
  ## 34: EVAPORATION FROM INTERCEPTION (ROCKS) 
  nms<-'EIn_rock'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df34<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 35: EVAPORATION FROM STANDING WATER 
  nms<-'EWAT'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours 
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df35<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 36: EVAPORATION FROM ICE 
  nms<-'EICE'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours    
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df36<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ## 37: EVAPORATION FROM SNOW 
  nms<-'ESN'
  fn<-fls[str_detect(fls,paste0('^',nms,'___.*','.grd','$'))]
  r<-stack(fn)[[NMONTHS]]      
  r<-r[[k]] 
  r<-r*nHours     
  r[is.na(r)]<-0
  extent(r)<-extent(dem)
  pts<-rasterToPoints(r)
  ele_v<-extract(dem,pts[,1:2])
  vals<-na.omit(values(r))
  df<-cbind.data.frame(ele=ele_v,vals=vals)
  df<-df[order(df$ele),]
  df37<-df
  rm(nms,fn,r,pts,ele_v,vals,df)
  
  
  ### MERGE DATA
  df<-cbind.data.frame(eb=df1[,1],ET=df1[,4],P=df2[,2],Imelt=df3[,2],
                       Smelt=df4[,2],Psno=df5[,2],Pliq=df6[,2],
                       V=df7[,2],Vice=df8[,2],ESI=df9[,2],Tr=df10[,2],
                       EGR=df11[,2],EIn=df12[,2],er=df13[,2],In=df14[,2],
                       FROCK=df15[,2],Qlatin=df16[,2],Qlatout=df17[,2],
                       Inveg=df18[,2],WIS=df19[,2],Rh=df20[,2],
                       f=df21[,2],Lk=df22[,2],qrunon=df23[,2],Rd=df24[,2],
                       ESo=df25[,2],ros=df26[,2],SWE=df27[,2],Rn=df28[,2],
                       Rsw=df29[,2],GPP=df30[,2],ANPP=df31[,2],An=df32[,2],
                       EInV=df33[,2],EInR=df34[,2],EWAT=df35[,2],EICE=df36[,2],
                       ESN=df37[,2],
                       asp=df1$asp,slo=df1$slo)
  rm(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10,df11,df12,df13,df14,df15,df16,
     df17,df18,df19,df20,df21,df22,df23,df24,df25,df26,df27,df28,df29,df30,
     df31,df32,df33,df34,df35,df36,df37)
  
  
  #########################################
  ## aggregate average per elevation band
  ##  (used for heatmaps)
  df$eb<-mroundu(df$eb,eb)
  # ## add >7000m a.s.l. to class below
  # ebs<-unique(df$eb)
  # if(max(ebs) > 7000){
  #   df[which(df$eb == max(ebs)),'eb']<-7000
  # }
  df$Area<-resol^2*10^-6
  df_peb<-aggregate(cbind(df$ET,df$P,df$Imelt,df$Smelt,df$Psno,df$Pliq,
                          df$V,df$Vice,df$ESI,df$Tr,df$EGR,df$EIn,df$er,
                          df$In,df$FROCK,df$Qlatin,df$Qlatout,df$Inveg,
                          df$WIS,df$Rh,df$f,df$Lk,df$qrunon,df$Rd,df$ESo,
                          df$ros,df$SWE,df$Rn,df$Rsw,df$GPP,df$ANPP,df$An,
                          df$EInV,df$EInR,df$EWAT,df$EICE,df$ESN),list(df$eb),
                    mean,na.rm=TRUE) 
  names(df_peb)<-c('eb','ET','Precipitation','Ice melt','Snow melt',
                   'Snowfall','Rain','Soil water','Soil ice','E snow-ice','Tr',
                   'E ground','E intercept.','Splash erosion','Total intercept. water',
                   'Water in fract. rocks','Incoming lat. subsurf. flow',
                   'Outgoing lat. subsurf. flow','Intercept. water on vegetation',
                   'Incoming water flux to the soil','Infiltrat. excess runoff',
                   'Infiltration','Bottom leakage soil to bedrock','Run-on',
                   'Saturat. excess runoff','E soil','Snow density','SWE',
                   'Net radiation','Incoming shortwave rad.','Gross primary prod.',
                   'Above ground net primary prod.','Net assimilation',
                   'E intercept. (vegetation)','E intercept. (rock)',
                   'E standing water','E ice','E snow')
  df_peb$Area<-aggregate(df$Area,list(df$eb),sum,na.rm=TRUE)[,2]
  df_peb$t<-k
  df_peb$mon<-mon
  df_peb$Date<-Dates[k]
  
  # write table
  fn<-path_data&'/TABLES/avg-p'&eb&'m__'&yearmon&'.txt'
  write.table(df_peb,fn,dec='.',sep='\t',quote=FALSE,
              col.names=TRUE,row.names=FALSE)
  
  df_peb_ls[[k]]<-df_peb
  
  
  ##############################################
  ## aggregate cumulative sum per elevation band
  ##  (only used for balance lineplots so far...)
  df_cumpeb<-aggregate(cbind(df$ET,df$P,df$Imelt,df$Smelt,df$Psno,df$Pliq,
                          df$V,df$Vice,df$ESI,df$Tr,df$EGR,df$EIn,df$er,
                          df$In,df$FROCK,df$Qlatin,df$Qlatout,df$Inveg,
                          df$WIS,df$Rh,df$f,df$Lk,df$qrunon,df$Rd,df$ESo),
                       list(df$eb),
                    mean,na.rm=TRUE) 
  df$n<-1
  df_cumpeb$n<-aggregate(df$n,list(df$eb),sum,na.rm=TRUE)[,2] 
  names(df_cumpeb)<-c('eb','ET','Precipitation','Ice melt','Snow melt',
                   'Snowfall','Rain','Soil water','Soil ice','E snow-ice','Tr',
                   'E ground','E intercept.','Splash erosion','Total intercept. water',
                   'Water in fract. rocks','Incoming lat. subsurf. flow',
                   'Outgoing lat. subsurf. flow','Intercept. water on vegetation',
                   'Incoming water flux to the soil','Infiltrat. excess runoff',
                   'Infiltration','Bottom leakage soil to bedrock','Run-on',
                   'Saturat. excess runoff','E soil','n')
  idx<-which(names(df_cumpeb) %in% c('eb','n'))
  df_cumpeb[,-idx]<-df_cumpeb[,-idx]*df_cumpeb$n
  df_cumpeb$t<-k
  df_cumpeb$mon<-mon
  df_cumpeb$Date<-Dates[k]
  
  # write table
  fn<-path_data&'/TABLES/cum-p'&eb&'m__'&yearmon&'.txt'
  write.table(df_cumpeb,fn,dec='.',sep='\t',quote=FALSE,
              col.names=TRUE,row.names=FALSE)
  
  df_cumpeb_ls[[k]]<-df_cumpeb
  rm(df,df_peb,df_cumpeb,fn)
}

df_peb<-as.data.frame(do.call('rbind',df_peb_ls))
df_cumpeb<-as.data.frame(do.call('rbind',df_cumpeb_ls))



##derive change in variables over time
ebs<-unique(df_peb$eb)
df_peb$'d liq. Soil water'<-NA
df_peb$'d Soil water'<-NA
df_peb$'d lat. Subsurf. flow'<-NA

df_cumpeb$'d liq. Soil water'<-NA
df_cumpeb$'d Soil water'<-NA
df_cumpeb$'d lat. Subsurf. flow'<-NA
for(i in 1:length(ebs)){
  eb<-ebs[i]
  
  ##############
  ### DF_PEB ###
  idx<-which(df_peb$eb == eb)
  
  ##dV_liq
  v<-df_peb[idx,'Soil water']
  v<-c(0,diff(v))
  df_peb[idx,'d liq. Soil water']<-v
  rm(v)
  
  ##dV_all
  v<-df_peb[idx,'Soil water']+df_peb[idx,'Soil ice']
  v<-c(0,diff(v))
  df_peb[idx,'d Soil water']<-v
  rm(v)
  
  ##dQlat
  v<-df_peb[idx,'Incoming lat. subsurf. flow']-df_peb[idx,'Outgoing lat. subsurf. flow']
  df_peb[idx,'d lat. Subsurf. flow']<-v
  rm(v)
  
  rm(idx)
  
  #################
  ### DF_CUMPEB ###
  idx<-which(df_cumpeb$eb == eb)
  
  ##dV_liq
  v<-df_cumpeb[idx,'Soil water']
  v<-c(0,diff(v))
  df_cumpeb[idx,'d liq. Soil water']<-v
  rm(v)
  
  ##dV_all
  v<-df_cumpeb[idx,'Soil water']+df_cumpeb[idx,'Soil ice']
  v<-c(0,diff(v))
  df_cumpeb[idx,'d Soil water']<-v
  rm(v)
  
  ##dQlat
  v<-df_cumpeb[idx,'Incoming lat. subsurf. flow']-df_cumpeb[idx,'Outgoing lat. subsurf. flow']
  df_cumpeb[idx,'d lat. Subsurf. flow']<-v
  rm(v)
  
  rm(idx,eb)
}



newdir<-path_data&'/PLOTS/VEG'
dir.create(newdir)
#===============================================================================
# +++ PLOT AVERAGE OVER ELEVATION (1 LINE PER YEAR) +++ 
# PANEL 1: GPP
# PANEL 2: ANPP
# PANEL 3: TRANSPIRATION
# PANEL 4: ASSIMILATION
#===============================================================================
mon_idx<-5:9 ##months to consider
xlab<-'Elevation [m a.s.l.]'
lwd<-1.5

### PANEL 1 ###
VAR<-'Gross primary prod.'
ylab<-expression(paste('[gC',m^-2,' day]'))
tit<-VAR
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms == VAR)]
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual mean per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),mean) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
###
pl1<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl1
rm(VAR,Date,eb,nms,ylab,rng,tit,DF,DF_ls)


### PANEL 2 ###
VAR<-'Above ground net primary prod.'
ylab<-expression(paste('[gC',m^-2,' day]'))
tit<-VAR
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms == VAR)]
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual mean per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),mean) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
###
pl2<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl2
rm(VAR,Date,eb,nms,ylab,rng,tit,DF,DF_ls)


### PANEL 3 ###
VAR<-'Tr'
ylab<-expression(paste('[mm]'))
tit<-'Transpiration'
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms == VAR)]
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual sum per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
###
pl3<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl3
rm(VAR,Date,eb,nms,ylab,rng,tit,DF,DF_ls)


### PANEL 4 ###
VAR<-'Net assimilation'
ylab<-expression(paste('[mol ',CO[2],' ',m^-2,']'))
tit<-VAR
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms == VAR)] * 10^-6 ###micromol to mol 
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual mean per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),mean) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),2),mroundu(max(DF$value),2))
###
pl4<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_text(colour='black',size=base_size),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl4
rm(VAR,Date,eb,nms,ylab,rng,tit,DF,DF_ls)

toptit<-textGrob(DOMAIN&' ('&month.abb[mon_idx[1]]&'-'&month.abb[mon_idx[length(mon_idx)]]&')',
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3,pl4),nrow=4,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/annual_pe_GPP-ANPP-Tr-An_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=30,units='cm')
graphics.off()
rm(fn,toptit,pl)



#===============================================================================
# +++ PLOT AVERAGE OVER ELEVATION (1 LINE PER YEAR) +++ 
# PANEL 1: P
# PANEL 2: ET
# PANEL 3: Tr
# PANEL 4: P-ET
#===============================================================================
mon_idx<-1:12 ##months to consider
xlab<-'Elevation [m a.s.l.]'
lwd<-1.5

### PANEL 1 ###
VAR<-'Precipitation'
ylab<-expression(paste('[mm]'))
tit<-VAR
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms == VAR)]
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual sum per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
# #### hack!!
# if(any(10:12 %in% mon_idx)){
#   print('discarded last (incomplete) year')
#   DF<-DF[-which(DF$Date == '2022'),]
# }
# ###
pl1<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl1
rm(VAR,Date,eb,nms,ylab,rng,tit,DF,DF_ls)


### PANEL 2 ###
VAR<-'ET'
ylab<-expression(paste('[mm]'))
tit<-VAR
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms == VAR)]
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual sum per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
# #### hack!!
# if(any(10:12 %in% mon_idx)){
#   print('discarded last (incomplete) year')
#   DF<-DF[-which(DF$Date == '2022'),]
# }
# ###
pl2<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl2
rm(VAR,Date,eb,nms,ylab,rng,tit,DF,DF_ls)


### PANEL 3 ###
VAR<-'Tr'
ylab<-expression(paste('[mm]'))
tit<-'Transpiration'
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms == VAR)]
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual sum per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
# #### hack!!
# if(any(10:12 %in% mon_idx)){
#   print('discarded last (incomplete) year')
#   DF<-DF[-which(DF$Date == '2022'),]
# }
# ###
pl3<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl3
rm(VAR,Date,eb,nms,ylab,rng,tit,DF,DF_ls)


### PANEL 4 ###
VAR1<-'Precipitation'
VAR2<-'ET'
ylab<-expression(paste('[mm]'))
tit<-'P-ET'
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms %in% c(VAR1,VAR2))]
DF<-DF$Precipitation-DF$ET
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual sum per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),100),mroundu(max(DF$value),100))
# #### hack!!
# if(any(10:12 %in% mon_idx)){
#   print('discarded last (incomplete) year')
#   DF<-DF[-which(DF$Date == '2022'),]
# }
# ###
pl4<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  geom_hline(yintercept=0,colour='black',linewidth=lwd/2)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_text(colour='black',size=base_size),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl4
rm(VAR1,VAR2,Date,eb,nms,ylab,rng,tit,DF,DF_ls)

toptit<-textGrob(DOMAIN&' ('&month.abb[mon_idx[1]]&'-'&month.abb[mon_idx[length(mon_idx)]]&')',
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3,pl4),nrow=4,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/annual_pe_Pr-ET-Tr-PrminusET_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=30,units='cm')
graphics.off()
rm(fn,toptit,pl)



#===============================================================================
# +++ PLOT AVERAGE OVER ELEVATION (1 LINE PER YEAR) +++ 
# PANEL 1: BLUE WATER
# PANEL 2: GREEN WATER
# PANEL 3: WHITE WATER
#===============================================================================
mon_idx<-1:12 ##months to consider
xlab<-'Elevation [m a.s.l.]'
lwd<-1.5
VAR<-c('E soil','Intercept. water on vegetation','E intercept. (rock)',
       'E standing water','Tr','E snow-ice','Precipitation') 
ylab<-expression(paste('[mm]'))

### PANEL 1 ###
tit<-'Blue water flux'
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms %in% VAR)]
DF$GR<-DF$Tr+DF$`Intercept. water on vegetation`+DF$`E soil`+DF$`E intercept. (rock)`+DF$`E standing water`
DF$WH<-DF$`E snow-ice`
DF$BL<-DF$Precipitation-DF$GR-DF$WH
# DF[DF$BL < 0,'BL']<-0
DF<-DF$BL
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual sum per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
# #### hack!!
# if(any(10:12 %in% mon_idx)){
#   print('discarded last (incomplete) year')
#   DF<-DF[-which(DF$Date == '2022'),]
# }
# ###
pl1<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  geom_hline(yintercept=0,linewidth=lwd/2)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl1
rm(Date,eb,nms,rng,tit,DF,DF_ls)


### PANEL 2 ###
tit<-'Green water flux'
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms %in% VAR)]
DF$GR<-DF$Tr+DF$`Intercept. water on vegetation`+DF$`E soil`+DF$`E intercept. (rock)`+DF$`E standing water`
DF<-DF$GR
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual sum per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
# #### hack!!
# if(any(10:12 %in% mon_idx)){
#   print('discarded last (incomplete) year')
#   DF<-DF[-which(DF$Date == '2022'),]
# }
# ###
pl2<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  geom_hline(yintercept=0,linewidth=lwd/2)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl2
rm(Date,eb,nms,rng,tit,DF,DF_ls)


### PANEL 3 ###
tit<-'White water flux'
DF<-df_peb
nms<-names(DF)
Date<-DF$Date
eb<-DF$eb
DF<-DF[,which(nms %in% VAR)]
DF$WH<-DF$`E snow-ice`
DF<-DF$WH
DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
DF$mon<-month(DF$Date)
DF$yr<-as.character(year(DF$Date))
DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
#derive annual sum per elevation band
DF_ls<-list()
for(i in 1:length(ebs)){
  eb<-ebs[i]
  df<-DF[which(DF$eb == eb),1:3]
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
  df$eb<-eb
  df$Date<-as.character(year(df$Date))
  DF_ls[[i]]<-df
  rm(df)
}
DF<-cbind.data.frame(do.call('rbind',DF_ls))
rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
# #### hack!!
# if(any(10:12 %in% mon_idx)){
#   print('discarded last (incomplete) year')
#   DF<-DF[-which(DF$Date == '2022'),]
# }
# ###
pl3<-ggplot(data=DF)+
  geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
  geom_hline(yintercept=0,linewidth=lwd/2)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=tit,x=xlab)+
  scale_color_nejm()+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_text(colour='black',size=base_size),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey'),
    panel.grid.major.y = element_line(size=0.75,colour='grey'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey'),
    panel.grid.major.x = element_line(size=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl3
rm(Date,eb,nms,rng,tit,DF,DF_ls)



# ### PANEL 4 ###
# tit<-'Blue/Green ratio'
# DF<-df_peb
# nms<-names(DF)
# Date<-DF$Date
# eb<-DF$eb
# DF<-DF[,which(nms %in% VAR)]
# DF$GR<-DF$Tr+DF$`Intercept. water on vegetation`+DF$`E soil`+DF$`E intercept. (rock)`+DF$`E standing water`
# DF$WH<-DF$`E snow-ice`
# DF$BL<-DF$Precipitation-DF$GR-DF$WH
# DF<-DF$BL/DF$GR
# DF<-cbind.data.frame(Date=Date,eb=eb,value=DF)
# DF$mon<-month(DF$Date)
# DF$yr<-as.character(year(DF$Date))
# DF<-DF[which(DF$mon %in% mon_idx),]   ##select only months Apr-Oct   
# #derive annual sum per elevation band
# DF_ls<-list()
# for(i in 1:length(ebs)){
#   eb<-ebs[i]
#   df<-DF[which(DF$eb == eb),1:3]
#   df<-aggregate(df[,-1],list(Date=cut(df$Date,'year')),sum) 
#   df$eb<-eb
#   df$Date<-as.character(year(df$Date))
#   DF_ls[[i]]<-df
#   rm(df)
# }
# DF<-cbind.data.frame(do.call('rbind',DF_ls))
# rng<-c(mroundd(min(DF$value),10),mroundu(max(DF$value),10))
# # #### hack!!
# # if(any(10:12 %in% mon_idx)){
# #   print('discarded last (incomplete) year')
# #   DF<-DF[-which(DF$Date == '2022'),]
# # }
# # ###
# pl4<-ggplot(data=DF)+
#   geom_line(aes(x=eb,y=value,group=Date,colour=Date),linewidth=lwd)+
#   geom_hline(yintercept=0,linewidth=lwd/2)+
#   scale_y_continuous(name=ylab,limits=rng)+
#   labs(title=tit,x=xlab)+
#   scale_color_nejm()+
#   theme_gray(base_size=base_size,base_family=base_family)+
#   theme(
#     legend.position="right",
#     legend.title = element_blank(),
#     axis.ticks = element_blank(),
#     axis.text.y = element_text(colour='black',size=base_size),
#     axis.title.y = element_text(colour='black',size=base_size),
#     axis.text.x = element_text(colour='black',size=base_size),
#     axis.title.x = element_blank(),
#     panel.background = element_rect(fill='white',colour='white'),
#     panel.grid.minor.y = element_line(size=0.1,colour='grey'),
#     panel.grid.major.y = element_line(size=0.75,colour='grey'),
#     panel.grid.minor.x = element_line(size=0.1,colour='grey'),
#     panel.grid.major.x = element_line(size=0.75,colour='grey'),
#     strip.background = element_blank(),
#     strip.text.y = element_blank(),
#     plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
#   ) 
# pl4
# rm(Date,eb,nms,rng,tit,DF,DF_ls)

toptit<-textGrob(DOMAIN&' ('&month.abb[mon_idx[1]]&'-'&month.abb[mon_idx[length(mon_idx)]]&')',
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3),nrow=3,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/annual_pe_BGW_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=27,units='cm')
graphics.off()
rm(fn,toptit,pl)


###
print('+++END SCRIPT <<PLOT ELE-MON>>+++   '&outnm)



