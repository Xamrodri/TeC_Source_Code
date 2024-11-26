################################################################################
# T&C OUTPUT :
#  COMPARE SOIL MOISTURE | PLOT
#
# 
# TODO: -
#       -
#
# NEW:  -
#       -
#
# 2023/09/25
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
# outnm<-'v20230906_Aaretal_2557days'
# outnm<-'v20230906_Aargauer_Rhein_2557days'
# outnm<-'v20230906_Berner_Oberland_2557days'
# outnm<-'v20230906_Bodensee_2557days'
# outnm<-'v20230906_Engadin_2557days'
# outnm<-'v20230906_Glatt_2557days'
# outnm<-'v20230906_Hinterrhein_2557days'
# outnm<-'v20230906_Jura_2557days'
# outnm<-'v20230906_Leman_2557days'
# outnm<-'v20230906_Limmat_2557days'
# outnm<-'v20230906_Neuchatel_2557days'
# outnm<-'v20230906_Oberwallis_2557days'
# outnm<-'v20230906_Reuss_2557days'
# outnm<-'v20230906_Rheintal_2557days'
# outnm<-'v20230906_Saane_2557days'
# outnm<-'v20230906_Schaffhausen_2557days'
# outnm<-'v20230906_Ticino_2557days'
# outnm<-'v20230906_Thur_2557days'
# outnm<-'v20230906_Toess_2557days'
# outnm<-'v20230906_Unterwallis_369days'
# outnm<-'v20230906_Vorderrhein_2557days'
# outnm<-'v20230906_Zuercher_Rhein_2557days'
###
# outnm<-'v20230906_Test-Aletsch_2557days'
# outnm<-'v20230906_Test-Grossbach_2557days'
# outnm<-'v20230906_Test-Poschiavino_2557days'
# outnm<-'v20230906_Test-Rein_da_Sumvitg_2557days'
outnm<-'v20230906_Test-Rietholzbach_2557days'
# outnm<-'v20230906_Test-Roggiasca_2557days'
# outnm<-'v20230906_Test-Rosegbach_2557days'

## define timespan & months to analyze
ts_init<-'2016-10-01 02:00:00'
ts_final<-'2022-09-30 23:00:00'
# ts_init<-'2015-10-01 02:00:00'
# ts_final<-'2022-08-01 23:00:00'


#===============================================================================
# PATHS
#===============================================================================
#X1
if(Sys.getenv('computername') == 'PASCAL-THINK'){ 
  root<-'D:'                                    ## Root
  path_onst<-'C:/Users/Pascal/switchdrive/WSL'  ## Path to online storage
}
#T490
if(Sys.getenv('computername') == 'CHARA'){ 
  root<-'E:'                                   ## Root
  path_onst<-'C:/Users/Buri/switchdrive/WSL'   ## Path to online storage
}

## Path to T&C data
path_data<-root&'/T&C/POSTPROCESSED/DISTRIBUTED/'&outnm

## Path to code scripts etc.
path_code<-path_onst&'/T&C/Code/POSTprocessing/LAN2019'

## Path to functions file
path_func<-path_onst&'/T&C/Code/Functions/Functions.r'

## Path to GIS data
path_gis<-root&'/QGIS'

## define folder of reference soil moisture data
refSMfo<-root&'/Meteodata/CH/SoilMoisture/from_AW'


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
library(lubridate)    #floor_date(), hour() 
library(networkD3)    #Sankey diagram
library(dplyr)        #Sankey diagram
library(webshot)      #save HTML screenshot
library(Metrics)      #rmse()
library(topmodel)     #NSeff()
library(stringr)      #str_replace
library(oce)          #despike()
library(zoom)         #zm()
library(extrafont)     #loadfonts()
library(tidyverse)
library(lubridate)
library(assertthat)
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
# TIME PERIOD OF INTEREST
#===============================================================================
ts_init<-as.POSIXct(ts_init,format='%Y-%m-%d %H:%M:%S')
ts_final<-as.POSIXct(ts_final,format='%Y-%m-%d %H:%M:%S')

PERIOD_SH<-format(ts_init,"%b%Y")&'-'&format(ts_final-86400,"%b%Y")
PERIOD_LO<-format(ts_init,"%d %b %Y")&' - '&format(ts_final-86400,"%d %b %Y")


#===============================================================================
# FIND SIMULATED SOIL MOISTURE STATIONS WITHIN DOMAIN
#===============================================================================
pois<-st_read(path_data&'/SHAPEFILES/POIs.shp',quiet=T)
sim_ids<-which(pois$ID2 == 'SM')
if(length(sim_ids) == 0 && nrow(pois) == 1){sim_ids<-1}  ##Test-domains
pois<-pois[sim_ids,]


#===============================================================================
# FIND SPATIAL RESOLUTION & PROJECTION
#===============================================================================
dem<-raster(path_data&'/RASTERS/INIT-COND/DTM.tif')
resol<-res(dem)[1]
projec<-projection(dem)


#===============================================================================
# READ OBSERVED SOIL MOISTURE DATA
#===============================================================================
setwd(refSMfo)
obs<-readRDS(file='02_table_soilmoisture.rds')
ch_stations<-st_transform(ch_stations,crs=projec)  ##reproject to CH
st_crs(ch_stations)<-st_crs(pois)
dist<-st_distance(pois,ch_stations)
idx<-apply(dist,1,which.min)
obs<-ch_stations[idx,]
## check
# dist[1,idx[1]]
# dist[2,idx[2]]
rm(dist,idx)

df_ls<-list()
missing_sta<-vector()
for(i in 1:nrow(obs)){
  STA_ID<-as.character(obs$id[i])
  STA_NM<-as.character(obs$name[i])
  STA_NM<-str_split(STA_NM,' - ')[[1]][1]
  
  #read data
  fn<-STA_ID&'.daily.mean.dat'
  files<-list.files()
  
  if(length(which(files == fn)) == 0){
    missing_sta[i]<-i
    print('Station  << '&STA_NM&' >>  missing!')
  }
  
  if(length(which(files == fn)) == 1){
    df<-read.table(fn,header=T,sep=' ',dec='.',na.strings='NaN') 
    df$date<-as.POSIXct(df$date,format='%Y-%m-%d',tz=TZ)
    
    #subset to initial timestep
    idxrm<-which(df$date < ts_init)
    if(length(idxrm) > 0){df<-df[-idxrm,]}
    # df$date[1:10]
    rm(idxrm)
    
    #subset to final timestep
    idxrm<-which(df$date > ts_final)
    if(length(idxrm) > 0){df<-df[-idxrm,]}
    # df$date[nrow(df)]
    rm(idxrm)
    
    colnames(df)<-c('Date',STA_NM)
    df_ls[[i]]<-df
    rm(df)
  }
}
if(length(missing_sta) > 0){
  obs<-obs[-missing_sta,]
  df_ls<-df_ls[[-missing_sta]]
  pois<-pois[-missing_sta,]
}
obsdat<-as.data.frame(do.call('cbind',df_ls))
rmidx<-which(names(obsdat) == 'Date')[-1]
if(length(rmidx) > 0){obsdat<-obsdat[,-rmidx]}
rm(df_ls,rmidx)

### HACK !!!
idx<-which(sapply(obsdat,is.character))
if(length(idx) > 0){
  ifelse(length(idx) == 1,
         obsdat[,idx]<-as.numeric(obsdat[,idx]),
         obsdat[,idx]<-lapply(obsdat[,idx],as.numeric)
  )
}
#####


#===============================================================================
# READ MODELLED RUNOFF DATA
#===============================================================================
setwd(path_data&'/TABLES')

df_ls<-list()
for(i in 1:nrow(pois)){
  STA_NM<-pois$ID[i]&'-'&pois$ID2[i]
  ELEV<-as.character(pois$elevation[i])
  
  #read data
  fn<-'TrackedPixel___'&STA_NM&'___'&ELEV&'m.txt'
  df<-read.table(fn,header=T,sep='\t',dec='.',na.strings='NaN') 
  ##correct for midnight hour
  df$Date<-as.character(df$Date)
  idx<-which(nchar(df$Date) < 19)
  df$Date[idx]<-df$Date[idx]&' 00:00:00'
  rm(idx)
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S',tz=TZ)
  df<-df[,which(names(df) %in% c('Date','QpointC'))]
  
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'day')),mean) #derive daily mean
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d')
  
  df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period
  
  ##convert to [m3 s-1]
  DTH<-1 ; DT<-3600;       ## dth= 1 h; dt= 3600 s 
  df$x<-(df$x/DTH)*resol^2/(DT*1000)   ##ROUTING_MODULE: l.163
  
  colnames(df)<-c('Date',STA_NM)
  df_ls[[i]]<-df
  rm(df)
}
simdat<-as.data.frame(do.call('cbind',df_ls))
rmidx<-which(names(simdat) == 'Date')[-1]
if(length(rmidx) > 0){simdat<-simdat[,-rmidx]}
rm(df_ls,rmidx)


newdir<-path_data&'/PLOTS'
dir.create(newdir)
dir.create(newdir&'/Q')
#===============================================================================
#===============================================================================
#===============================================================================
# PLOT DAILY RUNOFF @ POIs
#===============================================================================
#===============================================================================
#===============================================================================

#### LOOP OVER POIs #####
for(i in 1:nrow(pois)){
  ##extract names & tables
  nm<-as.character(pois$name[i])
  ele<-as.character(pois$elevation[i])
  
  df<-merge(obsdat[,c(1,i+1)],simdat[,c(1,i+1)],by='Date',all.y=T)
  names(df)<-c('DATE','OBS','SIM')
  
  ##############################
  ### PLOT RUNOFF COMPARISON ###
  ##############################
  
  ## PERFORMANCE ##
  ##calculate RMSE
  comp<-df[complete.cases(df),] 
  which(is.na(comp))
  RMSE<-round(rmse(comp$OBS,comp$SIM),2)
  ##calculate R^2
  R2<-round(cor(comp$OBS,comp$SIM)^2,3)
  ##calculate Nash-Sutcliffe
  NSE<-round(NSeff(comp$OBS,comp$SIM),3)
  ##calculate totals
  avg_obs<-round(mean(comp$OBS),2)
  avg_sim<-round(mean(comp$SIM),2)
  rm(comp)
  #################
  
  df<-melt(df,'DATE')
  tit<-'Runoff '&nm&' \n ('&ele&' m a.s.l.)'
  capt<-''
  xlab<-''
  ylab<-expression(paste('[',m^3,' ',s^-1,']'))
  capt<-paste('R2 = ',R2,'   |   RMSE = ',RMSE,
              'm3/s-1   |   NSE = ',NSE,'\n',
              'Average [m3/s-1]: ',avg_obs,' (obs.) | ',avg_sim,' (sim.)')
  leg<-''
  lims<-range(df$DATE)
  
  pl<-ggplot(df)+
    # annotate("rect",fill="red",alpha=0.2,xmin=dr1_init,xmax=dr1_final,
    #          ymin=-Inf,ymax=Inf)+
    geom_line(aes(x=DATE,y=value,group=variable,color=variable),linewidth=1)+
    scale_color_jco()+
    labs(title=tit,caption=capt,x=xlab,y=ylab,fill=leg)+
    theme_gray(base_size=base_size,base_family=base_family)+
    theme(
      legend.position = 'right',
      legend.direction = 'vertical',
      # legend.position = 'none',
      legend.title = element_blank(),
      legend.text = element_text(size=base_size),
      legend.key=element_rect(fill='white'),
      axis.ticks = element_blank(),
      axis.text.y = element_text(colour='black',size=base_size-2),
      axis.text.x = element_text(colour='black',size=base_size-2),
      panel.background = element_rect(fill='white',colour='white'),
      panel.grid.minor.y = element_line(size=0.1,colour='grey'),
      panel.grid.major.y = element_line(size=0.75,colour='grey'),
      panel.grid.major.x = element_line(size=0.75,colour='grey'),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      plot.margin=unit(c(1,1,0,0),'lines') #top/right/bottom/left
    ) 
  pl
  
  
  fn<-newdir&'/Q/Q-daily_'&nm&'_'&PERIOD_SH&'.png'
  ggsave(fn,pl,width=20,height=15,units='cm')
  graphics.off()
  rm(df,tit,capt,avg_obs,avg_sim,xlab,ylab,leg,lims,pl,fn)
  rm(nm,ele,DATE,SIM,OBS,RMSE,NSE,R2)
}



#===============================================================================
#===============================================================================
#===============================================================================
# PLOT WEEKLY RUNOFF @ POIs
#===============================================================================
#===============================================================================
#===============================================================================

#### LOOP OVER POIs #####
for(i in 1:nrow(pois)){
  ##extract names & tables
  nm<-as.character(pois$name[i])
  ele<-as.character(pois$elevation[i])
  
  df<-merge(obsdat[,c(1,i+1)],simdat[,c(1,i+1)],by='Date',all.y=T)
  names(df)<-c('DATE','OBS','SIM')
  
  df<-aggregate(df[,-1],list(DATE=cut(df$DATE,'week')),mean) #derive weekly mean
  df$DATE<-as.POSIXct(df$DATE,format='%Y-%m-%d')
  
  ##############################
  ### PLOT RUNOFF COMPARISON ###
  ##############################
  
  ## PERFORMANCE ##
  ##calculate RMSE
  comp<-df[complete.cases(df),] 
  RMSE<-round(rmse(comp$OBS,comp$SIM),2)
  ##calculate R^2
  R2<-round(cor(comp$OBS,comp$SIM)^2,3)
  ##calculate Nash-Sutcliffe
  NSE<-round(NSeff(comp$OBS,comp$SIM),3)
  ##calculate totals
  avg_obs<-round(mean(comp$OBS),2)
  avg_sim<-round(mean(comp$SIM),2)
  rm(comp)
  #################
  
  df<-melt(df,'DATE')
  tit<-'Runoff '&nm&' \n ('&ele&' m a.s.l.)'
  capt<-''
  xlab<-''
  ylab<-expression(paste('[',m^3,' ',s^-1,']'))
  capt<-paste('R2 = ',R2,'   |   RMSE = ',RMSE,
              'm3/s-1   |   NSE = ',NSE,'\n',
              'Average [m3/s-1]: ',avg_obs,' (obs.) | ',avg_sim,' (sim.)')
  leg<-''
  lims<-range(df$DATE)
  
  pl<-ggplot(df)+
    # annotate("rect",fill="red",alpha=0.2,xmin=dr1_init,xmax=dr1_final,
    #          ymin=-Inf,ymax=Inf)+
    geom_line(aes(x=DATE,y=value,group=variable,color=variable),linewidth=1)+
    scale_color_jco()+
    labs(title=tit,caption=capt,x=xlab,y=ylab,fill=leg)+
    theme_gray(base_size=base_size,base_family=base_family)+
    theme(
      legend.position = 'right',
      legend.direction = 'vertical',
      # legend.position = 'none',
      legend.title = element_blank(),
      legend.text = element_text(size=base_size),
      legend.key=element_rect(fill='white'),
      axis.ticks = element_blank(),
      axis.text.y = element_text(colour='black',size=base_size-2),
      axis.text.x = element_text(colour='black',size=base_size-2),
      panel.background = element_rect(fill='white',colour='white'),
      panel.grid.minor.y = element_line(size=0.1,colour='grey'),
      panel.grid.major.y = element_line(size=0.75,colour='grey'),
      panel.grid.major.x = element_line(size=0.75,colour='grey'),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      plot.margin=unit(c(1,1,0,0),'lines') #top/right/bottom/left
    ) 
  pl
  
  
  fn<-newdir&'/Q/Q-weekly_'&nm&'_'&PERIOD_SH&'.png'
  ggsave(fn,pl,width=20,height=15,units='cm')
  graphics.off()
  rm(df,tit,capt,avg_obs,avg_sim,xlab,ylab,leg,lims,pl,fn)
  rm(nm,ele,DATE,SIM,OBS,RMSE,NSE,R2)
}



#===============================================================================
#===============================================================================
#===============================================================================
# PLOT MONTHLY RUNOFF @ POIs
#===============================================================================
#===============================================================================
#===============================================================================

#### LOOP OVER POIs #####
for(i in 1:nrow(pois)){
  ##extract names & tables
  nm<-as.character(pois$name[i])
  ele<-as.character(pois$elevation[i])
  
  df<-merge(obsdat[,c(1,i+1)],simdat[,c(1,i+1)],by='Date',all.y=T)
  names(df)<-c('DATE','OBS','SIM')
  
  df<-aggregate(df[,-1],list(DATE=cut(df$DATE,'month')),mean) #derive monthly mean
  df$DATE<-as.POSIXct(df$DATE,format='%Y-%m-%d')
  
  ##############################
  ### PLOT RUNOFF COMPARISON ###
  ##############################
  
  ## PERFORMANCE ##
  ##calculate RMSE
  comp<-df[complete.cases(df),] 
  RMSE<-round(rmse(comp$OBS,comp$SIM),2)
  ##calculate R^2
  R2<-round(cor(comp$OBS,comp$SIM)^2,3)
  ##calculate Nash-Sutcliffe
  NSE<-round(NSeff(comp$OBS,comp$SIM),3)
  avg_obs<-round(mean(comp$OBS),2)
  avg_sim<-round(mean(comp$SIM),2)
  rm(comp)
  #################
  
  df<-melt(df,'DATE')
  tit<-'Runoff '&nm&' \n ('&ele&' m a.s.l.)'
  capt<-''
  xlab<-''
  ylab<-expression(paste('[',m^3,' ',s^-1,']'))
  capt<-paste('R2 = ',R2,'   |   RMSE = ',RMSE,
              'm3/s-1   |   NSE = ',NSE,'\n',
              'Average [m3/s-1]: ',avg_obs,' (obs.) | ',avg_sim,' (sim.)')
  leg<-''
  lims<-range(df$DATE)
  
  pl<-ggplot(df)+
    # annotate("rect",fill="red",alpha=0.2,xmin=dr1_init,xmax=dr1_final,
    #          ymin=-Inf,ymax=Inf)+
    geom_line(aes(x=DATE,y=value,group=variable,color=variable),linewidth=1)+
    scale_color_jco()+
    labs(title=tit,caption=capt,x=xlab,y=ylab,fill=leg)+
    theme_gray(base_size=base_size,base_family=base_family)+
    theme(
      legend.position = 'right',
      legend.direction = 'vertical',
      # legend.position = 'none',
      legend.title = element_blank(),
      legend.text = element_text(size=base_size),
      legend.key=element_rect(fill='white'),
      axis.ticks = element_blank(),
      axis.text.y = element_text(colour='black',size=base_size-2),
      axis.text.x = element_text(colour='black',size=base_size-2),
      panel.background = element_rect(fill='white',colour='white'),
      panel.grid.minor.y = element_line(size=0.1,colour='grey'),
      panel.grid.major.y = element_line(size=0.75,colour='grey'),
      panel.grid.major.x = element_line(size=0.75,colour='grey'),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      plot.margin=unit(c(1,1,0,0),'lines') #top/right/bottom/left
    ) 
  pl
  
  
  fn<-newdir&'/Q/Q-monthly_'&nm&'_'&PERIOD_SH&'.png'
  ggsave(fn,pl,width=20,height=15,units='cm')
  graphics.off()
  rm(df,tit,capt,avg_obs,avg_sim,xlab,ylab,leg,lims,pl,fn)
  rm(nm,ele,DATE,SIM,OBS,RMSE,NSE,R2)
}




###
print('+++END SCRIPT "PLOT VALIDATE SOILMOISTURE"+++   '&outnm)




