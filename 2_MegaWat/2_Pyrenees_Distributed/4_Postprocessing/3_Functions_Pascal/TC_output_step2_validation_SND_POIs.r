################################################################################
# T&C OUTPUT :
#  COMPARE SNOW DEPTH | PLOT
#
# 
# TODO: -
#       -
#
# NEW:  -
#       -
#
# 2023/10/31
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


## define timespan & months to analyze
ts_init<-'2015-10-01 02:00:00'
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
if(Sys.getenv('computername') == 'WSL26871'){ 
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

## define folder of reference snow data
refSNDfo<-root&'/Meteodata/CH/IMIS/Daily_data'


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
# FIND SIMULATED IMIS-STATIONS WITHIN DOMAIN
#===============================================================================
pois<-st_read(path_data&'/SHAPEFILES/POIs.shp',quiet=T)
sim_ids<-which(pois$ID2 %in% 'IMI')
if(length(sim_ids) == 0 && nrow(pois) == 1){sim_ids<-1}  ##Test-domains
pois<-pois[sim_ids,]


#===============================================================================
# READ OBSERVED SND DATA
#===============================================================================
setwd(refSNDfo)
df<-read.table('daily_snow_2016-2023.csv',sep=',',header=T)
df$Date<-as.POSIXct(df$measure_date,format='%Y-%m-%d %H:%M:%S',tz=TZ)

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
obsdat<-df
rm(df)


#===============================================================================
# READ MODELLED SND DATA
#===============================================================================
setwd(path_data&'/TABLES')
n<-nrow(pois)

df_ls<-list()
for(i in 1:n){
  STA_NM<-pois$ID[i]&'-'&pois$ID2[i]
  ELEV<-as.character(pois$elevation[i])
  
  print('Reading modelled SND-data at <<'&STA_NM&' ('&ELEV&'m)>>  ('&i&'/'&n&')')
  #read data
  fn<-'TrackedPixel___'&STA_NM&'___'&ELEV&'m.txt'
  df<-read.table(fn,header=T,sep='\t',dec='.',na.strings='NaN') 
  ##correct for midnight hour
  df$Date<-as.character(df$Date)
  idx<-which(nchar(df$Date) < 19)
  df$Date[idx]<-df$Date[idx]&' 00:00:00'
  rm(idx)
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S',tz=TZ)
  df<-df[,which(names(df) %in% c('Date','SND'))]
  
  df$SND<-df$SND*100 ##m to cm
  
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'day')),median) #derive daily median
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d')
  
  df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period
  
  colnames(df)<-c('Date',STA_NM)
  df_ls[[i]]<-df
  rm(df)
}
simdat<-as.data.frame(do.call('cbind',df_ls))
rmidx<-which(names(simdat) == 'Date')[-1]
if(length(rmidx) > 0){simdat<-simdat[,-rmidx]}
rm(df_ls,rmidx)
names(simdat)<-str_replace(names(simdat),'-IMI','')

newdir<-path_data&'/PLOTS'
dir.create(newdir)
dir.create(newdir&'/SND')
#===============================================================================
#===============================================================================
#===============================================================================
# PLOT DAILY SND @ POIs
#===============================================================================
#===============================================================================
#===============================================================================

#### LOOP OVER POIs #####
for(i in 1:n){
  ##extract names & tables
  nm<-as.character(pois$name[i])
  id<-as.character(pois$ID[i])
  ele<-as.character(pois$elevation[i])
  
  print('Plotting daily SND-validation at <<'&nm&' ('&ele&'m)>>  ('&i&'/'&n&')')
  
  ##extract station data
  obs<-obsdat[which(obsdat$station_code == id),c('HS','Date')]
  
  ##extract modelled data
  sim<-simdat[,which(names(simdat) %in% c('Date',id))]
  
  ifelse(nrow(obsdat) > 0,
         df<-merge(obs,sim,by='Date',all.y=T),
         df<-cbind.data.frame(sim$Date,NA,sim[,2])
  )
  names(df)<-c('DATE','OBS','SIM')
  
  ##############################
  ### PLOT SND COMPARISON ###
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
  if(is.nan(avg_sim)){avg_sim<-round(mean(df$SIM,na.rm=T),2)}
  
  df<-melt(df,'DATE')
  tit<-'Snow depth '&nm&' \n ('&ele&' m a.s.l.)'
  capt<-''
  xlab<-''
  ylab<-'[cm]'
  capt<-paste('R2 = ',R2,'   |   RMSE = ',RMSE,
              'cm   |   NSE = ',NSE,'\n',
              'Average [cm]: ',avg_obs,' (obs.) | ',avg_sim,' (sim.)')
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
  
  
  fn<-newdir&'/SND/SND-daily_'&nm&'_'&PERIOD_SH&'.png'
  ggsave(fn,pl,width=20,height=15,units='cm')
  graphics.off()
  rm(df,tit,capt,avg_obs,avg_sim,xlab,ylab,leg,lims,pl,fn)
  rm(nm,ele,RMSE,NSE,R2)
}



#===============================================================================
#===============================================================================
#===============================================================================
# PLOT WEEKLY SND @ POIs
#===============================================================================
#===============================================================================
#===============================================================================
#### LOOP OVER POIs #####
for(i in 1:n){
  ##extract names & tables
  nm<-as.character(pois$name[i])
  id<-as.character(pois$ID[i])
  ele<-as.character(pois$elevation[i])
  
  print('Plotting weekly SND-validation at <<'&nm&' ('&ele&'m)>>  ('&i&'/'&n&')')
  
  ##extract station data
  obs<-obsdat[which(obsdat$station_code == id),c('HS','Date')]
  
  ##extract modelled data
  sim<-simdat[,which(names(simdat) %in% c('Date',id))]
  
  ifelse(nrow(obsdat) > 0,
         df<-merge(obs,sim,by='Date',all.y=T),
         df<-cbind.data.frame(sim$Date,NA,sim[,2])
  )
  names(df)<-c('DATE','OBS','SIM')
  
  df<-aggregate(df[,-1],list(DATE=cut(df$DATE,'week')),median) #derive weekly median
  df$DATE<-as.POSIXct(df$DATE,format='%Y-%m-%d')
  
  ##############################
  ### PLOT SND COMPARISON ###
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
  if(is.nan(avg_sim)){avg_sim<-round(mean(df$SIM,na.rm=T),2)}
  
  df<-melt(df,'DATE')
  tit<-'Snow depth '&nm&' \n ('&ele&' m a.s.l.)'
  capt<-''
  xlab<-''
  ylab<-'[cm]'
  capt<-paste('R2 = ',R2,'   |   RMSE = ',RMSE,
              'cm   |   NSE = ',NSE,'\n',
              'Average [cm]: ',avg_obs,' (obs.) | ',avg_sim,' (sim.)')
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
  
  
  fn<-newdir&'/SND/SND-weekly_'&nm&'_'&PERIOD_SH&'.png'
  ggsave(fn,pl,width=20,height=15,units='cm')
  graphics.off()
  rm(df,tit,capt,avg_obs,avg_sim,xlab,ylab,leg,lims,pl,fn)
  rm(nm,ele,RMSE,NSE,R2)
}



#===============================================================================
#===============================================================================
#===============================================================================
# PLOT MONTHLY SND @ POIs
#===============================================================================
#===============================================================================
#===============================================================================
#### LOOP OVER POIs #####
for(i in 1:n){
  ##extract names & tables
  nm<-as.character(pois$name[i])
  id<-as.character(pois$ID[i])
  ele<-as.character(pois$elevation[i])
  
  print('Plotting monthly SND-validation at <<'&nm&' ('&ele&'m)>>  ('&i&'/'&n&')')
  
  ##extract station data
  obs<-obsdat[which(obsdat$station_code == id),c('HS','Date')]
  
  ##extract modelled data
  sim<-simdat[,which(names(simdat) %in% c('Date',id))]
  
  ifelse(nrow(obsdat) > 0,
         df<-merge(obs,sim,by='Date',all.y=T),
         df<-cbind.data.frame(sim$Date,NA,sim[,2])
  )
  names(df)<-c('DATE','OBS','SIM')
  
  df<-aggregate(df[,-1],list(DATE=cut(df$DATE,'month')),median) #derive weekly median
  df$DATE<-as.POSIXct(df$DATE,format='%Y-%m-%d')
  
  ##############################
  ### PLOT SND COMPARISON ###
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
  if(is.nan(avg_sim)){avg_sim<-round(mean(df$SIM,na.rm=T),2)}
  
  df<-melt(df,'DATE')
  tit<-'Snow depth '&nm&' \n ('&ele&' m a.s.l.)'
  capt<-''
  xlab<-''
  ylab<-'[cm]'
  capt<-paste('R2 = ',R2,'   |   RMSE = ',RMSE,
              'cm   |   NSE = ',NSE,'\n',
              'Average [cm]: ',avg_obs,' (obs.) | ',avg_sim,' (sim.)')
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
  
  
  fn<-newdir&'/SND/SND-monthly_'&nm&'_'&PERIOD_SH&'.png'
  ggsave(fn,pl,width=20,height=15,units='cm')
  graphics.off()
  rm(df,tit,capt,avg_obs,avg_sim,xlab,ylab,leg,lims,pl,fn)
  rm(nm,ele,RMSE,NSE,R2)
}


print(outnm)

