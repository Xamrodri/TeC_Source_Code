################################################################################
# T&C OUTPUT :
#  COMPARE SNOW DEPTH/TAIR/WS/SWOUT | PLOT
#
# 
# TODO: -
#       -
#
# NEW:  -
#       -
#
# 2023/11/03
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
# outnm<-'v20231006_Aaretal_404days'
# outnm<-'v20231006_Aargauer_Rhein_2557days'
# outnm<-'v20231006_Berner_Oberland_779days'
# outnm<-'v20231006_Bodensee_2557days'
outnm<-'v20231006_Engadin_1906days'
# outnm<-'v20231006_Glatt_2557days'
# outnm<-'v20231006_Hinterrhein_1853days'
# outnm<-'v20231006_Jura_1562days'
# outnm<-'v20231006_Leman_932days'
# outnm<-'v20231006_Limmat_841days'
# outnm<-'v20231006_Neuchatel_813days'
# outnm<-'v20231006_Oberwallis_1498days'
# outnm<-'v20231006_Reuss_1257days'
# outnm<-'v20231006_Rheintal_402days'
# outnm<-'v20231006_Saane_566days'
# outnm<-'v20231006_Schaffhausen_2557days'
# outnm<-'v20231006_Ticino_130days'
# outnm<-'v20231006_Thur_841days'
# outnm<-'v20231006_Toess_1130days'
# outnm<-'v20231006_Unterwallis_376days'
# outnm<-'v20231006_Vorderrhein_1123days'
# outnm<-'v20231006_Zuercher_Rhein_2557days'
###
# outnm<-'v20231006_Test-Aletsch_2557days'
# outnm<-'v20231006_Test-Grossbach_2557days'
# outnm<-'v20231006_Test-Poschiavino_2557days'
# outnm<-'v20231006_Test-Rein_da_Sumvitg_2557days'
# outnm<-'v20231006_Test-Rietholzbach_2557days'
# outnm<-'v20231006_Test-Roggiasca_2557days'
# outnm<-'v20231006_Test-Rosegbach_2557days'


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

## define folder of reference IMIS data
refIMISfo<-root&'/Meteodata/CH/IMIS'


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
library(gridExtra)
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


#===============================================================================
# FIND SIMULATED IMIS-STATIONS WITHIN DOMAIN
#===============================================================================
pois<-st_read(path_data&'/SHAPEFILES/POIs.shp',quiet=T)
sim_ids<-which(pois$ID2 %in% 'IMI')
if(length(sim_ids) == 0 && nrow(pois) == 1){sim_ids<-1}  ##Test-domains
pois<-pois[sim_ids,]


#===============================================================================
# EXTRACT DOMAIN
#===============================================================================
DOMAIN<-str_split(outnm,'_')[[1]]
DOMAIN<-DOMAIN[-1]
DOMAIN<-DOMAIN[-length(DOMAIN)]
if(length(DOMAIN) > 1){DOMAIN<-paste(DOMAIN,collapse=' ')}


#===============================================================================
# READ MODELLED SNOW DATA
#===============================================================================
setwd(path_data&'/TABLES')
n<-nrow(pois)

simdat<-list()
for(i in 1:n){
  STA_NM<-pois$ID[i]&'-'&pois$ID2[i]
  ELEV<-as.character(pois$elevation[i])
  
  print(DOMAIN&': reading modelled snow-data at <<'&STA_NM&' ('&ELEV&'m)>>  ('&i&'/'&n&')')
  #read data
  fn<-'TrackedPixel___'&STA_NM&'___'&ELEV&'m.txt'
  df<-read.table(fn,header=T,sep='\t',dec='.',na.strings='NaN') 
  ##correct for midnight hour
  df$Date<-as.character(df$Date)
  idx<-which(nchar(df$Date) < 19)
  df$Date[idx]<-df$Date[idx]&' 00:00:00'
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S',tz=TZ)
  
  df<-df[,which(names(df) %in% c('Date','Ta_S','Ws_S','SND','snow_albedo','Rsw_space'))]
  df[df$snow_albedo < 0.1,'snow_albedo']<-NA #discard night hours (albedo =0)
  df$SWout<-df$Rsw_space*df$snow_albedo
  df$SND<-df$SND*100 ##m to cm
  
  df<-aggregate(df[,-1],list(Date=cut(df$Date,'day')),mean,na.rm=T) #derive daily mean
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d')
  
  df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period
  simdat[[i]]<-cbind.data.frame(Date=df$Date,Ta=df$Ta_S,Ws=df$Ws_S,SND=df$SND,alb=df$snow_albedo,SWout=df$SWout)
  rm(STA_NM,ELEV,fn,df,idx)
}
names(simdat)<-pois$ID


#===============================================================================
# RE-ADJUST TIMESPAN TO ANALYZE (IN CASE SIMULATION PERIOD IS SHORTER)
#===============================================================================
if(simdat[[1]]$Date[1] > ts_init){ts_init<-simdat[[1]]$Date[1]}
if(simdat[[1]]$Date[nrow(simdat[[1]])] < ts_final){ts_final<-simdat[[1]]$Date[nrow(simdat[[1]])]}
PERIOD_SH<-format(ts_init,"%b%Y")&'-'&format(ts_final-86400,"%b%Y")
PERIOD_LO<-format(ts_init,"%d %b %Y")&' - '&format(ts_final-86400,"%d %b %Y")


#===============================================================================
# READ OBSERVED IMIS SNOW/METEO DATA
#===============================================================================
setwd(refIMISfo&'/Halfhourly_data')
files<-list.files(pattern='.csv')
obsdat<-list()
for(i in 1:n){
  ID<-pois$ID[i]
  ELEV<-pois$elevation[i]
  print(DOMAIN&': reading observed IMIS-meteodata at <<'&ID&' ('&ELEV&'m)>>  ('&i&'/'&n&')')
  fn<-ID&'.csv'
  
  if(fn %in% files){
    df<-read.table(fn,sep=',',header=T)
    if(length(df$TA_30MIN_MEAN) == 0){df$TA_30MIN_MEAN<-NA}
    if(length(df$VW_30MIN_MEAN) == 0){df$VW_30MIN_MEAN<-NA}
    if(length(df$HS) == 0){df$HS<-NA}
    if(length(df$RSWR_30MIN_MEAN) == 0){df$RSWR_30MIN_MEAN<-NA}
    df<-cbind.data.frame(Date=as.POSIXct(df$measure_date,format='%Y-%m-%d %H:%M:%S',tz=TZ),
                         Ta=df$TA_30MIN_MEAN,
                         Ws=df$VW_30MIN_MEAN,
                         SND=df$HS,
                         SWout=df$RSWR_30MIN_MEAN
    )
    
    df<-aggregate(df[,-1],list(Date=cut(df$Date,'day')),mean,na.rm=T) #derive daily mean
    df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d')
    
    df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period
    ##QC
    df[which(df$Ta < -100),'Ta']<-NA
    df[which(df$Ws < 0),'Ws']<-NA
    df[which(df$SND < 0),'SND']<-NA
    df[which(df$SWout < 0),'SWout']<-NA
  }
  
  if(fn %notin% files){df<-NA}
    
  obsdat[[i]]<-df
  rm(ID,ELEV,fn,df)
}


newdir<-path_data&'/PLOTS'
dir.create(newdir)
dir.create(newdir&'/SNOW')
#===============================================================================
#===============================================================================
#===============================================================================
# PLOT DAILY SNOW VALIDATION @ POIs
#===============================================================================
#===============================================================================
#===============================================================================
#### LOOP OVER POIs #####
for(i in 1:n){
  ##extract names & tables
  nm<-as.character(pois$name[i])
  id<-as.character(pois$ID[i])
  ele<-as.character(pois$elevation[i])
  
  print(DOMAIN&': plotting daily Snow-validation at <<'&nm&' ('&ele&'m)>>  ('&i&'/'&n&')')
  
  obs<-obsdat[[i]]
  sim<-simdat[[i]]
  vars<-names(sim)[-1]
    
  if(is.null(nrow(obs))){
    obs<-sim
    obs[,-1]<-obs[,-1]*0
    obs[obs == 0]<-NA
    obs<-obs[,-which(names(obs) == 'alb')]
  }
  df<-merge(obs,sim,by='Date',all.y=T)
  names(df)<-str_replace(names(df),'.x','_obs')
  names(df)<-str_replace(names(df),'.y','_sim')
  

  pl_ls<-list()
  for(j in 1:length(vars)){    ###"aes_string()" needed for ggplot in loop!
    VAR<-vars[j]
    if(VAR != 'alb'){
      DF<-cbind.data.frame(Date=df$Date,
                           OBS=df[,which(names(df) == paste0(VAR,'_obs'))],
                           SIM=df[,which(names(df) == paste0(VAR&'_sim'))])
    }
    if(VAR == 'alb'){
      DF<-cbind.data.frame(Date=df$Date,
                           OBS=NA,
                           SIM=df[,which(names(df) == VAR)])
    }
    
    ##QC
    if(VAR == 'Ta'){tit<-'Air temperature'; un<-'°C'; ylab<-'Daily mean ['&un&']'; DF[DF < -100]<-NA}
    if(VAR == 'Ws'){tit<-'Wind speed'; un<-'m/s'; ylab<-'Daily mean ['&un&']'; DF[DF <= 0]<-NA}
    if(VAR == 'SND'){tit<-'Snow depth'; un<-'cm'; ylab<-'Daily mean ['&un&']'; DF[DF < 0]<-NA}
    if(VAR == 'alb'){tit<-'Albedo'; un<-'-'; ylab<-'Daily mean ['&un&']'; DF[DF < 0]<-NA}
    if(VAR == 'SWout'){tit<-'Reflected shortwave radiation'; un<-'W/m2'; ylab<-'Daily mean ['&un&']'; DF[DF < 0]<-NA}
    
    ## PERFORMANCE ##
    ##calculate RMSE
    comp<-DF[complete.cases(DF),] 
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
    if(is.nan(avg_sim)){avg_sim<-round(mean(DF$SIM,na.rm=T),2)}
    capt<-paste0('R2 = ',R2,'   |   RMSE = ',RMSE,' ',
                 un,'   |   NSE = ',NSE,
                 '   |   Average [',un,']: ',avg_obs,' (obs.) | ',avg_sim,' (sim.)')
    
    DF<-melt(DF,'Date')
    
    pl<-ggplot(DF)+
      geom_line(aes_string(x=as.Date(DF[,'Date']),y=DF[,'value'],colour=DF[,'variable']),linewidth=0.5)+
      labs(title=tit,subtitle='',caption=capt,x='',y=ylab,fill='')+
      scale_x_date(date_labels="%b")+
      scale_color_manual(values=c('black','steelblue3'))+
      theme_gray(base_size=base_size,base_family=base_family)+
      theme(
        legend.position = 'right',
        legend.direction = 'vertical',
        # legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size=base_size),
        legend.key=element_rect(fill='white'),
        axis.ticks = element_blank(),
        axis.text.y = element_text(colour='black',size=base_size),
        axis.text.x = element_text(colour='black',size=base_size),
        axis.title= element_text(colour='black',size=base_size),
        panel.background = element_rect(fill='white',colour='white'),
        panel.grid.minor.y = element_line(linewidth=0.1,colour='grey'),
        panel.grid.major.y = element_line(linewidth=0.75,colour='grey'),
        panel.grid.major.x = element_line(linewidth=0.75,colour='grey'),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(colour='black',size=base_size),
        plot.subtitle = element_text(colour='black',size=base_size),
        plot.margin=unit(c(1,1,0,0),'lines') #top/right/bottom/left
      )
    print(pl)
    Sys.sleep(1)
    pl_ls[[j]]<-pl
    rm(pl)
  }
  
  toptit<-nm&' ('&id&'),  '&ele&' m a.s.l.  ('&PERIOD_LO&')'
  # pl<-marrangeGrob(pl_ls,nrow=5,ncol=1,top=toptit) ##n=5
  pl<-marrangeGrob(pl_ls,nrow=5,ncol=1,top=grid::textGrob(toptit,gp=grid::gpar(fontsize=base_size))) ##n=5
  
  pl
  
  fn<-newdir&'/SNOW/validation-daily___'&id&'___ea-Rsw-Ta-Ws-Pr_'&PERIOD_SH&'.png'
  ggsave(fn,pl,width=30,height=40,units='cm') ##n=18
  graphics.off()
  
}




print('+++END SCRIPT <<VALIDATION IMIS POIS>>+++   '&outnm)



