################################################################################
# T&C OUTPUT :
#  COMPARE METEO | PLOT
#
# 
# TODO: -
#       -
#
# NEW:  -
#       -
#
# 2023/11/01
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

## define folder of reference runoff data
refMeteofo<-root&'/Meteodata/CH/SwissMetNet_HourlyData_fromMZ'


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
base_size<-18
base_family<-'Arial'


#===============================================================================
# TIME PERIOD OF INTEREST
#===============================================================================
ts_init<-as.POSIXct(ts_init,format='%Y-%m-%d %H:%M:%S')
ts_final<-as.POSIXct(ts_final,format='%Y-%m-%d %H:%M:%S')


#===============================================================================
# FIND SIMULATED METEO-STATIONS WITHIN DOMAIN
#===============================================================================
pois<-st_read(path_data&'/SHAPEFILES/POIs.shp',quiet=T)
sim_ids<-which(pois$ID2 %in% 'MET')
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
# READ MODELLED METEO DATA
#===============================================================================
setwd(path_data&'/TABLES')
n<-nrow(pois)

simdat<-list()
for(i in 1:n){
  STA_NM<-pois$ID[i]&'-'&pois$ID2[i]
  ELEV<-as.character(pois$elevation[i])
  
  print(DOMAIN&': reading modelled Meteo-data at station <<'&STA_NM&' ('&ELEV&'m)>>  ('&i&'/'&n&')')
  #read data
  fn<-'TrackedPixel___'&STA_NM&'___'&ELEV&'m.txt'
  df<-read.table(fn,header=T,sep='\t',dec='.',na.strings='NaN') 
  ##correct for midnight hour
  df$Date<-as.character(df$Date)
  idx<-which(nchar(df$Date) < 19)
  df$Date[idx]<-df$Date[idx]&' 00:00:00'
  rm(idx)
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S',tz=TZ)
  df<-df[,which(names(df) %in% c('Date','ea_S','Ta_S','Ws_S','Pr_S','Rsw_space'))]
  #derive daily avg.
  df_avg<-aggregate(df[,names(df) %in% c('ea_S','Ta_S','Ws_S','Rsw_space')],list(Date=cut(df$Date,'day')),mean) 
  df_avg$Date<-as.POSIXct(df_avg$Date,format='%Y-%m-%d')
  #derive daily sum
  df_sum<-aggregate(df[,'Pr_S'],list(Date=cut(df$Date,'day')),sum)
  df<-cbind.data.frame(df_avg,Pr_S=df_sum[,-1])
  #crop to selected period
  df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] 
  
  colnames(df)<-str_replace(colnames(df),'_S','')
  colnames(df)<-str_replace(colnames(df),'_space','')
  simdat[[i]]<-df
  rm(df)
}
names(simdat)<-pois$name
# plot(simdat[[1]]$Ta,type='l') ##check


#===============================================================================
# RE-ADJUST TIMESPAN TO ANALYZE (IN CASE SIMULATION PERIOD IS SHORTER)
#===============================================================================
if(simdat[[1]]$Date[1] > ts_init){ts_init<-simdat[[1]]$Date[1]}
if(simdat[[1]]$Date[nrow(simdat[[1]])] < ts_final){ts_final<-simdat[[1]]$Date[nrow(simdat[[1]])]}
PERIOD_SH<-format(ts_init,"%b%Y")&'-'&format(ts_final-86400,"%b%Y")
PERIOD_LO<-format(ts_init,"%d %b %Y")&' - '&format(ts_final-86400,"%d %b %Y")


#===============================================================================
# READ OBSERVED SWISSMETNET-DATA
#===============================================================================
setwd(refMeteofo)
yrs<-year(ts_init):year(ts_final)

vars<-c('Ddruck',   ##vapor pressure [bar]
        'global',   ##incoming SW radiation [W/m2]
        'temp',     ##air temperature [°C]
        'wind',     ##wind speed [m/s]
        'precip')   ##precipitation [mm]
vars2<-c('ea','Rsw','Ta','Ws','Pr')

obsdat<-list()
for(i in 1:n){ ##loop over POIs
  POI<-names(simdat)[i]
  print(DOMAIN&': reading observed Meteo-data at station <<'&POI&'>>  ('&i&'/'&n&')')
  sim_loc<-cbind.data.frame(x=mroundd(pois$chx[i],250),
                            y=mroundd(pois$chy[i],250),
                            # z=mroundd(pois$elevation[i],100),
                            nm=POI)
  
  dat2_ls<-list()
  for(j in 1:length(vars)){ ##loop over variables
    VAR<-vars[j]
    VAR2<-vars2[j]
    print('++++++++++++++++++++++ '&VAR2)
    
    dat1_ls<-list()
    for(k in 1:length(yrs)){ ##loop over years
      YR<-as.character(yrs[k])
      
      dat<-read.table(YR&'_stundenmittel_'&VAR&'.dat',header=F,skip=1)
      colnames(dat)<-dat[4,]
      #####
      obs_loc<-cbind.data.frame(x=mroundd(as.numeric(dat[2,-c(1:4)])+2*10^6,250),
                                y=mroundd(as.numeric(dat[3,-c(1:4)])+1*10^6,250),
                                # z=mroundd(as.numeric(dat[1,-c(1:4)]),100),
                                nm=colnames(dat)[-c(1:4)])
      #####
      idx<-which(obs_loc$x == sim_loc$x & obs_loc$y == sim_loc$y)
      dat<-dat[-c(1:4),c(1:4,4+idx)]
      dat<-as.data.frame(sapply(dat[,],as.numeric))
      if(length(idx) == 0){dat$na<-NA}
      dat$Date<-as.POSIXct(paste(dat$YY,dat$MM,dat$DD,dat$HH),format='%Y %m %d %H')
      dat<-dat[,-c(1:4)]
      names(dat)<-c(VAR2,'Date')
      dat1_ls[[k]]<-dat
      rm(YR,obs_loc,idx,dat)
    }
    dat1<-as.data.frame(do.call('rbind',dat1_ls))
    dat2_ls[[j]]<-dat1
    rm(VAR,VAR2,dat1,dat1_ls)
  }
  dat2<-as.data.frame(do.call('cbind',dat2_ls))
  idx<-which(names(dat2) == 'Date')
  dat2<-cbind.data.frame(Date=dat2[,idx[1]],dat2[,-idx])

  #derive daily avg.
  df_avg<-aggregate(dat2[,names(dat2) %in% c('ea','Ta','Ws','Rsw')],list(Date=cut(dat2$Date,'day')),mean) 
  df_avg$Date<-as.POSIXct(df_avg$Date,format='%Y-%m-%d')
  #derive daily sum
  df_sum<-aggregate(dat2[,'Pr'],list(Date=cut(dat2$Date,'day')),sum)
  dat2<-cbind.data.frame(df_avg,Pr=df_sum[,-1])
  #crop to selected period
  dat2<-dat2[which(dat2[,1] >= ts_init & dat2[,1] <= ts_final),] 
  #define NAs
  dat2[dat2 == -9999]<-NA
  obsdat[[i]]<-dat2
  rm(POI,idx,df_avg,df_sum,dat2,dat2_ls)
}
names(obsdat)<-pois$name
# plot(obsdat[[1]]$Ta,type='l') ##check

newdir<-path_data&'/PLOTS'
dir.create(newdir)
dir.create(newdir&'/METEO')
#===============================================================================
#===============================================================================
#===============================================================================
# PLOT DAILY METEO VALIDATION @ POIs
#===============================================================================
#===============================================================================
#===============================================================================
#### LOOP OVER POIs #####
for(i in 1:n){
  ##extract names & tables
  nm<-as.character(pois$name[i])
  id<-as.character(pois$ID[i])
  ele<-as.character(pois$elevation[i])
  
  print(DOMAIN&': plotting daily Meteo-validation at <<'&nm&' ('&ele&'m)>>  ('&i&'/'&n&')')
  
  obs<-obsdat[[i]]
  sim<-simdat[[i]]
  df<-merge(obs,sim,by='Date',all.y=T)
  names(df)<-str_replace(names(df),'.x','_obs')
  names(df)<-str_replace(names(df),'.y','_sim')
  
  pl_ls<-list()
  for(j in 1:length(vars2)){    ###"aes_string()" needed for ggplot in loop!
    VAR<-vars2[j]
    DF<-cbind.data.frame(Date=df$Date,
                         OBS=df[,which(names(df) == paste0(VAR,'_obs'))],
                         SIM=df[,which(names(df) == paste0(VAR&'_sim'))])
    
    ##QC
    

    if(VAR == 'ea'){tit<-'Vapour pressure'; un<-'Pa'; ylab<-'Daily mean ['&un&']'; DF[DF < 0]<-NA; DF$OBS<-DF$OBS*100}
    if(VAR == 'Rsw'){tit<-'Incoming shortwave radiation'; un<-'W/m2'; ylab<-'Daily mean ['&un&']'; DF[DF < 0]<-NA}
    if(VAR == 'Ta'){tit<-'Air temperature (2m)'; un<-'°C'; ylab<-'Daily mean ['&un&']'; DF[DF < -100]<-NA}
    if(VAR == 'Ws'){tit<-'Wind speed'; un<-'m/s'; ylab<-'Daily mean ['&un&']'; DF[DF < 0]<-NA}
    if(VAR == 'Pr'){tit<-'Total precipitation'; un<-'mm'; ylab<-'Daily sum ['&un&']'; DF[DF < 0]<-NA}
    
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
    if(VAR == 'Pr'){avg_obs<-round(avg_obs*365.25,0); avg_sim<-round(avg_sim*365.25,0)}
    
    capt<-paste0('R2 = ',R2,'   |   RMSE = ',RMSE,' ',
                un,'   |   NSE = ',NSE,
                '   |   Average [',un,']: ',avg_obs,' (obs.) | ',avg_sim,' (sim.)')
    if(VAR == 'Pr'){capt<-str_replace(capt,'Average','Avg. annual sum')}

    DF<-melt(DF,'Date')

    pl<-ggplot(DF)+
      geom_line(aes_string(x=as.Date(DF[,'Date']),y=DF[,'value'],colour=DF[,'variable']),linewidth=0.5)+
      labs(title=tit,subtitle='',caption=capt,x='',y=ylab,fill='')+
      scale_x_date(date_labels="%b")+
      scale_color_manual(values=c('black','darkorange'))+
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
  
  fn<-newdir&'/METEO/validation-daily___'&id&'___ea-Rsw-Ta-Ws-Pr_'&PERIOD_SH&'.png'
  ggsave(fn,pl,width=30,height=40,units='cm') ##n=18
  graphics.off()
  
}


print('+++END SCRIPT <<VALIDATION SWISSMETNET POIS>>+++   '&outnm)



