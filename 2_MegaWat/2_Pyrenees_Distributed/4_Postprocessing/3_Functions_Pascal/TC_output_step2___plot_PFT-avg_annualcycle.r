################################################################################
# T&C OUTPUT :
#  PLOT SPATAL AVERAGE TIMESERIES PER PFT OF GIVEN DOMAIN
#
# 
# TODO: -
#       -
#
# NEW:  -
#       -
#
# 2023/10/05
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
outnm<-'v20231006_Glatt_2557days'
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
# outnm<-'v20231006_Toess_2557days'
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
ts_init<-'2016-01-01 02:00:00'
ts_final<-'2022-09-30 23:00:00'
# ts_init<-'2015-10-01 02:00:00'
# ts_final<-'2016-08-01 23:00:00'

# ## define drought events (for plotting)
# dr1_init<-'2018-01-01 01:00:00'
# dr1_final<-'2018-12-31 23:00:00'
# ##
# dr2_init<-'2022-01-01 01:00:00'
# dr2_final<-'2022-09-30 23:00:00'


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
library(lubridate)    #floor_date(), hour() 
library(networkD3)    #Sankey diagram
library(dplyr)        #Sankey diagram
library(webshot)      #save HTML screenshot
library(Metrics)      #rmse()
library(topmodel)     #NSeff()
library(stringr)      #str_replace
library(oce)          #despike()
library(zoom)         #zm()
library(extrafont)    #loadfonts()
library(lubridate)
library(assertthat)
library(grid)        #textGrob()
library(gridExtra)   #marrangeGrob()
library(MetBrewer)   #met.brewer()   
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

# dr1_init<-as.POSIXct(dr1_init,format='%Y-%m-%d %H:%M:%S')
# dr1_final<-as.POSIXct(dr1_final,format='%Y-%m-%d %H:%M:%S')

PERIOD_SH<-format(ts_init,"%b%Y")&'-'&format(ts_final-86400,"%b%Y")
PERIOD_LO<-format(ts_init,"%d %b %Y")&' - '&format(ts_final-86400,"%d %b %Y")


#===============================================================================
# EXTRACT DOMAIN NAME
#===============================================================================
DOMAIN<-str_split(outnm,'_')[[1]]
DOMAIN<-DOMAIN[-c(1,length(DOMAIN))]


setwd(path_data&'/RASTERS/INIT-COND')
#===============================================================================
# READ LC, OSAT & OHY
#===============================================================================
lc<-raster('VEG_CODE.tif')
ohy<-raster('ohy.tif')
osat<-raster('osat.tif')


#===============================================================================
# READ MODELLED PFT-AVERAGE PER DOMAIN DATA
#===============================================================================
setwd(path_data&'/TABLES')
fn<-list.files(pattern='SPAVG-PFT')

nm_v<-vector()
df_ls<-list()
dfweekly_ls<-list()
for(i in 1:length(fn)){
  nm<-str_split(fn[i],'SPAVG-')[[1]][2]
  nm_v[i]<-str_split(nm,'.txt')[[1]][1]
  rm(nm)
  
  df<-read.table(fn[i],header=T,sep='\t',dec='.')
  ##plot(df$AgeL_H[1:17000],type='l')
  ##plot(df$LAIdead_H [1:17000],type='l')
  ##plot(df$LAI_H [1:17000],type='l')
  ##plot(df$PHE_S_H [1:17000],type='l')
  
  ##correct for midnight hour
  df$Date<-as.character(df$Date)
  idx<-which(nchar(df$Date) < 19)
  df$Date[idx]<-df$Date[idx]&' 00:00:00'
  rm(idx)
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S',tz=TZ)
  
  #derive daily mean
  # sumidx<-which(names(df) %in% c('Pr','Pr_liq','Pr_sno','T'))
  df_avg<-aggregate(df[,-1],list(Date=cut(df$Date,'day')),mean) 
  df_avg$Date<-as.POSIXct(df_avg$Date,format='%Y-%m-%d')
  # #derive daily sum
  # df_sum<-aggregate(df[,sumidx],list(Date=cut(df$Date,'day')),sum) 
  # df<-cbind.data.frame(df_avg,df_sum[,-1])
  df<-df_avg
  df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period
  df_ls[[i]]<-df
  rm(df_avg)
  
  #derive weekly mean
  # sumidx<-which(names(df) %in% c('Pr','Pr_liq','Pr_sno','T'))
  df_avg<-aggregate(df[,-1],list(Date=cut(df$Date,'week')),mean) 
  df_avg$Date<-as.POSIXct(df_avg$Date,format='%Y-%m-%d')
  df<-df_avg
  df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period
  dfweekly_ls[[i]]<-df
  rm(df_avg,df)
}
rm(fn)
names(df_ls)<-nm_v
names(dfweekly_ls)<-nm_v
PFTs<-str_replace(nm_v,'PFT1___',''); PFTs<-str_replace(PFTs,'PFT2___','')
PFTs<-unique(PFTs)


#===============================================================================
# READ MODELLED LC-AVERAGE PER DOMAIN DATA
#===============================================================================
df2_ls<-list()
for(i in 1:length(PFTs)){
  fn<-'SPAVG-LC___'&PFTs[i]&'.txt'
  print(DOMAIN&'   +++ reading <<'&fn&'>> +++')
  df<-read.table(fn,header=T,sep='\t',dec='.')
  ##plot(df$O[1:17000],type='l')
  ##plot(df$V [1:17000],type='l')
  
  ##correct for midnight hour
  df$Date<-as.character(df$Date)
  idx<-which(nchar(df$Date) < 19)
  df$Date[idx]<-df$Date[idx]&' 00:00:00'
  rm(idx)
  df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S',tz=TZ)
  
  #derive daily mean
  # sumidx<-which(names(df) %in% c('Pr','Pr_liq','Pr_sno','T'))
  df_avg<-aggregate(df[,-1],list(Date=cut(df$Date,'day')),mean) 
  df_avg$Date<-as.POSIXct(df_avg$Date,format='%Y-%m-%d')
  # #derive daily sum
  # df_sum<-aggregate(df[,sumidx],list(Date=cut(df$Date,'day')),sum) 
  # df<-cbind.data.frame(df_avg,df_sum[,-1])
  df<-df_avg
  df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period
  
  ##############################################
  # DERIVE EFFECTIVE SOIL SATURATION (PER PFT)
  #
  # Se=(O-Ohy)/(Osat-Ohy)
  #   Osat = saturated water content [-]   
  #   Ohy = residual or hygroscopic moisture content [-] 
  #
  # %%% #1) Evergreen Tree Forest
  # %%% #2) Mixed Evergreen (50%) Decidous (50%)
  # %%% #3) Decidous
  # %%% #4) Meadow - Grasses
  # %%% #5) Rock
  # %%% #6) Water
  # %%% #7) Ice (bedrock)
  if(PFTs[i] == 'Evergreen'){lccode<-1}
  if(PFTs[i] == 'Mixed-Evergreen-Decidous'){lccode<-2}
  if(PFTs[i] == 'Decidous'){lccode<-3}
  if(PFTs[i] == 'Grass'){lccode<-4}
  OHY<-ohy;    OHY[lc != lccode]<-NA;   OHY<-mean(values(OHY),na.rm=T)
  OSAT<-osat;  OSAT[lc != lccode]<-NA;  OSAT<-mean(values(OSAT),na.rm=T)
  
  df$Se<-(df$O-OHY)/(OSAT-OHY)
  ##############################################
  
  df2_ls[[i]]<-df
  rm(df_avg,df)
}
names(df2_ls)<-PFTs


newdir<-path_data&'/PLOTS'
dir.create(newdir)
dir.create(newdir&'/VEG')
#===============================================================================
# +++ PLOT LAI AVERAGE PER DOY (1 LINE PER YEAR) +++ 
# PANEL 1: DECIDOUS
# PANEL 2: EVERGREEN
# PANEL 3: GRASS
# PANEL 4: MIXED EVERGREEN-DECIDOUS
#===============================================================================
Date<-df_ls[[1]][,'Date']
DF1<-lapply(df_ls,'[',,'LAI_H'); DF1<-as.data.frame(do.call('cbind',DF1))
DF2<-lapply(df_ls,'[',,'LAI_L'); DF2<-as.data.frame(do.call('cbind',DF2))
DF<-DF1+DF2
DF<-cbind.data.frame(Date,DF)
idx<-which(grepl('mix',names(DF),ignore.case=T))
DF$mixed<-rowSums(DF[,idx])
names(DF)[ncol(DF)]<-names(DF)[idx[1]]
DF<-DF[,-idx]
names(DF)<-str_replace(names(DF),'PFT1___','')
rng<-c(0,mroundu(max(DF[,-1]),0.5))
ylab<-expression(paste('LAI  ','[',m^2,m^-2,']'))
xlab<-'DOY'
lwd<-1.5

DF$yr<-year(DF$Date)
DF$doy<-yday(DF$Date)
rm(Date,DF1,DF2,idx)

### PANEL 1: DECIDOUS ###
PFT<-'Decidous'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl1<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 2: EVERGREEN ###
PFT<-'Evergreen'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl2<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 3: GRASS ###
PFT<-'Grass'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl3<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 4: MIXED EVERGREEN-DECIDOUS ###
PFT<-'Mixed-Evergreen-Decidous'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl4<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)


toptit<-textGrob(DOMAIN,
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3,pl4),nrow=4,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/VEG/daily_PFT-avg_LAI_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=30,units='cm')
graphics.off()
rm(fn,DF,pl)


#===============================================================================
# +++ PLOT LAI AVERAGE PER DOY (1 LINE PER YEAR) +++ 
# PANEL 1: DECIDOUS
# PANEL 2: EVERGREEN
# PANEL 3: GRASS
# PANEL 4: MIXED EVERGREEN-DECIDOUS
#===============================================================================
Date<-df_ls[[1]][,'Date']
DF1<-lapply(df_ls,'[',,'LAIdead_H'); DF1<-as.data.frame(do.call('cbind',DF1))
DF2<-lapply(df_ls,'[',,'LAIdead_L'); DF2<-as.data.frame(do.call('cbind',DF2))
DF<-DF1+DF2
DF<-cbind.data.frame(Date,DF)
idx<-which(grepl('mix',names(DF),ignore.case=T))
DF$mixed<-rowSums(DF[,idx])
names(DF)[ncol(DF)]<-names(DF)[idx[1]]
DF<-DF[,-idx]
names(DF)<-str_replace(names(DF),'PFT1___','')
rng<-c(0,mroundu(max(DF[,-1]),0.5))
ylab<-expression(paste('Dead LAI  ','[',m^2,m^-2,']'))
xlab<-'DOY'
lwd<-1.5

DF$yr<-year(DF$Date)
DF$doy<-yday(DF$Date)
rm(Date,DF1,DF2,idx)

### PANEL 1: DECIDOUS ###
PFT<-'Decidous'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl1<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 2: EVERGREEN ###
PFT<-'Evergreen'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl2<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 3: GRASS ###
PFT<-'Grass'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl3<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 4: MIXED EVERGREEN-DECIDOUS ###
PFT<-'Mixed-Evergreen-Decidous'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl4<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)


toptit<-textGrob(DOMAIN,
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3,pl4),nrow=4,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/VEG/daily_PFT-avg_LAIdead_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=30,units='cm')
graphics.off()
rm(fn,DF,pl)


#===============================================================================
# +++ PLOT GPP AVERAGE PER DOY (1 LINE PER YEAR) +++ 
# PANEL 1: DECIDOUS
# PANEL 2: EVERGREEN
# PANEL 3: GRASS
# PANEL 4: MIXED EVERGREEN-DECIDOUS
#===============================================================================
Date<-dfweekly_ls[[1]][,'Date']
DF1a<-lapply(dfweekly_ls,'[',,'NPP_H'); DF1a<-as.data.frame(do.call('cbind',DF1a))
DF2a<-lapply(dfweekly_ls,'[',,'NPP_L'); DF2a<-as.data.frame(do.call('cbind',DF2a))
DF1b<-lapply(dfweekly_ls,'[',,'RA_H'); DF1b<-as.data.frame(do.call('cbind',DF1b))
DF2b<-lapply(dfweekly_ls,'[',,'RA_L'); DF2b<-as.data.frame(do.call('cbind',DF2b))
DF<-DF1a+DF1b+DF1b+DF2b
DF<-cbind.data.frame(Date,DF)
idx<-which(grepl('mix',names(DF),ignore.case=T))
DF$mixed<-rowSums(DF[,idx])
names(DF)[ncol(DF)]<-names(DF)[idx[1]]
DF<-DF[,-idx]
names(DF)<-str_replace(names(DF),'PFT1___','')
rng<-c(0,mroundu(max(DF[,-1]),5))
ylab<-expression(paste('(weekly) GPP ','[gC ',m^-2,day^-1,']'))
xlab<-'DOY'
lwd<-1.5

DF$yr<-year(DF$Date)
DF$doy<-yday(DF$Date)
rm(Date,DF1a,DF2a,DF1b,DF2b,idx)

### PANEL 1: DECIDOUS ###
PFT<-'Decidous'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl1<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 2: EVERGREEN ###
PFT<-'Evergreen'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl2<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 3: GRASS ###
PFT<-'Grass'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl3<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 4: MIXED EVERGREEN-DECIDOUS ###
PFT<-'Mixed-Evergreen-Decidous'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl4<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)


toptit<-textGrob(DOMAIN,
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3,pl4),nrow=4,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/VEG/daily_PFT-avg_GPP_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=30,units='cm')
graphics.off()
rm(fn,DF,pl)



#===============================================================================
# +++ PLOT Se AVERAGE PER DOY (1 LINE PER YEAR) +++ 
# PANEL 1: DECIDOUS
# PANEL 2: EVERGREEN
# PANEL 3: GRASS
# PANEL 4: MIXED EVERGREEN-DECIDOUS
#===============================================================================
Date<-df2_ls[[1]][,'Date']
DF<-lapply(df2_ls,'[',,'Se'); DF<-as.data.frame(do.call('cbind',DF))
DF<-cbind.data.frame(Date,DF)
rng<-c(mroundd(min(DF[,-1],na.rm=T),0.1),mroundu(max(DF[,-1],na.rm=T),0.1))
if(rng[1] < 0){rng[1]<-0};  if(rng[2] > 1){rng[2]<-1}
ylab<-'Effective soil saturation [-]'
xlab<-'DOY'
lwd<-1.5

DF$yr<-year(DF$Date)
DF$doy<-yday(DF$Date)
rm(Date)

### PANEL 1: DECIDOUS ###
PFT<-'Decidous'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl1<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 2: EVERGREEN ###
PFT<-'Evergreen'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl2<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 3: GRASS ###
PFT<-'Grass'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl3<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)

### PANEL 4: MIXED EVERGREEN-DECIDOUS ###
PFT<-'Mixed-Evergreen-Decidous'
DAT<-DF[,which(names(DF) %in% c(PFT,'yr','doy'))]
DAT[,'yr']<-as.character(DAT[,'yr'])
colnames(DAT)[1]<-'value'
###
pl4<-ggplot(data=DAT)+
  geom_line(aes(x=doy,y=value,group=yr,colour=yr),linewidth=lwd)+
  scale_y_continuous(name=ylab,limits=rng)+
  labs(title=PFT,x=xlab)+
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
rm(PFT,DAT)


toptit<-textGrob(DOMAIN,
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3,pl4),nrow=4,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/VEG/daily_PFT-avg_Se_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=30,units='cm')
graphics.off()
rm(fn,DF,pl)




###
print('+++END SCRIPT <<PLOT PFT-AVG ANNUAL CYCLE>>+++   '&outnm)


