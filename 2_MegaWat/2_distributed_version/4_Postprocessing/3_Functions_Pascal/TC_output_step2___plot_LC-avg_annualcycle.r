################################################################################
# T&C OUTPUT :
#  PLOT SPATAL AVERAGE TIMESERIES PEr LAND COVEr CLASS
#
# 
# TODO: -
#       -
#
# NEW:  -
#       -
#
# 2023/09/14
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
# outnm<-'v20231006_Aaretal_2557days'
# outnm<-'v20231006_Aargauer_Rhein_2557days'
# outnm<-'v20231006_Berner_Oberland_2557days'
# outnm<-'v20231006_Bodensee_2557days'
# outnm<-'v20231006_Engadin_2557days'
# outnm<-'v20231006_Glatt_2557days'
# outnm<-'v20231006_Hinterrhein_2557days'
# outnm<-'v20231006_Jura_2557days'
# outnm<-'v20231006_Leman_2557days'
# outnm<-'v20231006_Limmat_2557days'
# outnm<-'v20231006_Neuchatel_2557days'
# outnm<-'v20231006_Oberwallis_2557days'
# outnm<-'v20231006_Reuss_2557days'
# outnm<-'v20231006_Rheintal_2557days'
# outnm<-'v20231006_Saane_2557days'
outnm<-'v20231006_Schaffhausen_2557days'
# outnm<-'v20231006_Ticino_2557days'
# outnm<-'v20231006_Thur_2557days'
# outnm<-'v20231006_Toess_2557days'
# outnm<-'v20231006_Unterwallis_369days'
# outnm<-'v20231006_Vorderrhein_2557days'
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
ts_init<-'2016-10-01 02:00:00'
ts_final<-'2022-09-30 23:00:00'
# ts_init<-'2015-10-01 02:00:00'
# ts_final<-'2016-08-01 23:00:00'

# ## define drought events (for plotting)
# dr1_init<-'2018-01-01 01:00:00'
# dr1_final<-'2018-12-31 23:00:00'
# ##
# dr2_init<-'2022-01-01 01:00:00'
# dr2_final<-'2022-09-30 23:00:00'

## define LC classes to compare among each other (names from "SPAVG-LC___XXX.txt" tables)
LCcomp_vars<-c('Evergreen','Mixed-Evergreen-Decidous','Decidous','Grass')


#===============================================================================
# PATHS
#===============================================================================
### Root
# root<-'N:/gebhyd/8_Him/Personal_folders/Pascal'
root<-'E:'

## Path to online storage
path_onst<-'C:/Users/Buri/switchdrive/WSL' 

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
library(extrafont)     #loadfonts()
library(lubridate)
library(assertthat)
library(grid)       #textGrob()
library(gridExtra)  #marrangeGrob()
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


#===============================================================================
# READ MODELLED LANDCOVER-AVERAGE DATA (SPAVG-LC___XXX)
#===============================================================================
setwd(path_data&'/TABLES')
fn<-list.files(pattern='SPAVG-LC___')
df_lc_ls<-list()
for(i in 1:length(fn)){
  FN<-fn[i]
  df_lc<-read.table(FN,header=T,sep='\t',dec='.')
  
  print('+++ reading <<'&FN&'>> +++')
  ##correct for midnight hour
  df_lc$Date<-as.character(df_lc$Date)
  idx<-which(nchar(df_lc$Date) < 19)
  df_lc$Date[idx]<-df_lc$Date[idx]&' 00:00:00'
  rm(idx)
  df_lc$Date<-as.POSIXct(df_lc$Date,format='%Y-%m-%d %H:%M:%S',tz=TZ)
  
  #derive monthly mean
  sumidx<-which(names(df_lc) %in% c('Pr','Pr_liq','Pr_sno','T'))
  df_avg<-aggregate(df_lc[,-c(1,sumidx)],list(Date=cut(df_lc$Date,'month')),mean) 
  df_avg$Date<-as.POSIXct(df_avg$Date,format='%Y-%m-%d')
  #derive monthly sum
  df_sum<-aggregate(df_lc[,sumidx],list(Date=cut(df_lc$Date,'month')),sum) 
  df_lc<-cbind.data.frame(df_avg,df_sum[,-1])
  rm(df_avg,df_sum)
  df_lc<-df_lc[which(df_lc[,1] >= ts_init & df_lc[,1] <= ts_final),] #crop to selected period
  df_lc_ls[[i]]<-df_lc
  
  ##name
  FN<-str_split(FN,'___')[[1]][2]
  FN<-str_split(FN,'\\.')[[1]][1]
  names(df_lc_ls)[i]<-FN
  rm(df_lc,FN)
}
rm(fn)


#===============================================================================
# READ MODELLED LANDCOVER-AVERAGE DATA (....LC_code_X)
#===============================================================================


#===============================================================================
# READ MODELLED PFT-AVERAGE DATA (....PFT_code_X)
#===============================================================================









newdir<-path_data&'/PLOTS'
dir.create(newdir)
dir.create(newdir&'/METEO')
#===============================================================================
# +++ PLOT LAND COVER AVERAGE +++ 
# PANEL 1: T
# PANEL 2: OF
# PANEL 3: V
#===============================================================================
T_ls<-list()
OF_ls<-list()
V_ls<-list()
for(i in 1:length(LCcomp_vars)){
  idx<-which(LCcomp_vars[i] == names(df_lc_ls))
  DF<-df_lc_ls[[idx]]
  if(i == 1){
    dates<-DF$Date
  }
  T_ls[[i]]<-DF[,which(names(DF) == 'T')]     ## T
  OF_ls[[i]]<-DF[,which(names(DF) == 'OF')]   ## OF
  V_ls[[i]]<-DF[,which(names(DF) == 'V')]     ## V
  rm(idx,DF)
}


### PANEL 1: T
DF<-cbind.data.frame(Date=dates,do.call('cbind',T_ls))
colnames(DF)[-1]<-LCcomp_vars
ylim<-mround(range(DF[,-1],na.rm=T),20)
###
ylab<-expression(paste("Transp.  [mm ",mo^-1,"]"))
brks<-seq(ylim[1],ylim[2],20)
leg<-''
tit<-''
sub<-''
capt<-''

DF<-melt(DF,'Date')
###
pl1<-ggplot(data=DF)+
  geom_line(aes(x=Date,y=value,colour=variable),linewidth=1.5)+
  scale_y_continuous(name=ylab,breaks=brks,
                     limits=ylim)+
  scale_color_jco()+
  labs(title=tit)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="top",
    legend.title = element_blank(),
    legend.text = element_text(colour='black',size=base_size),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(linewidth=0.1,colour='grey'),
    panel.grid.major.y = element_line(linewidth=0.75,colour='grey'),
    panel.grid.minor.x = element_line(linewidth=0.1,colour='grey'),
    panel.grid.major.x = element_line(linewidth=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl1


### PANEL 2: OF
DF<-cbind.data.frame(Date=dates,do.call('cbind',OF_ls))
colnames(DF)[-1]<-LCcomp_vars
ylim<-mround(range(DF[,-1],na.rm=T),0.05)
###
ylab<-'Soil moisture (top) [-]'
brks<-seq(ylim[1],ylim[2],0.05)
leg<-''
tit<-''
sub<-''
capt<-''

DF<-melt(DF,'Date')
###
pl2<-ggplot(data=DF)+
  geom_line(aes(x=Date,y=value,colour=variable),linewidth=1.5)+
  scale_y_continuous(name=ylab,breaks=brks,
                     limits=ylim)+
  scale_color_jco()+
  labs(title=tit)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(linewidth=0.1,colour='grey'),
    panel.grid.major.y = element_line(linewidth=0.75,colour='grey'),
    panel.grid.minor.x = element_line(linewidth=0.1,colour='grey'),
    panel.grid.major.x = element_line(linewidth=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl2


### PANEL 3: V
DF<-cbind.data.frame(Date=dates,do.call('cbind',V_ls))
colnames(DF)[-1]<-LCcomp_vars
ylim<-mround(range(DF[,-1],na.rm=T),0.5)
###
ylab<-'Soil water [mm]'
brks<-seq(ylim[1],ylim[2],0.5)
leg<-''
tit<-''
sub<-''
capt<-''

DF<-melt(DF,'Date')
###
pl3<-ggplot(data=DF)+
  geom_line(aes(x=Date,y=value,colour=variable),linewidth=1.5)+
  scale_y_continuous(name=ylab,breaks=brks,
                     limits=ylim)+
  scale_color_jco()+
  labs(title=tit)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="none",
    legend.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='white',colour='white'),
    panel.grid.minor.y = element_line(linewidth=0.1,colour='grey'),
    panel.grid.major.y = element_line(linewidth=0.75,colour='grey'),
    panel.grid.minor.x = element_line(linewidth=0.1,colour='grey'),
    panel.grid.major.x = element_line(linewidth=0.75,colour='grey'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl3


toptit<-textGrob(DOMAIN,
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3),nrow=3,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/METEO/monthly_LC_avg_T-OF-V_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=20,units='cm') 
graphics.off()
rm(fn,DF,pl)




###
print('+++END SCRIPT "PLOT LC-AVG ANNUAL CYCLE"+++   '&outnm)





