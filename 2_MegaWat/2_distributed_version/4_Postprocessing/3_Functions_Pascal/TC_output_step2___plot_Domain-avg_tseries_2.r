################################################################################
# T&C OUTPUT :
#  PLOT SPATAL AVERAGE TIMESERIES OVER DOMAIN
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
# outnm<-'v20231006_Aaretal_404days'
# outnm<-'v20231006_Aargauer_Rhein_2557days'
# outnm<-'v20231006_Berner_Oberland_779days'
# outnm<-'v20231006_Bodensee_2557days'
outnm<-'v20231006_Engadin_1906days'
# outnm<-'v20231006_Glatt_2557days'
outnm<-'v20231006_Hinterrhein_1853days'
outnm<-'v20231006_Jura_1562days'
# outnm<-'v20231006_Leman_932days'
# outnm<-'v20231006_Limmat_841days'
# outnm<-'v20231006_Neuchatel_813days'
outnm<-'v20231006_Oberwallis_1498days'
outnm<-'v20231006_Reuss_1257days'
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


setwd(path_data&'/RASTERS/INIT-COND')
#===============================================================================
# READ LC & DEM
#===============================================================================
lc<-raster('VEG_CODE.tif')
dem<-raster('DTM.tif')
ele_rng<-range(values(dem),na.rm=T)


#===============================================================================
# READ MODELLED DOMAIN-AVERAGE DATA
#===============================================================================
fn<-path_data&'/TABLES/SPAVG.txt'
df<-read.table(fn,header=T,sep='\t',dec='.')
rm(fn)

##correct for midnight hour
df$Date<-as.character(df$Date)
idx<-which(nchar(df$Date) < 19)
df$Date[idx]<-df$Date[idx]&' 00:00:00'
rm(idx)
df$Date<-as.POSIXct(df$Date,format='%Y-%m-%d %H:%M:%S',tz=TZ)

#derive monthly mean
sumidx<-which(names(df) %in% c('Pr','Pr_liq','Pr_sno','T'))
df_avg<-aggregate(df[,-c(1,sumidx)],list(Date=cut(df$Date,'month')),mean) 
df_avg$Date<-as.POSIXct(df_avg$Date,format='%Y-%m-%d')
#derive monthly sum
df_sum<-aggregate(df[,sumidx],list(Date=cut(df$Date,'month')),sum) 
df<-cbind.data.frame(df_avg,df_sum[,-1])
rm(df_avg,df_sum)
df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period


newdir<-path_data&'/PLOTS'
dir.create(newdir)
dir.create(newdir&'/METEO')
#===============================================================================
# +++ PLOT DOMAIN AVERAGE +++ 
# PANEL 1: Ta & Pr
# PANEL 2: Ds & Rsw
# PANEL 3: Csno & SWE
#===============================================================================
DF<-df
DF<-DF[,which(names(DF) %in% c('Date','Pr','Ta','Csno','SWE','Ds','Rsw'))]

### PANEL 1: Ta & Pr ###
ylim<-c(mroundd(min(DF$Ta,na.rm=T),5),mroundu(max(DF$Ta,na.rm=T),5))
ylim2<-c(0,mroundu(max(DF$Pr,na.rm=T),200))
fact1<-max(ylim)/max(ylim2) #factor for secondary axis
###
ylab=expression(paste(T['air'],"  [°C]"))
ylab2<-expression(paste("Precipitation  [mm ",mo^-1,"]"))
brks<-seq(ylim[1],ylim[2],5)
brks2<-seq(ylim2[1],ylim2[2],200)
col1<-'firebrick3'
col2<-'navyblue'
leg<-''
tit<-''
sub<-''
capt<-'Avg.: '&round(mean(DF$Ta,na.rm=T),1)&' °C  |  '&round(mean(DF$Pr,na.rm=T),1)&' mm mo-1'
###
pl1<-ggplot(data=DF)+
  geom_bar(aes(x=Date,y=Pr*fact1),stat='identity',colour="white",fill=col2)+ #colour="#636363")+
  geom_hline(yintercept=0,color='black',linewidth=0.5)+
  geom_line(aes(x=Date,y=Ta),color=col1,linewidth=2)+
  scale_y_continuous(name=ylab,breaks=brks,
                     limits=ylim,
                     sec.axis=sec_axis(trans=~./fact1,name=ylab2,breaks=brks2))+
  labs(title=tit,caption=capt)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="none",
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour=col1,size=base_size),
    axis.title.y = element_text(colour=col1,size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    axis.text.y.right = element_text(colour=col2,size=base_size),
    axis.title.y.right = element_text(colour=col2,size=base_size),
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


### PANEL 2: Ds & Rsws 
ylim<-c(mroundd(min(DF$Ds,na.rm=T),200),mroundu(max(DF$Ds,na.rm=T),200))
ylim2<-c(mroundd(min(DF$Rsw,na.rm=T),50),mroundu(max(DF$Rsw,na.rm=T),50))
fact2<-max(ylim)/max(ylim2) #factor for secondary axis
###
ylab="VPD  [Pa]"
ylab2<-expression(paste("SW  [W ",m^-2,"]"))
brks<-seq(ylim[1],ylim[2],200)
brks2<-seq(ylim2[1],ylim2[2],50)
col1<-'royalblue3'
col2<-'darkorange1'
leg<-''
tit<-''
sub<-''
capt<-'Avg.: '&round(mean(DF$Ds,na.rm=T),1)&' Pa  |  '&round(mean(DF$Rsw,na.rm=T),1)&' W m-2'
###
pl2<-ggplot(data=DF)+
  geom_line(aes(x=Date,y=Ds),color=col1,linewidth=2)+
  geom_line(aes(x=Date,y=Rsw*fact2),color=col2,linewidth=2)+
  scale_y_continuous(name=ylab,breaks=brks,
                     limits=ylim,
                     sec.axis=sec_axis(trans=~./fact2,name=ylab2,breaks=brks2))+
  labs(title=tit,caption=capt)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="none",
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour=col1,size=base_size),
    axis.title.y = element_text(colour=col1,size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    axis.text.y.right = element_text(colour=col2,size=base_size),
    axis.title.y.right = element_text(colour=col2,size=base_size),
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


### PANEL 3: Csno & SWE
ylim<-c(0,1)
ylim2<-c(mroundd(min(DF$SWE,na.rm=T),200),mroundu(max(DF$SWE,na.rm=T),200))
fact3<-max(ylim2)/max(ylim) #factor for secondary axis
###
ylab="Snow covered fraction [-]"
ylab2<-expression(paste("SWE  [mm ",mo^-1,"]"))
brks<-seq(ylim[1],ylim[2],0.2)
brks2<-seq(ylim2[1],ylim2[2],200)
col1<-'saddlebrown'
col2<-'chartreuse4'
leg<-''
tit<-''
sub<-''
capt<-'Avg.: '&round(mean(DF$Csno,na.rm=T),2)&' [-]  |  '&round(mean(DF$SWE,na.rm=T),1)&' mm'
###
pl3<-ggplot(data=DF)+
  geom_line(aes(x=Date,y=Csno),color=col1,linewidth=2)+
  geom_line(aes(x=Date,y=SWE/fact3),color=col2,linewidth=2)+
  scale_y_continuous(name=ylab,breaks=brks,
                     limits=ylim,
                     sec.axis=sec_axis(trans=~.*fact3,name=ylab2,breaks=brks2))+
  labs(title=tit,caption=capt)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="none",
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour=col1,size=base_size),
    axis.title.y = element_text(colour=col1,size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    axis.text.y.right = element_text(colour=col2,size=base_size),
    axis.title.y.right = element_text(colour=col2,size=base_size),
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
# pl_ls[[3]]<-pl3
# rm(pl3)

toptit<-textGrob(DOMAIN&' ('&ele_rng[1]&'-'&ele_rng[2]&' m a.s.l.)',
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2,pl3),nrow=3,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/METEO/monthly_avg_Ta-P-Ds-Rsw-Csno-SWE_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=25,height=25,units='cm') 
graphics.off()
rm(fn,DF,pl)



dir.create(newdir&'/BGW')
#===============================================================================
# +++ PLOT DOMAIN AVERAGE +++ 
# 
# PANEL 1: ABSOLUTE BGW
# PANEL 2: RELATIVE BGW
#===============================================================================
DF<-cbind.data.frame(Date=df$Date,G=df$EG+df$EIn+df$EWAT+df$T,W=df$EICE+df$ESN)
DF[which(DF$W < 0),'W']<-0
DF$B<-df$Pr-DF$G-DF$W
DF[which(DF$B < 0),'B']<-0
###
DF$TOT<-rowSums(DF[,-1])
DF$B_rel<-DF$B/DF$TOT
DF$G_rel<-DF$G/DF$TOT
DF$W_rel<-DF$W/DF$TOT


### PANEL 1: ABSOLUTE BGW
ylim<-c(0,mroundu(max(DF[,2:4],na.rm=T),50))
###
ylab<-expression(paste("[mm ",mo^-1,"]"))
brks<-seq(ylim[1],ylim[2],50)
leg<-''
tit<-''
sub<-''
capt<-'Avg. [mm mo-1]: '&round(mean(DF$B,na.rm=T),1)&' (Blue) |  '&round(mean(DF$G,na.rm=T),1)&' (Green)  |  '&round(mean(DF$W,na.rm=T),1)&' (White)'
###
pl1<-ggplot(data=DF)+
  geom_line(aes(x=Date,y=B),color='dodgerblue4',linewidth=2)+
  geom_line(aes(x=Date,y=G),color='chartreuse4',linewidth=2)+
  geom_line(aes(x=Date,y=W),color='white',linewidth=2)+
  scale_y_continuous(name=ylab,breaks=brks,
                     limits=ylim)+
  labs(title=tit,caption=capt)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    # legend.position="none",
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_text(colour='black',size=base_size),
    axis.title.x = element_blank(),
    panel.background = element_rect(fill='grey',colour='grey'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey30'),
    panel.grid.major.y = element_line(size=0.75,colour='grey30'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey30'),
    panel.grid.major.x = element_line(size=0.75,colour='grey30'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl1
dev.off()

### PANEL 2: RELATIVE BGW
DF2<-melt(DF[,c(1,6:8)],'Date')
ylim<-c(0,1.01)
###
ylab<-'[-]'
brks<-seq(0,1,0.2)
leg<-''
tit<-''
sub<-''
###
pl2<-ggplot(DF2,aes(x=Date,y=value,fill=variable)) + 
  geom_area(alpha=0.7)+
  scale_y_continuous(name=ylab,breaks=brks,
                     limits=ylim)+
  scale_fill_manual(name='',labels=c('Blue','Green','White'),values=c('dodgerblue4','chartreuse4','white'))+
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
    panel.background = element_rect(fill='grey',colour='grey'),
    panel.grid.minor.y = element_line(size=0.1,colour='grey30'),
    panel.grid.major.y = element_line(size=0.75,colour='grey30'),
    panel.grid.minor.x = element_line(size=0.1,colour='grey30'),
    panel.grid.major.x = element_line(size=0.75,colour='grey30'),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    plot.margin=unit(c(0,1,0,0.2),'lines') #top/right/bottom/left
  ) 
pl2


toptit<-textGrob(DOMAIN,
                 gp=gpar(fontsize=base_size+5,font=2))
pl<-marrangeGrob(list(pl1,pl2),nrow=2,ncol=1,top=toptit) ##n=12
pl

fn<-newdir&'/BGW/monthly_avg'&PERIOD_SH&'.png'
ggsave(fn,pl,width=20,height=28,units='cm') 
graphics.off()
rm(fn,DF,pl)



###
print('+++END SCRIPT <<PLOT DOMAIN-AVG TSERIES>>+++   '&outnm)


