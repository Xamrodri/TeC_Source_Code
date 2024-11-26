################################################################################
# T&C OUTPUT STEP 2:
#  DERIVE WATER BALANCE COMPONENTS CONTRIBUTION TO STREAMNETWORK-POINTS
#  (after running "TC_output_step2___derive_streams_uparea.R")
# 
# TODO: -
#
#
# 2024/01/17
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
## define number of pourpoints to consider (if 'NA' take the max. available)
nppt<-NA
# nppt<-10

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
library(terra)
library(sf)
library(raster)
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
# READ DEM
#===============================================================================
## DEM
dem<-rast(path_data&'/RASTERS/INIT-COND/DTM.tif')
resol<-res(dem)[1]
ext_dem<-extent(raster::raster(dem))
slo<-terra::terrain(dem,v='slope',unit='degrees')
asp<-terra::terrain(dem,v='aspect',unit='degrees')


#===============================================================================
# RETRIEVE DATES FROM SIMULATION
#===============================================================================
#read catchment tables
setwd(path_data&'/TABLES')
spavg_df<-read.table('SPAVG.txt',header=T,sep='\t')
spavg_df$Date<-as.POSIXct(spavg_df$Date,format='%Y-%m-%d %H:%M:%S')


#===============================================================================
# TIME PERIOD OF INTEREST
#===============================================================================
ts_init<-as.POSIXct(ts_init,format='%Y-%m-%d %H:%M:%S')
ts_final<-as.POSIXct(ts_final,format='%Y-%m-%d %H:%M:%S')

t1<-spavg_df$Date[1]
tn<-spavg_df$Date[nrow(spavg_df)]
tn<-floor_date(tn,'month')         ####hack!!
ini<-round(as.numeric(ts_init - t1)/30.42)+1
fin<-round(as.numeric(tn - t1)/30.42)
NMONTHS<-ini:fin
n<-length(NMONTHS)
rm(t1,ini,fin)

#index for selected period
idx_ts<-which(spavg_df$Date >= ts_init & spavg_df$Date <= tn)
dates_all<-spavg_df$Date
dates_period<-spavg_df$Date[idx_ts]
dates_period<-dates_period[1:(length(dates_period)-3)]

PERIOD_SH<-format(ts_init,"%b%Y")&'-'&format(tn-86400,"%b%Y")
PERIOD_LO<-format(ts_init,"%d %b %Y")&' - '&format(tn-86400,"%d %b %Y")

year_month_v<-unique(paste0(year(dates_period),'-',month(dates_period)))
Dates<-seq(ts_init,tn+86400,'1 month')
nDays<-difftime(Dates[-1],Dates[-length(Dates)],unit="day")
nHours<-difftime(Dates[-1],Dates[-length(Dates)],unit="hour")


#===============================================================================
# READ POURPOINTS & UPSTREAM AREAS ALONG STREAMNETWORK
#===============================================================================
setwd(path_data&'/RASTERS/UPAREA')
fn<-list.files(pattern='uparea_id___n')
if(is.na(nppt)){fn<-fn[length(fn)]}   #take the file with the highest 'n'
if(!is.na(nppt)){
  idx<-sapply(str_extract_all(fn,"\\d+"),tail,1)
  fn<-fn[which(idx == nppt)]
  rm(idx)
}
df<-read.table(fn,header=T)

fn<-str_replace(fn,'_id___','_br___')
fn<-str_replace(fn,'.check','.grd')
br<-raster::brick(fn)
br<-terra::rast(br)
rm(fn)

sta<-ncol(df)+1
Npt<-nrow(df)


#===============================================================================
# DERIVE MONTHLY CONTRIBUTIONS TO POURPOINTS FROM EACH VARIABLE
#===============================================================================
setwd(path_data&'/RASTERS/MONTHLY')
fls<-list.files()

# CKt = (V_tgtm1 - V_tg) + (Vice_tgtm1 - Vice_tg) + Pr_tg - EG_tg*dth - T_tg*dth - EIn_tg*dth - Lk_tg...
# - ESN_tg*dth - EIn_urb_tg*dth - EWAT_tg*dth - EIn_rock_tg*dth - EICE_tg*dth ...
# + (SWE_tgtm1 -SWE_tg) + (In_tgtm1 -In_tg) ...
# +  (ICE_tgtm1 -ICE_tg) +  (WAT_tgtm1 -WAT_tg)  +  (FROCK_tgtm1 -FROCK_tg) ...
# + Qlat_in_tgtm1 + q_runon_tgtm1 + Q_channel_tgtm1  ...
# - Qlat_in_tg - Q_exit - Qsub_exit - q_runon_tg - Q_channel_tg -Swe_exit;

##define variable names of interest (to write out)
nm_v<-c('V','Vice','Pr','EG','T','EIn','Lk_tg','ESN','EIn_urb','EWAT','EIn_rock',
        'EICE','SWE_tg','In','ICE','WAT','FROCK','Qlat_in','q_runon','Q_channel',
        'Q_exit','Qsub_exit','Swe_exit',
        'Pr_liq','Pr_sno','Imelt','Smelt')
N<-length(nm_v)
#################################
for(i in 1:N){
  nm<-nm_v[i]
  print(DOMAIN&': calculating avg. << '&nm&' >> contribution per ppt (N='&Npt&'), variable '&i&'/'&N)
  # fn<-fls[str_detect(fls,paste0('^',nm,'___.*','.grd','$'))]
  fn<-fls[str_detect(fls,paste0('^',nm,'___.*','.grd','$'))]
  # r<-rast(fn)[[NMONTHS]]  
  r<-brick(fn)[[NMONTHS]] 
  # r<-rast(fn,lyrs=1:10)
  r[is.na(r)]<-0          ###NAs from RhiresD?
  DF<-df
  for(j in 1:n){
    print('+++++++++++++++++++++++++++++++++++++++++++ month '&j&'/'&n)
    # DF<-cbind.data.frame(DF,global(r[[i]]*br,'mean',na.rm=T)[,1])
    DF<-cbind.data.frame(DF,global(r[[i]]*br,'mean',na.rm=T)[,1])
  }
  DF<-round(DF,2)
  colnames(DF)[sta:ncol(DF)]<-year_month_v
  write.table(DF,path_data&'/TABLES/'&'uparea-contrib_'&nm&'___'&PERIOD_SH&'___n'&Npt&'.txt',
              dec='.',sep='\t',quote=FALSE,col.names=TRUE,row.names=FALSE)
  rm(nm,r,DF)
}



###
print('+++END SCRIPT <<DERIVE STREAM CONTRIBUTIONS>>+++   '&outnm)



