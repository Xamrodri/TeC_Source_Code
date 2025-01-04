################################################################################
# T&C OUTPUT STEP 2:
#  PLOT WATER BALANCE COMPONENTS CONTRIBUTION TO STREAMNETWORK-POINTS
#  (after running "TC_output_step2___derive_streams_uparea.R"
#   & "TC_output_step2___derive_streams_contributions.R")
# 
# TODO: -
#
#
# 2023/12/20
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
ts_final<-'2022-09-30 23:00:00'    


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
path_gis<-root&'/QGIS/TC_CH'


##define filename of hillshade
fn_hs<-path_gis&'/DEM/HS_AW3D250_EPSG-2056.tif'


#===============================================================================
# LOAD PACKAGES
#===============================================================================
# define packages to load
library(terra)
library(sf)
library(sp)
library(raster)
library(reshape2)
library(ggplot2)
library(plyr)
library(lattice)
library(ggnewscale)
library(ggspatial)    #annotation_scale()
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

hs<-terra::rast(fn_hs)
crs(hs)<-crs(dem)
hs<-crop(hs,dem)
hs<-mask(hs,dem)
# plot(hs,col=grey(1:100/100),legend=F,main="",axes=FALSE)


hs_spdf <- as(raster(hs), "SpatialPixelsDataFrame")
hs_spdf <- as.data.frame(hs_spdf)
colnames(hs_spdf) <- c("value", "x", "y")


dem_spdf <- as(raster(dem), "SpatialPixelsDataFrame")
dem_spdf <- as.data.frame(dem_spdf)
colnames(dem_spdf) <- c("value", "x", "y")


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
years_v<-unique(year(dates_period))
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
N<-nrow(df)
rm(fn)


# tt<-calc(stack(br),sum,na.rm=T)
# plot(tt)


#===============================================================================
# READ MONTHLY CONTRIBUTIONS TO POURPOINTS FROM EACH VARIABLE
#===============================================================================
setwd(path_data&'/TABLES')
fls<-list.files()

##################################
nm<-'Pr_liq'
RA<-read.table(path_data&'/TABLES/'&'uparea-contrib_'&nm&'___'&PERIOD_SH&'___n'&N&'.txt',
            header=TRUE)
rm(nm)

##################################
nm<-'Smelt'
SM<-read.table(path_data&'/TABLES/'&'uparea-contrib_'&nm&'___'&PERIOD_SH&'___n'&N&'.txt',
                   header=TRUE)
rm(nm)

##################################
nm<-'Imelt'
IM<-read.table(path_data&'/TABLES/'&'uparea-contrib_'&nm&'___'&PERIOD_SH&'___n'&N&'.txt',
                   header=TRUE)
rm(nm)


#===============================================================================
# BOXPLOTS RELATIVE ICE MELT (PER YEAR & ELEVATION BAND)
#===============================================================================
new_dir<-path_data&'/PLOTS/STREAMCONTRIB'
dir.create(new_dir)

yr<-years_v
mo<-5:10  ##summer period

DF2_ls<-list()
for(i in 1:length(yr)){
  YR<-yr[i]
  DF1_ls<-list()
  for(j in 1:length(mo)){
    NM<-'X'&YR&'.'&mo[j]
    if(NM %in% names(IM)){
      DF1<-cbind.data.frame(df,ra=RA[,NM],sm=SM[,NM],im=IM[,NM])
      DF1$ru<-DF1$ra+DF1$sm+DF1$im ## total runoff
      DF1$rel_im<-DF1$im/DF1$ru  ## rel. ice melt contribution
      DF1$rel_sm<-DF1$sm/DF1$ru  ## rel. snow melt contribution
      DF1$rel_ra<-DF1$ra/DF1$ru  ## rel. rain contribution
      
      # ##check
      # ppts<-vect(DF1,geom=c("x","y"))
      # spplot(ppts, c("rel_ra", "rel_sm", "rel_im"))
      
      DF1$yr<-YR
      DF1$mo<-mo[j]
      DF1_ls[[j]]<-DF1
      rm(DF1)
    }
    if(NM %notin% names(IM)){DF1_ls[[j]]<-NULL}
    rm(NM)
  }
  DF2_ls[[i]]<-as.data.frame(do.call('rbind',DF1_ls))
  rm(YR,DF1_ls)
}

DF<-as.data.frame(do.call('rbind',DF2_ls))
DF$mo_nm<-month.abb[DF$mo]
DF$EB<-round_any(DF$DTM,500,floor)
DF$EB_lab<-as.character(DF$EB)&'-'&as.character(DF$EB+499)&' m a.s.l.'
DF$EB<-factor(DF$EB)
DF$yr<-factor(DF$yr)
##sorting elevation band factors
DF$EB_lab<-factor(DF$EB_lab)
idx<-unlist(lapply(str_split(levels(DF$EB_lab),'-'),function(x) x[-2]))
idx<-sort(as.numeric(idx),index.return=T)$ix
DF$EB_lab<-factor(DF$EB_lab,levels=levels(DF$EB_lab)[idx])

DF$mo_nm<-factor(DF$mo_nm,month.abb[mo])
pl<-ggplot(DF,aes(x=EB,y=rel_im,fill=mo_nm))+
  geom_boxplot()+
  # scale_fill_discrete(breaks=month.abb[mo])+
  scale_fill_viridis(discrete=T, name="",option='inferno',breaks=month.abb[mo])+
  facet_grid(yr ~ EB_lab, scales="free")+
  labs(x="",y='Ice melt contribution [-]',title=DOMAIN)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.text = element_text(colour='black',size=base_size),
    axis.ticks = element_blank(),
    axis.text.y = element_text(colour='black',size=base_size),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(colour='black',size=base_size),
    strip.text.y = element_text(colour='black',size=base_size),
    plot.margin=unit(c(0.2,1,0.2,0.2),'lines') #top/right/bottom/left
  ) 
pl

fn<-new_dir&'/Boxplots_IMcontrib_mo'&range(mo)[1]&'-'&range(mo)[2]&'___sim'&PERIOD_SH&'.png'
ggsave(fn,pl,width=30,height=25,units='cm') 
graphics.off()
rm(fn,pl)


#===============================================================================
# MAPPLOTS RELATIVE ICE MELT (PER YEAR & MONTH)
#===============================================================================
df<-DF
df$my<-paste0(df$mo_nm,'_',df$yr)

pl<-ggplot() +
  # geom_tile(data=hs_spdf,aes(x=x,y=y,fill=value),show.legend=F)+
  # scale_fill_gradient(low="black",high="white")+
  # new_scale_fill()+
  geom_point(data=df,aes(x=x,y=y,color=rel_im))+
  scale_color_viridis(name="Ice melt \n contrib. \n [-]",option='inferno')+
  facet_grid(yr ~ mo_nm, scales="free")+
  labs(x='',y='',title=DOMAIN)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    aspect.ratio = 1,
    legend.position="right",
    legend.text = element_text(colour='black',size=base_size),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_text(colour='black',size=base_size),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(colour='black',size=base_size),
    strip.text.y = element_text(colour='black',size=base_size),
    plot.margin=unit(c(0.2,1,0.2,0.2),'lines') #top/right/bottom/left
  ) 
pl

fn<-new_dir&'/Mapplots_IMcontrib_mo'&range(mo)[1]&'-'&range(mo)[2]&'___sim'&PERIOD_SH&'.png'
ggsave(fn,pl,width=30,height=25,units='cm') 
graphics.off()
rm(fn,pl)


#===============================================================================
# MAPPLOT DEM & POINTS
#===============================================================================
attr<-df[1:N,1:6]

pl<-ggplot() +
  geom_tile(data=hs_spdf,aes(x=x,y=y,fill=value),show.legend=F)+
  scale_fill_gradient(low="black",high="white")+
  new_scale_fill()+
  geom_tile(data=dem_spdf,aes(x=x,y=y,fill=value),alpha=0.4)+
  scale_fill_viridis(name="Elevation \n [m a.s.l.]",option='viridis')+
  geom_point(data=attr,aes(x=x,y=y),colour='black')+
  coord_equal()+
  annotation_scale(text_cex=1.5,text_family=base_family)+
  labs(x='Eastings [m]',y='Northings [m]',title=DOMAIN)+
  theme_gray(base_size=base_size,base_family=base_family)+
  theme(
    legend.position="right",
    legend.text = element_text(colour='black',size=base_size),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_text(colour='black',size=base_size),
    strip.text.y = element_text(colour='black',size=base_size),
    plot.margin=unit(c(0.2,1,0.2,0.2),'lines') #top/right/bottom/left
  ) 
pl

fn<-new_dir&'/DEM_and_streampoints.png'
ggsave(fn,pl,width=15,height=15,units='cm') 
graphics.off()
rm(fn,pl)





###
print('+++END SCRIPT <<PLOT STREAM CONTRIBUTIONS>>+++   '&outnm)



