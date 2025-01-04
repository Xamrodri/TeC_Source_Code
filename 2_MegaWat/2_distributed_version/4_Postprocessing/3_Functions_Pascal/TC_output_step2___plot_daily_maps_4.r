################################################################################
# T&C OUTPUT:
#  READ DAILY MAPS FROM RASTER-FILES & PLOT 
#
# 
# TODO: -
# 
# NEW:  
#       
#
#
# 2023/08/29
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


#===============================================================================
# PATHS
#===============================================================================
### Root
# root<-'N:/gebhyd/8_Him/Personal_folders/Pascal'
root<-'E:'

## Path to online storage
path_onst<-'C:/Users/Buri/switchdrive/WSL' 

## Path to stored data
path_data<-root&'/T&C/OUTPUTSTORAGE/DISTRIBUTED/'&outnm

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
library(plotly)
library(doParallel)   #foreach()
library(doSNOW)

# define no. of digits printed in console
options('scipen'=100, 'digits'=4)

## define corresponding timezone
# TZ<-'Europe/Zurich'  #considers summer time (CET/CEST)
TZ<-'Africa/Algiers'  #does not have daylight saving

# set time zone & change local time from german to english for plotting
# (check in sessionInfo() ) 
Sys.setenv(TZ=TZ)
Sys.setlocale('LC_TIME', 'C')

# access functions file
source(path_func)


#===============================================================================
# READ SIMULATION FILE << AAA_TC_XXX.m >> OR << BBB_TC_XXX_reinit.m >>
#===============================================================================
setwd(path_data)
files<-list.files()
##read AAA-file
simfile<-files[startsWith(files,'AAA_')]
simfile<-readLines(simfile,warn=FALSE)
reinit=FALSE


#===============================================================================
# GET SIMULATION NAMES
#===============================================================================
VERSION<-grep(pattern="^VERSION",simfile)
VERSION<-simfile[VERSION]
VERSION<-str_split(VERSION,"'")[[1]][2]

DOMAIN<-grep(pattern="^DOMAIN",simfile)
DOMAIN<-simfile[DOMAIN]
DOMAIN<-str_split(DOMAIN,"'")[[1]][2]

simnm<-VERSION&'_'&DOMAIN
simnm1<-simnm

# ##check for BBB-file (re-initialization)
# simnm2<-files[startsWith(files,'BBB_')]
# if(length(simnm2) == 1){
#   reinit=TRUE
#   simfile2<-readLines(simnm2,warn=FALSE)
#   ##get secondary simulation name
#   ## VERSION
#   VERSION<-grep(pattern="VERSION",x=simfile2)
#   VERSION<-simfile2[VERSION[1]]
#   VERSION<-str_match(VERSION,"'\\s*(.*?)\\s*'")[,2]
#   ## DOMAIN
#   DOMAIN<-grep(pattern="DOMAIN",x=simfile2)
#   DOMAIN<-simfile2[DOMAIN[1]]
#   DOMAIN<-str_match(DOMAIN,"'\\s*(.*?)\\s*'")[,2]
#   simnm2<-VERSION&'_'&DOMAIN
#   simnm<-simnm2
# }
# rm(files)

if(reinit == F){simf<-simfile}
if(reinit == T){simf<-simfile2}


#===============================================================================
# GET INPUT FOLDER NAME
#===============================================================================
idx<-grep(pattern="^GEODATAINPUT = \\[",simf)
fo<-simf[idx]
fo<-str_split(fo,"___v")[[1]][2]
inputfo<-path_input&'/GeodataInput_'&DOMAIN&'___v'&as.character(gsub("\\D","",fo))
rm(idx,fo)


#===============================================================================
# LOAD GEODATA & GEOREF INPUT-FILES << Geodata_XXX.mat >> & << GEOREF_XXX.mat >>
#===============================================================================
setwd(inputfo)
files<-list.files()
geodatafile<-readMat(files[which(str_starts(files,'Geodata_'))])

### RASTER DIMENSIONS & coordinates
m_cell<-geodatafile[['m']][1,1]
n_cell<-geodatafile[['n']][1,1]
resol<-geodatafile[['cellsize']][1,1]
x<-geodatafile[['x']][,1]
y<-geodatafile[['y']][,1]

### PROJECTION
georef<-readLines(files[which(str_starts(files,'GEOREF_'))])
epsg<-str_split(georef[[1]][1],'-')[[1]]
epsg<-as.numeric(epsg[length(epsg)])
projec<-make_EPSG()
projec<-projec[which(projec$code == epsg),'prj4']

### EXTENT
ext<-c(min(x)-resol/2,max(x)+resol/2,min(y)-resol/2,max(y)+resol/2)

### BASIN MASK
msk<-geodatafile[['MASK']]
msk[msk < 1]<-NA
msk<-apply(msk,2,rev)   ## Matlab: flipud(DTM)
#check
# plot(raster(msk))
msk<-raster(nrow=m_cell,ncol=n_cell,vals=msk)
projection(msk)<-as.character(projec)
extent(msk)<-ext

### SOIL THICKNESS
soilth<-geodatafile[['SOIL.TH']]
soilth[soilth < 0]<-NA
soilth<-apply(soilth,2,rev)   ## Matlab: flipud(DTM)
#check
# plot(raster(soilth))
soilth<-raster(nrow=m_cell,ncol=n_cell,vals=soilth)
projection(soilth)<-as.character(projec)
extent(soilth)<-ext


#===============================================================================
# READ HILLSHADE
#===============================================================================
hs<-raster(fn_hs)
hs<-crop(hs,ext)


#===============================================================================
# CREATE NEW FOLDERS
#===============================================================================
new_dir<-path_store&'/'&outnm
dir.create(new_dir&'/PLOTS')
dir.create(new_dir&'/PLOTS/DAILY')



setwd(new_dir&'/RASTERS/DAILY')

#===============================================================================
# PLOT ONE VARIABLE & MAP PER DAY
#===============================================================================

###### V ########
VOI<-'OF'
# VOI<-'O450'

files<-list.files()
br<-brick(files[which(startsWith(files,VOI))[1]])


###testing
# levelplot(br[[1:5]])


# ##APR-OCT 2016
# br<-br[[184:367]]

# ##JAN 2018-DEC 2020
# br<-br[[824:1919]]

# ## 1 APR- 31 OCT 2018
br<-br[[914:1127]]

# ## 1OCT 2020-30SEP 2021
# br<-br[[1828:2192]]

## 1OCT 2021-30SEP 2022
# br<-br[[2193:2557]]


dates<-str_split(names(br),'X')
dates<-do.call('cbind',dates)[2,]
dates<-str_replace(dates,'[.]',' ')
dates<-str_replace(dates,'[.]',' ')
DATES<-as.POSIXct(dates,format='%d %b %Y')

dir.create(new_dir&'/PLOTS/DAILY/'&VOI)


mx<-max(values(br))
mx<-mroundu(mx,0.1)

pal<- rev(plasma(10))

for(j in 1:nlayers(br)){
  r<-br[[j]]
  r[msk == 0]<-NA
  dat<-dates[j]
  tit<-VOI&' [-] \n '&dat
  brks<-c(0.01,seq(0.1,mx,0.1))
  
  print('+++ plotting daily <<'&VOI&'>> for '&DOMAIN&' ('&dat&') +++')
  ### Plot maps with fixed color scale (for GIFs...)
  fn<-new_dir&'/PLOTS/DAILY/'&VOI&'/fxd_'&VOI&'_ts'&j&'.png'
  png(file=fn,units='in',width=5,height=5,res=300)
  plot(hs,col=grey(1:100/100),legend=F,main=tit,axes=FALSE)
  plot(r,axes=FALSE,box=FALSE,col=pal,breaks=brks,add=T,alpha=0.6,colNA=NA)
  # plot(gla_p,border='white',add=T)
  # plot(deb_p,border='red',add=T)
  dev.off()
  rm(fn)
}




#===============================================================================
# PLOT MULTIPLE VARIABLES & MAPS PER DAY
#===============================================================================
VOIS<-c('OF','O125','O450')

files<-list.files()
br_ls<-list()
mx_v<-list()
for(i in 1:length(VOIS)){
  
  print('+++ subsetting daily <<'&VOIS[i]&'>> for '&DOMAIN&' +++')
  br<-brick(files[which(startsWith(files,VOIS[i]))[1]])
  
  ###testing
  # levelplot(br[[1:5]])
  
  # br<-br[[184:367]]       ##APR-OCT 2016
  # br<-br[[824:1919]]      ##JAN 2018-DEC 2020
  # br<-br[[914:1127]]        ##1 APR- 31 OCT 2018
  br<-br[[1279:1492]]        ##1 APR- 31 OCT 2019
  # br<-br[[1828:2192]]     ##1 OCT 2020-30 SEP 2021
  # br<-br[[2193:2557]]     ##1 OCT 2021-30 SEP 2022
  
  if(i == 1){
  dates<-str_split(names(br),'X')
  dates<-do.call('cbind',dates)[2,]
  dates<-str_replace(dates,'[.]',' ')
  dates<-str_replace(dates,'[.]',' ')
  DATES<-as.POSIXct(dates,format='%d %b %Y')
  }
  
  mx_v[i]<-max(values(br))

  br_ls[[i]]<-br
  rm(br)
}
names(br_ls)<-VOIS

voi_nms<-paste(VOIS,sep="",collapse="-")

dir.create(new_dir&'/PLOTS/DAILY/MULTIPLE')

mx<-mroundu(max(unlist(mx_v)),0.1)

pal<- rev(plasma(10))

for(j in 1:nlayers(br_ls[[1]])){
  r<-brick(br_ls[[1]][[j]],br_ls[[2]][[j]],br_ls[[3]][[j]])
  names(r)<-VOIS
  r[msk == 0]<-NA
  dat<-dates[j]
  tit<-'Soil moisture [-] \n '&dat
  brks<-c(0.01,seq(0.1,mx,0.1))
  
  print('+++ plotting daily <<'&voi_nms&'>> for '&DOMAIN&' ('&dat&') +++')
  ### Plot maps with fixed color scale (for GIFs...)
  fn<-new_dir&'/PLOTS/DAILY/MULTIPLE/fxd_'&voi_nms&'_ts'&j&'.png'
  png(file=fn,units='in',width=10,height=5,res=300)
  
  
  print(levelplot(r,main=tit))

  dev.off()
  rm(fn)
}






#===============================================================================
# PLOT ONE VARIABLE & MAP PER DAY
#===============================================================================
# 
# 
# ###### CSNO ########
# VOI<-'Csno'
# 
# files<-list.files()
# br<-brick(files[which(startsWith(files,VOI))[1]])
# 
# ##APR-OCT 2016
# br<-br[[184:367]]
# 
# dates<-str_split(names(br),'X')
# dates<-do.call('cbind',dates)[2,]
# dates<-str_replace(dates,'[.]',' ')
# dates<-str_replace(dates,'[.]',' ')
# DATES<-as.POSIXct(dates,format='%d %b %Y')
# 
# dir.create(new_dir&'/PLOTS/DAILY/'&VOI)
# 
# mx<-1
# 
# pal<- rev(viridis(10))
# 
# for(j in 1:nlayers(br)){
#   r<-br[[j]]
#   r[msk == 0]<-NA
#   
#   if(VOI == 'Csno'){
#     r<-r/24
#     r[r < 1]<-NA
#     }
#   
#   
#   dat<-dates[j]
#   tit<-'Snow cover \n '&dat
#   brks<-seq(0,1,0.1)
#   
#   ### Plot maps with fixed color scale (for GIFs...)
#   fn<-new_dir&'/PLOTS/DAILY/'&VOI&'/fxd_'&VOI&'_ts'&j&'.png'
#   png(file=fn,units='in',width=5,height=5,res=300)
#   plot(hs,col=grey(1:100/100),legend=F,main=tit,axes=FALSE)
#   plot(r,axes=FALSE,box=FALSE,col=pal,breaks=brks,add=T,alpha=0.6,colNA=NA,legend=F)
#   # plot(gla_p,border='white',add=T)
#   # plot(deb_p,border='red',add=T)
#   dev.off()
#   rm(fn)
# }



###
print('+++END SCRIPT "PLOT DAILY MAPS"+++   '&outnm)


