################################################################################
# T&C-MULTIPOINT OUTPUT:
#  ANALYZE SOILPROFILE SITES CLIMATOLOGY
#
# 
# TODO: -
# 
# NEW:  
#       
#
#
# 15/08/2023
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
outnm<-'v20230727_SP'

# choose site
SITE<-'BLH'  ### BÜELENHORN
SITE<-'SWZ'  ### SCHWARZHORN


#===============================================================================
# PATHS
#===============================================================================
### Root
# root<-'N:/gebhyd/8_Him/Personal_folders/Pascal'
root<-'E:'

## Path to online storage
path_onst<-'C:/Users/Buri/switchdrive/WSL' 

## Path to stored data
path_data<-root&'/T&C/OUTPUTSTORAGE/MULTIPOINT/'&outnm

## Path to input data
path_input<-root&'/T&C/INPUTS'

## Path to store data
path_store<-root&'/T&C/POSTPROCESSED/MULTIPOINT'

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
library(scales)       #date_breaks()

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

# plotting settings
base_size<-21
base_family<-'Arial'

#===============================================================================
# ADD VARIABLE NAMES & CORRESPONDING UNITS
#===============================================================================
setwd(path_code)
setwd('..')
var_df<-read.table('TCmp_Variables.txt',header=T,sep='\t')
var_df[]<-lapply(var_df,as.character)


#===============================================================================
# READ SIMULATION FILE << AAA_TCmp_XXX.m >>
#===============================================================================
setwd(path_data)
files<-list.files()
##read AAA-file
simfile<-files[startsWith(files,'AAA_')]
simfile<-readLines(simfile,warn=FALSE)


#===============================================================================
# READ PARAMETERS FILE << PARAMETERS_SOIL_XXX.m >>
#===============================================================================
paramfile<-grep(pattern="^fnParameter",simfile)
paramfile<-simfile[paramfile]
paramfile<-str_split(paramfile,"'")[[1]][2]
paramfile<-readLines(path_input&'/'&paramfile,warn=FALSE)


#===============================================================================
# READ SOIL DISCRETIZATION
#===============================================================================
# find line where Zs (soil layer depths) is defined
rawfile<-paramfile
idx1<-grep(pattern="^Zs",x=rawfile)
rawfile<-rawfile[idx1[1]]
rawfile<-str_split(rawfile,'%')[[1]][1]
Zs<-as.numeric(str_extract_all(rawfile,"[0-9]+")[[1]])
rm(idx1,rawfile)

# find middle depth of layers
dz<-diff(Zs) ### [mm] Thickness of the Layers
Dz<-vector()
for(ii in 1:(length(Zs)-1)){
  ifelse(ii > 1,
         Dz[ii]<-(dz[ii]+ dz[ii-1])/2, ## Delta Depth Between Middle Layer[mm]
         Dz[ii]<-dz[1]/2               ## Delta Depth Between First Middle Layer and soil surface [mm]
  )
}
DZ<-cumsum(Dz)
rm(Dz)

Z_df<-cbind.data.frame(ly=as.character(1:(length(Zs)-1)),Zs_lb=Zs[-1],dz=dz,center=DZ)


#===============================================================================
# GET POI-INFORMATION
#===============================================================================
fn<-grep(pattern="^fnPOIs",simfile)
fn<-simfile[fn]
fn<-str_split(fn,"'")[[1]][2]

pois_df<-read.table(path_input&'/'&fn,header=T)
rm(fn)


#===============================================================================
# CREATE NEW FOLDERS
#===============================================================================
new_dir<-path_store&'/'&outnm
dir.create(new_dir)
dir.create(new_dir&'/PLOTS')
dir.create(new_dir&'/PLOTS/'&SITE)

setwd(path_data&'/OUTPUTS')
#===============================================================================
# CHHOSE SITE
#===============================================================================
fn_ls<-list.files(pattern='_results.mat')
fn_ls<-fn_ls[startsWith(fn_ls,SITE)]
n<-length(fn_ls)

tab_ls<-list()
for(i in 1:n){
  ##read table
  idx<-fn_ls[i]
  ID<-str_split(idx,'-')[[1]][1]
  idx<-which(pois_df$ID == ID)
  fn<-pois_df$ID[idx]&'-'&pois_df$ID2[idx]&'_'&pois_df$elevation[idx]&'m'
  print('+++ Reading table <<'&fn&'>> ('&i&'/'&n&') +++')
  fn<-new_dir&'/TABLES/'&fn&'___daily_avg.txt'
  tab<-read.table(fn,header=T,dec='.',sep='\t')
  tab$Date<-as.POSIXct(tab$Date,format='%Y-%m-%d')
  
  tab_ls[[i]]<-tab
  names(tab_ls)[[i]]<-ID
  rm(idx,ID,fn,tab)
}


#===============================================================================
# PLOT HEATMAPS FOR SINGLE VARIABLES (ROW PER FIELDSITE, COLUMN PER MONTH)
#===============================================================================

voi_v<-c('Csno','Ds','ea','EG','ESN','ET','f','G','H','NPP','O','Pr.liq','Pr.sno',
         'Rn','SND','T','Ta','Tdamp','V','Vice')
cols<-c("magma","inferno","plasma","viridis","cividis","rocket","mako","turbo")
cols<-c(cols,cols,cols,cols,cols,cols,cols)


for(i in 1:length(voi_v)){
  VOI<-voi_v[i]
  
  ##subset data
  if(VOI %in% c('O','V','Vice')){    #soil layers
    df<-lapply(tab_ls, function(x) x[, startsWith(names(x),VOI&'_ly')])
    sites<-names(tab_ls)
    df_tot<-list()
    rmidx<-vector()
    for(j in 1:length(df)){
      DF<-df[[j]]
      if(ncol(DF) == 0){
        df_tot[[j]]<-NULL
        rmidx<-c(rmidx,j)
      }
      if(ncol(DF) > 0){
        ly<-str_split(names(DF),'_ly',simplify=T)[,2]
        idx<-match(Z_df$ly,ly)
        idx<-stats::na.omit(idx)
        w<-Z_df$dz[idx]
        w<-w/sum(w)     #check: sum(w) should be 1!
        DF<-sweep(DF, 2, w, FUN="*")   #each layer multiplied with layer-specific weight
        DF<-rowSums(DF)                #and summed up
        df_tot[[j]]<-DF
      }
      rm(DF)
    }
    sites<-sites[-rmidx]
    rm(rmidx)
  }
  if(VOI %notin% c('O','V','Vice')){   
    df<-lapply(tab_ls, function(x) x[, names(x) == VOI])
    if(length(df[[1]]) == 0){
      df<-lapply(tab_ls, function(x) x[, names(x) == paste0(VOI,'.H')])  ##high vegetation
    }
    if(length(df[[1]]) == 0){
      df<-lapply(tab_ls, function(x) x[, names(x) == paste0(VOI,'.L')])  ##low vegetation
    }
  }
  
  ##combine sites & add date
  df<-as.data.frame(do.call('cbind',df))
  df<-cbind(Date=tab_ls[[1]][,'Date'],df)
  if(exists('df_tot')){
    df_tot<-delete.NULLs(df_tot)
    df_tot<-as.data.frame(do.call('cbind',df_tot))
    names(df_tot)<-sites
    df_tot<-cbind(Date=tab_ls[[1]][,'Date'],df_tot)
  }
  
  ##aggregate
  if(VOI %in% c('EG','ESN','ET','f','Pr.liq','Pr.sno','T')){
    df<-stats::aggregate(df[,-1],list(Date=cut(df$Date,'month')),sum) ##sum up from mm/h to mm/mo
    df[,-1]<-df[,-1]*24
  }
  if(VOI %notin% c('EG','ESN','ET','f','Pr.liq','Pr.sno','T')){
    df<-stats::aggregate(df,list(Date=cut(df$Date,'month')),mean)[,-1]         ##average from daily to monthly mean
    if(exists('df_tot')){
      df_tot<-stats::aggregate(df_tot,list(Date=cut(df_tot$Date,'month')),mean)[,-1] 
    }
  }
  
  ##add soil layer depths
  if(grepl('_ly',names(df)[2])){
    for(j in 1:nrow(Z_df)){
      ly<-Z_df$ly[j]
      ce<-Z_df$center[j]/10  #mm -> cm
      idx<-which(grepl('.'&VOI&'_ly'&ly,names(df)))
      names(df)[idx]<-str_replace(names(df)[idx],'.'&VOI&'_ly'&ly,'_'&ce&'cm')
      rm(idx)
    }
    idx<-match(names(df_tot),pois_df$ID)
    idx<-na.omit(idx)
    names(df_tot)[-1]<-names(df_tot)[-1]&'_'&pois_df$SOILTH[idx]*100&'cm'
  }
  
  df$Date<-year(df$Date)*100+month(df$Date)
  df<-melt(df,'Date')
  if(exists('df_tot')){
    df_tot<-melt(df_tot,'Date')
  }
  # ggplot(df)+
  #   geom_line(aes(x=Date,y=value,colour=variable))
  
  xlab<-''
  ylab<-''
  tit<-''
  if(VOI == 'ANPP'){leg<-'Above Ground Net Primary Production \n (monthly avg.) [gC m-2 day-1]'}
  if(VOI == 'Csno'){leg<-'Snow cover fraction \n(monthly avg.) [-]'}
  if(VOI == 'Ds'){leg<-'Vapour Pressure Deficit \n (monthly avg.) [Pa]'}
  if(VOI == 'ea'){leg<-'Vapour Pressure \n (monthly avg.) [Pa]'}
  if(VOI == 'EG'){leg<-'Bare Ground Evaporation\n (monthly sum) [mm]'}
  if(VOI == 'ESN'){leg<-'Snow sublimation \n (monthly sum) [mm]'}
  if(VOI == 'ET'){leg<-'Total Evapotranspiration \n (monthly sum) [mm]'}
  if(VOI == 'f'){leg<-'Infiltration \n (monthly avg.) [mm]'}
  if(VOI == 'G'){leg<-'Ground Heat Flux \n (monthly avg.) [W m-2]'}
  if(VOI == 'H'){leg<-'Sensible Heat Flux \n (monthly avg.) [W m-2]'}
  if(VOI == 'NPP'){leg<-'Net Primary Production \n (monthly avg.) [gC m-2 day-1]'}
  if(VOI == 'O'){leg<-'Soil Moisture \n (monthly avg.) [-]'}
  if(VOI == 'Pr.liq'){leg<-'Rainfall \n (monthly sum) [mm]'}
  if(VOI == 'Pr.sno'){leg<-'Snowfall \n (monthly sum) [mm]'}
  if(VOI == 'Rn'){leg<-'Net Radiation \n (monthly avg.) [W m-2]'}
  if(VOI == 'SND'){leg<-'Snow Depth \n (monthly avg.) [m]'}
  if(VOI == 'T'){leg<-'Transpiration \n (monthly sum) [mm]'}
  if(VOI == 'Ta'){leg<-'2m Air Temperature \n (monthly avg.) [°C]'}
  if(VOI == 'Tdamp'){leg<-'Soil Temperature at dampening depth \n (monthly avg.) [°C]'}
  if(VOI == 'V'){leg<-'Liquid water content \n (monthly avg.) [mm]'}
  if(VOI == 'Vice'){leg<-'Frozen water content \n (monthly avg.) [mm]'}
  
  capt<-''
  brk<-c('201601','201701','201801','201901','202001','202101')
  lab<- c('2016',  '2017',  '2018',  '2019',  '2020',  '2021')
  
  ##drop non-soil sites for plotting specific variables that are zero for non-soil sites
  if(VOI %in% c('f','EG')){
    rmidx<-pois_df$ID[which(pois_df$SOILTH == 0)]
    rmidx<-which(df$variable %in% rmidx)
    df<-df[-rmidx,]
    rm(rmidx)
  }
  
  
  if(VOI %in% c('ESN','G','H','Ta','Tdamp')){
    zli<-mroundu(max(abs(range(df$value))),1)
    pl<-ggplot(df,aes(x=as.character(Date),y=variable))+ 
      geom_tile(aes(fill=value))+ 
      scale_fill_stepsn(colors=c('#2166ac','#67a9cf','#d1e5f0','#f7f7f7',
                                 '#fddbc7','#ef8a62','#b2182b'),
                        n.breaks=8,limits=c(-zli,zli),show.limits=F)+
      scale_x_discrete(breaks=brk,labels=lab)+
      # scale_x_date(breaks=date_breaks("1 year"),labels=date_format("%Y"),
      #              expand=c(0,0))+
      labs(x=xlab,y=ylab,fill=leg,caption=capt,tile=tit)+ 
      theme_gray(base_size=base_size,base_family=base_family)+
      theme(legend.title = element_text(size = base_size), 
            legend.position = "top",
            legend.key.height= unit(0.5, 'cm'),
            legend.key.width= unit(2.2, 'cm'),
            panel.background = element_blank(), 
            axis.text.x = element_text(colour = "black", size = base_size, angle = 90),
            axis.text.y = element_text(colour = "black", size = base_size), 
            axis.title = element_text(size = base_size+1), 
            legend.text = element_text(size = base_size), legend.key = element_blank())
    
    pl
  }
  
  if(VOI %in% c('Csno','Ds','ea','EG','ET','f','NPP','Pr.liq','Pr.sno',
                       'Rn','SND','T')){
    pl<-ggplot(df,aes(x=as.character(Date),y=variable))+ 
      geom_tile(aes(fill=value))+ 
      scale_fill_viridis(option=cols[i],direction=1)+
      scale_x_discrete(breaks=brk,labels=lab)+
      # scale_x_date(breaks=date_breaks("1 year"),labels=date_format("%Y"),
      #              expand=c(0,0))+
      labs(x=xlab,y=ylab,fill=leg,caption=capt,tile=tit)+ 
      theme_gray(base_size=base_size,base_family=base_family)+
      theme(legend.title = element_text(size = base_size), 
            legend.position = "top",
            legend.key.height= unit(0.5, 'cm'),
            legend.key.width= unit(2.2, 'cm'),
            panel.background = element_blank(), 
            axis.text.x = element_text(colour = "black", size = base_size, angle = 90),
            axis.text.y = element_text(colour = "black", size = base_size), 
            axis.title = element_text(size = base_size+1), 
            legend.text = element_text(size = base_size), legend.key = element_blank())
    
    pl
  }
  
  if(VOI %in% c('O','V','Vice')){
    ## total soil profile
    pl<-ggplot(df_tot,aes(x=as.character(Date),y=variable))+ 
      geom_tile(aes(fill=value))+ 
      scale_fill_viridis(option=cols[i],direction=1)+
      scale_x_discrete(breaks=brk,labels=lab)+
      # scale_x_date(breaks=date_breaks("1 year"),labels=date_format("%Y"),
      #              expand=c(0,0))+
      labs(x=xlab,y=ylab,fill=leg,caption=capt,tile=tit)+ 
      theme_gray(base_size=base_size,base_family=base_family)+
      theme(legend.title = element_text(size = base_size), 
            legend.position = "top",
            legend.key.height= unit(0.5, 'cm'),
            legend.key.width= unit(2.2, 'cm'),
            panel.background = element_blank(), 
            axis.text.x = element_text(colour = "black", size = base_size, angle = 90),
            axis.text.y = element_text(colour = "black", size = base_size), 
            axis.title = element_text(size = base_size+1), 
            legend.text = element_text(size = base_size), legend.key = element_blank())
    
    pl
    
    rm(df_tot)
    
    
    ## individual soil layers
    DEPTHS<-c(0.5,12.5,45,85)# cm (have to correspond to existing layer center depths)
    for(j in 1:length(DEPTHS)){
      DEPTH<-as.character(DEPTHS[j])
      DF<-df[which(grepl(DEPTH&'cm',df$variable)),]
      DF$variable<-str_replace(DF$variable,'_'&DEPTH&'cm','')
      LEG<-str_replace(leg,'\n','at '&DEPTH&' cm \n')
      pl_ly<-ggplot(DF,aes(x=as.character(Date),y=variable))+ 
        geom_tile(aes(fill=value))+ 
        scale_fill_viridis(option=cols[i],direction=1)+
        scale_x_discrete(breaks=brk,labels=lab)+
        # scale_x_date(breaks=date_breaks("1 year"),labels=date_format("%Y"),
        #              expand=c(0,0))+
        labs(x=xlab,y=ylab,fill=LEG,caption=capt,tile=tit)+ 
        theme_gray(base_size=base_size,base_family=base_family)+
        theme(legend.title = element_text(size = base_size), 
              legend.position = "top",
              legend.key.height= unit(0.5, 'cm'),
              legend.key.width= unit(2.2, 'cm'),
              panel.background = element_blank(), 
              axis.text.x = element_text(colour = "black", size = base_size, angle = 90),
              axis.text.y = element_text(colour = "black", size = base_size), 
              axis.title = element_text(size = base_size+1), 
              legend.text = element_text(size = base_size), legend.key = element_blank())
      
      pl_ly
      
      fn<-new_dir&'/PLOTS/'&SITE&'/'&VOI&'-heatmap_ps-pmo_at'&DEPTH&'cm.png'
      ggsave(fn,pl_ly,width=25,height=15,units='cm')  ## save plot
      graphics.off()
      
      rm(pl_ly,DF,LEG,fn)
    }
  }
  
  fn<-new_dir&'/PLOTS/'&SITE&'/'&VOI&'-heatmap_ps-pmo.png'
  ggsave(fn,pl,width=25,height=15,units='cm')  ## save plot
  graphics.off()
  
  rm(pl,fn)
}





# ##### SNOW COVER PER SITE
# df<-lapply(tab_ls, function(x) x[, names(x) == 'Csno'] )
# df<-as.data.frame(do.call('cbind',df))
# df$Date<-tab_ls[[1]][,'Date']
# 
# df<-melt(df,'Date')
# ggplot(df)+
#   geom_line(aes(x=Date,y=value,colour=variable))
# 
# 
# ##### TA PER SITE
# df<-lapply(tab_ls, function(x) x[, names(x) == 'Ta'] )
# df<-as.data.frame(do.call('cbind',df))
# df$Date<-tab_ls[[1]][,'Date']
# 
# df<-melt(df,'Date')
# ggplot(df)+
#   geom_line(aes(x=Date,y=value,colour=variable))
# 
# 
# ##### NPP.L PER SITE
# df<-lapply(tab_ls, function(x) x[, names(x) == 'NPP.L'] )
# df<-as.data.frame(do.call('cbind',df))
# df$Date<-tab_ls[[1]][,'Date']
# 
# df<-melt(df,'Date')
# ggplot(df)+
#   geom_line(aes(x=Date,y=value,colour=variable))

