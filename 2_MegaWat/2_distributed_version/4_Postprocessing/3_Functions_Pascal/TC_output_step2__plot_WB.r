################################################################################
# T&C OUTPUT:
#  PLOT WATER BALANCe PARTITION
#
# 
# NEW: -
#
# TODO: -
#
# 2024/01/04
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

## define timespan & months to analyze
ts_init<-'2016-10-01 00:00:00'
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
msk<-raster('msk.tif')
ele_rng<-range(values(dem),na.rm=T)
resol<-res(dem)[1]


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
dates<-cbind.data.frame(Date=df$Date,count=1)
#derive monthly mean
sumidx<-which(names(df) %in% c('Pr','Pr_liq','Pr_sno','T'))
df_avg<-aggregate(df[,-c(1,sumidx)],list(Date=cut(df$Date,'month')),mean) 
df_avg$Date<-as.POSIXct(df_avg$Date,format='%Y-%m-%d')
#derive monthly sum
df_sum<-aggregate(df[,sumidx],list(Date=cut(df$Date,'month')),sum) 
df<-cbind.data.frame(df_avg,df_sum[,-1])
rm(df_avg,df_sum)
df<-df[which(df[,1] >= ts_init & df[,1] <= ts_final),] #crop to selected period


#===============================================================================
# ADD VARIABLE NAMES & CORRESPONDING UNITS
#===============================================================================
setwd(path_code)
setwd('..')
var_df<-read.table('TC_Variables.txt',header=T,sep='\t')
var_df[]<-lapply(var_df,as.character)


newdir<-path_data&'/PLOTS'
dir.create(newdir)
dir.create(newdir&'/WB')


#===============================================================================
#===============================================================================
# WATER BALANCE PLOTTING
#===============================================================================
#===============================================================================
nms_pos<-c('Smelt','Imelt','Pr.liq','Pr.sno')
nms_neg<-c('EG','T','EIn','ESN','EWAT','EIn.rock','EICE')
nms_sto<-c('V.clnd','Vice.clnd','SWE')

idx_pos<-match(nms_pos,var_df$var)
idx_neg<-match(nms_neg,var_df$var) 
idx_sto<-match(nms_sto,var_df$var) 
var_pos<-var_df$variable[idx_pos]
var_neg<-var_df$variable[idx_neg]
var_sto<-var_df$variable[idx_sto]
units_pos<-var_df$unit[idx_pos]
units_neg<-var_df$unit[idx_neg]
units_sto<-var_df$unit[idx_sto]
rm(idx_pos,idx_neg,idx_sto)

nms_pos<-str_replace(nms_pos,'\\.','_')
nms_neg<-str_replace(nms_neg,'\\.','_')
nms_sto<-str_replace(nms_sto,'\\.','_')
df_pos<-df[,match(c('Date',nms_pos),names(df))] #use match() to keep order
names(df_pos)[-1]<-var_pos
df_neg<-df[,match(c('Date',nms_neg),names(df))] 
names(df_neg)[-1]<-var_neg
df_sto<-df[,match(c('Date',nms_sto),names(df))] 
names(df_sto)[-1]<-var_sto

##add basin runoff ('QpointC')
setwd(path_data&'/TABLES')
fn<-list.files(pattern='TrackedPixel___Outlet-')
Q_ls<-list()
for(i in 1:length(fn)){
  Q<-read.table(fn[i],header=T,sep='\t')
  Q$Date<-as.POSIXct(Q$Date,format='%Y-%m-%d %H:%M:%S')
  Q<-Q[,which(names(Q) %in% c('Date','QpointC'))]     ##Channel discharge [mm]
  # Q[,2]<-Q[,2]*resol^2/(3600*1000)##[m3/s]
  ifelse(i == 1,
         Q_ls[[i]]<-Q,
         Q_ls[[i]]<-Q[,2]
  )
  rm(Q)
}
df_Q<-as.data.frame(do.call('cbind',Q_ls))
if(ncol(df_Q) > 2){df_Q<-cbind.data.frame(Date=df_Q[,1],rowSums(df_Q[,-1]))}

################################
##rescaling: Precip.?, Runoff?
# df_pos$Precipitation<-df_pos$Precipitation/sum(values(msk))
df_Q[,2]<-df_Q[,2]/sum(values(msk))
################################

##derive monthly Q values
df_Q<-aggregate(df_Q[,-1],list(Date=cut(df_Q$Date,'month')),sum)
df_Q$Date<-as.POSIXct(df_Q$Date,format='%Y-%m-%d')
df_Q<-df_Q[which(df_Q[,1] >= ts_init & df_Q[,1] <= ts_final),] #crop to selected period


# ###check variables
# VAR<-'Imelt'
# VAR<-'Smelt'
# VAR<-'EIn_rock'
# VAR<-'EICE'
# VAR<-'ESN'
# VAR<-'EWAT'
# VAR<-'EIn_H'
# VAR<-'EIn_L'
# VAR<-'EG'
# setwd(path_data&'/RASTERS/MONTHLY')
# tt<-brick(list.files(pattern='^'&VAR&'___')[1])
# levelplot(tt[[6:11]]*24*30.25)
# levelplot(tt[[6:11]])
# rm(VAR,tt)


##adjust values for hours per month 
hpm<-aggregate(dates[,'count'],list(Date=cut(dates$Date,'month')),sum) 
names(hpm)<-c('Date','count')
hpm$Date<-as.POSIXct(hpm$Date,format='%Y-%m-%d')
hpm<-hpm[which(hpm[,1] >= ts_init & hpm[,1] <= ts_final),] #crop to selected period

df_pos$`Snow melt`<-df_pos$`Snow melt`*hpm$count
df_pos$`Ice melt`<-df_pos$`Ice melt`*hpm$count

df_neg$`Evapor. from Bare soil`<-df_neg$`Evapor. from Bare soil`*hpm$count
df_neg$`Evapor. from interc. water Veg.`<-df_neg$`Evapor. from interc. water Veg.`*hpm$count
# df_neg$`Transpir. Veg.`<-df_neg$`Transpir. Veg.`*hpm$count
df_neg$`Evapor. from the snowpack at the ground`<-df_neg$`Evapor. from the snowpack at the ground`*hpm$count
df_neg$`Evapor. from water and ponds`<-df_neg$`Evapor. from water and ponds`*hpm$count
df_neg$`Evapor. from rocks`<-df_neg$`Evapor. from rocks`*hpm$count
df_neg$`Evapor./sublimation from Ice`<-df_neg$`Evapor./sublimation from Ice`*hpm$count


#===============================================================================
# PLOT TOTAL ABSOLUTE VALUES PER HYDROLOGICAL YEAR:
#   ET(NON-FROZEN), IMELT, RAIN, Q, SMELT, SWE-CHANGE, SOILWATER-CHANGE
#   + 
#   (SNOWFALL, ESNOW/ICE, EG, T, EIn)
#===============================================================================
yrs<-unique(year(df$Date))[-1]  ##hydrological year
#derive annual totals
pl_ls<-list()
n<-length(yrs)
for(i in 1:n){
  #select timesteps within hydrological year
  idx1<-which(year(df_pos$Date) == (yrs[i]-1) & month(df_pos$Date) %in% 10:12) #Oct-Dec
  idx2<-which(year(df_pos$Date) == (yrs[i]) & month(df_pos$Date) %in% 1:9)     #Jan-Sep
  idx<-c(idx1,idx2); rm(idx1,idx2)
  
  v<-colSums(df_pos[idx,-1])
  v<-v[-4]##remove snowfall
  v<-c(v,'ET non-frozen'=sum(colSums(cbind(df_neg$`Evapor. from Bare soil`[idx],
                                           df_neg$`Transpir. Veg.`[idx],
                                           df_neg$`Evapor. from interc. water Veg.`[idx],
                                           df_neg$`Evapor. from water and ponds`[idx],
                                           df_neg$`Evapor. from rocks`[idx]))),
       Runoff=sum(df_Q[idx,-1]))
  soilstor<-(df_sto[max(idx),2]-df_sto[min(idx),2])+(df_sto[max(idx),3]-df_sto[min(idx),3])
  names(soilstor)<-'Soil storage change'
  snostor<-(df_sto[max(idx),4]-df_sto[min(idx),4])
  names(snostor)<-'Snow storage change'
  v<-c(v,soilstor,snostor)
  idx_neg<-match(c('ET non-frozen','Runoff'),names(v))
  v[idx_neg]<-v[idx_neg]*-1
  # >   names(v)
  # [1] "Snow melt"           "Ice melt"            "Liquid Precip."     
  # [4] "ET non-frozen"       "Runoff"              "Soil storage change"
  # [7] "Snow storage change"
  names(v)<-c("Snow melt","Ice melt","Rain",
              "ET non-frozen","Runoff","dSoil storage",
              "dSnow storage")
  inp<-round(sum(c(v['Snow melt'],v['Ice melt'],v['Rain'])),0)
  out<-round(sum(c(v['ET non-frozen'],v['Runoff'])),0)*-1
  dsoil<-round(v['dSoil storage'],0)
  dsnow<-round(v['dSnow storage'],0)
  
  #################################
  ##add specific variables 
  v<-c(v,'Snowfall'=sum(df_pos$`Solid (snow) Precip.`[idx]))
  v<-c(v,'Transp'=sum(df_neg$`Transpir. Veg.`[idx]))
  v<-c(v,'E bare'=sum(df_neg$`Evapor. from Bare soil`[idx]))
  v<-c(v,'E interc'=sum(df_neg$`Evapor. from interc. water Veg.`[idx]))
  v<-c(v,'E snow'=sum(df_neg$`Evapor. from the snowpack at the ground`[idx]))
  v<-c(v,'E ice'=sum(df_neg$`Evapor./sublimation from Ice`[idx]))
  v<-c(v,'E water'=sum(df_neg$`Evapor. from water and ponds`[idx]))
  v<-c(v,'E rock'=sum(df_neg$`Evapor. from rocks`[idx]))
  #################################
  
  #reshape
  v<-as.data.frame(v)
  v$x<-row.names(v)
  v$v<-round(v$v,0)
  
  # show_col(pal_nejm("default")(7))
  cols<-pal_nejm("default")(7)
  v_ord<-v
  v_ord$lab<-c('004_SMELT','003_IMELT','002_PLIQ','001_ETLIQ','007_Q','005_dSW','006_dSWE',
               '008_SNOWFALL','009_TRANSP','010_EBARE','011_EIN','012_ESN','013_EICE',
               '014_EWAT','015_EROCK')
  v_ord$lab<-factor(v_ord$lab, levels = v_ord$lab)
  v_ord$x<-factor(v_ord$x, levels = v_ord$x)
  cols<-c(cols[c(4,3,2,1,7,5,6)],rep('#636363',8))
  
  #define labels
  tit<-'hydrol. year '&yrs[i]
  capt<-''
  xlab<-''
  ylab<-'[mm]'
  leg<-''
  capt<-'Input (snow/ice melt, rain; '&inp&' mm) - Output (ET non-frozen, runoff; '&out&' mm) '&
    '+ dSoil ('&dsoil&' mm) + dSnow ('&dsnow&' mm) = '&inp-out+dsoil+dsnow&' mm'
  #plotting
  if(i < n){
  pl<-ggplot(data=v_ord,aes(x=lab,y=v))+
    geom_bar(stat='identity',colour="#636363",fill=cols)+
    geom_text(aes(label=v),position=position_dodge(width=0.9),vjust=-0.25,size=5)+
    labs(title=tit,caption=capt,x=xlab,y=ylab,fill=leg)+
    scale_x_discrete(labels=v_ord$x)+
    geom_hline(yintercept=0,linewidth=1,colour='black')+
    theme_gray(base_size=base_size,base_family=base_family)+
    theme(
      legend.position = 'none',
      # legend.direction = 'vertical',
      # legend.title = element_blank(),
      # legend.text = element_text(size=base_size),
      # legend.key=element_rect(fill='white'),
      axis.ticks = element_blank(),
      axis.text.y = element_text(colour='black',size=base_size),
      axis.text.x = element_blank(),
      panel.background = element_rect(fill='white',colour='white'),
      panel.grid.minor.y = element_line(linewidth=0.1,colour='grey'),
      panel.grid.major.y = element_line(linewidth=0.75,colour='grey'),
      panel.grid.major.x = element_line(linewidth=0.75,colour='grey'),
      strip.background = element_blank(),
      strip.text.y = element_blank(),
      plot.title = element_text(colour='black',size=base_size-2),
      plot.margin=unit(c(0,0,0,0),'lines') #top/right/bottom/left
    ) 
  }
  if(i == n){
    pl<-ggplot(data=v_ord,aes(x=lab,y=v))+
      geom_bar(stat='identity',colour="#636363",fill=cols)+
      geom_text(aes(label=v),position=position_dodge(width=0.9),vjust=-0.25,size=5)+
      labs(title=tit,caption=capt,x=xlab,y=ylab,fill=leg)+
      scale_x_discrete(labels=v_ord$x)+
      geom_hline(yintercept=0,linewidth=1,colour='black')+
      theme_gray(base_size=base_size,base_family=base_family)+
      theme(
        legend.position = 'none',
        # legend.direction = 'vertical',
        # legend.title = element_blank(),
        # legend.text = element_text(size=base_size),
        # legend.key=element_rect(fill='white'),
        axis.ticks = element_blank(),
        axis.text.y = element_text(colour='black',size=base_size),
        axis.text.x = element_text(colour='black',size=base_size,angle=90),
        panel.background = element_rect(fill='white',colour='white'),
        panel.grid.minor.y = element_line(linewidth=0.1,colour='grey'),
        panel.grid.major.y = element_line(linewidth=0.75,colour='grey'),
        panel.grid.major.x = element_line(linewidth=0.75,colour='grey'),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        plot.title = element_text(colour='black',size=base_size),
        plot.margin=unit(c(0,0,0,0),'lines') #top/right/bottom/left
      ) 
  }
  print(pl)
  Sys.sleep(1)
  pl_ls[[i]]<-pl
  rm(pl)
}

pl<-marrangeGrob(pl_ls,nrow=length(yrs),ncol=1,
                 top=grid::textGrob(DOMAIN,gp=grid::gpar(fontsize=base_size+2)))

# grid.arrange(pl_ls, ncol=1,top=grid::textGrob(DOMAIN,gp=grid::gpar(fontsize=base_size)))

pl

fn<-newdir&'/WB/WB_annualtotals_'&PERIOD_SH&'.png'
ggsave(fn,pl,width=30,height=60,units='cm')
graphics.off()
print(fn)










