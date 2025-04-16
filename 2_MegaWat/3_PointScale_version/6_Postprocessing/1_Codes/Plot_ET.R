################################################################################
# POSTPROCESSING PLOTS
#
#===============================================================================
#
# 2025/01/06
#
#
# Maximiliano Rodriguez | Cryosphere and Mountain Hydrosphere | 
# Institute of Science and Technology (ISTA) |
# Am Campus 1 | 3400 Klosterneuburg | mrodrigu@ist.ac.at
#
# Code explanation:
#   Plot for Evaporation from MODIS images
#
# Units
# MODIS product
#   ET: Kg m-2 8days-1
#
# TC model
#   ET: mm day-1
#
#
################################################################################


#=============================================================================== 
# CLEAN
#=============================================================================== 

rm(list = ls())
#dev.off()

#=============================================================================== 
# POINT
#===============================================================================

site = 'Lago_di_Corbara' 

#=============================================================================== 
# Libraries
#=============================================================================== 

library(readxl) #for opening excel files
library(dplyr) # for data manipulation
library(lubridate) #for dates
library(extrafont) #Fonts
library('SPEI') # Calculate evapotranspiration

#=============================================================================== 
# Definition of characters
#=============================================================================== 
# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

#=============================================================================== 
# PATHS
#=============================================================================== 

# TC results
path_model = 'C:/Users/mrodrigu/Desktop/19_ISTA/1_TC/3_Model_Source/2_MegaWat/'
path_results =  path_model&'3_PointScale_version/4_Outputs/'&site&'_daily_results.txt'

# MODIS product
path_MODIS = 'C:/Users/mrodrigu/Desktop/19_ISTA/14_Evaluation/2_Apennines/1_MOD16A2_061_2008/4_Series_unified/'
path_data =  path_MODIS&'MOD16A2_Monte.xlsx'

#=============================================================================== 
# RESULTS
#===============================================================================

TC_data <- read.csv(path_results, header = TRUE)
str(TC_data) #Type of data

# Add dates
TC_data <- TC_data %>%
  mutate(Julian_day = row_number(), .after = "Date")

# Average each 8 days

ET_data <- TC_data %>%
  mutate(Julian = ceiling(row_number() / 8)) %>%
  mutate(Julian = Julian*8-7) %>%
  group_by(Julian) %>%
  summarize(ET_8days = sum(ET))

#MODIS DATA
MODIS_data <- read_excel(path_data)


#=============================================================================== 
# CALCULATIONS FOR CHECK
#===============================================================================

#Annual evaporation MODIS
# Kg m-2 y-1 = mm y-1
ET_year = sum(MODIS_data$Monte_ET_500m)

#=============================================================================== 
# CALCULATIONS theoretical calculation
#===============================================================================

path_stations <- 'C:/Users/mrodrigu/Desktop/19_ISTA/7_Forcing/1_Raw_data/2_Apennines/2_Stations/1_CNR/9_Temperature/5_Table/'
load(path_stations&'Temp_min.RData')
Temp_min <- DAT
load(path_stations&'Temp_max.RData')
Temp_max <- DAT
rm(DAT)   

#Temp_min <- Temp_min %>% select(1, order(colnames(.)[-1]) + 1)
Orv_scalo_min <- Temp_min %>% select(Date, Orvieto_Scalo)
Orv_scalo_max <- Temp_max %>% select(Date, Orvieto_Scalo)

Orv_scalo_min_monthly <- Orv_scalo_min %>%
                         mutate(month = month(Date)) %>% # Create a month column
                         group_by(month) %>%
                         summarize(temp_min = mean(Orvieto_Scalo, na.rm = TRUE)) # Calculate monthly average

Orv_scalo_max_monthly <- Orv_scalo_max %>%
                         mutate(month = month(Date)) %>% # Create a month column
                         group_by(month) %>%
                         summarize(temp_max = mean(Orvieto_Scalo, na.rm = TRUE)) # Calculate monthly average


# PET according to Thornthwaite
#tho <- thornthwaite(TC_data$Ta,42.7)

# Hargreaves
har <- hargreaves(Orv_scalo_min_monthly$temp_min, Orv_scalo_max_monthly$temp_max,lat=42.7)

# Penman, based on sun hours, ignore NAs
#pen <- penman(TMIN,TMAX,AWND,tsun=TSUN,lat=37.6475,z=402.6,na.rm=TRUE)

#=============================================================================== 
# PLOTS
#===============================================================================

path_plot = "C:/Users/mrodrigu/Desktop/19_ISTA/1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/6_Postprocessing/2_Plots/"
tiff(filename = path_plot&"ET.tif", width = 5, height = 5, units = "in", res = 600 ,pointsize = 12)

#x = 1:nrow(TC_data_sel)
yl = c(0,40) #y axis size

plot(ET_data$Julian,ET_data$ET_8days,xlim = c(0,365),ylim = yl, type="o", pch=16,
     xlab="",ylab="", frame.plot = F,axes=F, col = 'blue', xaxs = "i", yaxs = "i")
par(new = T)
plot(MODIS_data$Day,MODIS_data$Monte_ET_500m,xlim = c(0,365),ylim = yl, type="o", pch=16,
     xlab="",ylab="", frame.plot = F,axes=F, col = 'red', xaxs = "i", yaxs = "i")

#lb_pos = seq(1,nrow(ET_data),8)

axis(1,las=2, at=  ET_data$Julian,  ##seq(1,nrow(xx2),1), 
     labels = ET_data$Julian, cex.axis = 0.7) #c(xx2[,"Interval"]), cex.axis = 1.5)

axis(2,las=1, cex.axis = 1) 

mtext("Julian calendar", side = 1, line = 2,cex = 1.0)
mtext("ET [mm 8day-1]", side = 2, line=2,cex = 1.0, adj = 0.5)
mtext("Monte Terminillo (Lazio)", side = 3, line=1,cex = 1.0)

legend("topright",legend = c("TC-modelled","MODIS product"),
       col = c('blue','red'), title = "", ncol = 1,
       lty = c(1,1), bty = "n", cex = 1)

dev.off()

