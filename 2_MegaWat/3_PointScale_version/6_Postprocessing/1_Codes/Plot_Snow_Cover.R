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
#   Plot for Snow depth 
################################################################################

#=============================================================================== 
# CLEAN
#=============================================================================== 

rm(list = ls())
#dev.off()

#=============================================================================== 
# Libraries
#=============================================================================== 

library(readxl) #for opening excel files
library(dplyr) # for data manipulation
library(lubridate) #for dates
library(extrafont) #Fonts
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
path_results =  path_model&'3_PointScale_version/4_Outputs/Monte_Terminillo_daily_results.txt'

TC_data <- read.csv(path_results, header = TRUE)
str(TC_data) #Type of data

# Observed data
path_station =  path_model&'3_PointScale_version/7_Evaluation_data/Snow_Depth.xlsx'

obs_data <- read_excel(path_station)
str(obs_data) #Type of data
#=============================================================================== 
# DATES FOR PLOTS
#===============================================================================

# Plotting period
x1 <- as.POSIXct("2008-01-01 00:00:00", format="%Y-%m-%d %H:%M:%S", tz="Europe/Rome")
x2 <- as.POSIXct("2008-12-30 23:00:00", format="%Y-%m-%d %H:%M:%S", tz="Europe/Rome")

#Valid time zones
#valid_timezones <- OlsonNames()

# Checking time zone of obs_data dataframe
current_tz <- attr(obs_data$Date, "tzone")

# Changing time zone to Europe
obs_data$Date <- with_tz(obs_data$Date, tzone = "Europe/Rome")

# Changing type of TC_data date column.
TC_data$Date <- as.POSIXct(TC_data$Date, format="%d-%b-%Y", tz="Europe/Rome")

#=============================================================================== 
# Change of units for consistency
#===============================================================================
obs_data[,'Monte_Terminillo_m'] = obs_data[,'Monte_Terminillo']/100

#=============================================================================== 
# Data selection for plots
#===============================================================================

obs_data_sel <- obs_data %>% filter(Date >= x1 & Date <= x2) 
TC_data_sel <- TC_data %>% filter(Date >= x1 & Date <= x2)

#=============================================================================== 
# PLOTS
#===============================================================================

path_plot = "C:/Users/mrodrigu/Desktop/19_ISTA/1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/6_Postprocessing/2_Plots/"
tiff(filename = path_plot&"Snow_Depth.tif", width = 5, height = 5, units = "in", res = 600 ,pointsize = 12)

x = 1:nrow(TC_data_sel)
yl = c(0,2) #y axis size

plot(x,TC_data_sel$SND,xlim = c(0,365),ylim = yl, type="o", pch=16,
     xlab="",ylab="", frame.plot = F,axes=F, col = 'blue', xaxs = "i", yaxs = "i")
par(new = T)
plot(x,obs_data_sel$Monte_Terminillo_m,xlim = c(0,365),ylim = yl, type="o", pch=16,
     xlab="",ylab="", frame.plot = F,axes=F, col = 'red', xaxs = "i", yaxs = "i")

lb_pos = seq(1,nrow(TC_data_sel),30)

axis(1,las=2, at= lb_pos,  ##seq(1,nrow(xx2),1), 
     labels = TC_data_sel[lb_pos,"Date"], cex.axis = 0.7) #c(xx2[,"Interval"]), cex.axis = 1.5)

axis(2,las=1, cex.axis = 1) 
mtext("Snow depth [m]", side = 2, line=2,cex = 1.0, adj = 0.5)
mtext("Monte Terminillo (Lazio)", side = 3, line=1,cex = 1.0)

legend("topright",legend = c("TC-modelled","Observed"),
       col = c('blue','red'), title = "", ncol = 1,
       lty = c(1,1), bty = "n", cex = 1)

dev.off()
