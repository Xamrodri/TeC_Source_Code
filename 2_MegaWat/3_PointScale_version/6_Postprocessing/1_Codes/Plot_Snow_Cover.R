#===============================================================================
# POSTPROCESSING PLOTS
#
#===============================================================================
#
# 2025/01/06
#
#
#
# Maximiliano Rodriguez | High Mountain Glaciers and Hydrology | 
#  Swiss Federal Institute for Forest, Snow and Landscape Research, WSL |
#  Office MG E 35 | Zürcherstrasse 111, 8903 Birmensdorf | pascal.buri@wsl.ch
#
################################################################################

rm(list = ls())
dev.off()

# define &-sign for pasting string-elements
'&' <- function(...) UseMethod('&')
'&.default' <- .Primitive('&')
'&.character' <- function(...) paste(...,sep='')

path_model = 'M:/19_ISTA/1_TC/3_Model_Source/2_MegaWat/'
path_results =  path_model&'3_PointScale_version/4_Outputs/Apennine_hill_results.txt'

my_data <- read.csv(path_results, header = TRUE)

x = 1:nrow(my_data)

plot(x,my_data[,'SND'],xlim = c(0,100),ylim = c(0,10), type="o", pch=16,
     xlab="",ylab="", frame.plot = F,axes=F, col = 'blue', xaxs = "i", yaxs = "i")

axis(1,las=2, at= seq(1,nrow(my_data),1),  ##seq(1,nrow(xx2),1), 
     labels = my_data[,"Date"]) #c(xx2[,"Interval"]), cex.axis = 1.5)

axis(2,las=1, cex.axis = 1) 
mtext("Snow depth [m]", side = 2, line=2,cex = 1.0, adj = 0.5)


