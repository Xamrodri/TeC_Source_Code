###############################################################################
# PARTITION LIQUID-SOLID PRECIPITATION (DING ET AL., 2014)
#
# Ding et al., 2014 JH  (http://dx.doi.org/10.1016/j.jhydrol.2014.03.038)
#
# Ta:     vector of air temperature [°C]
# rH:     vector of air relative humidity [%]
# Pres:   vector of atmospheric pressure [Pa]

PPARTITION_DING2014<-function(Ta,rH,Pres) {
  
  ## Constants
  g<-9.81 # [m/s^2] gravity acceleration
  P_Ref<-1013.25 # [Pa] reference pressure
  Rd<-287.05 # [J/kg K] dry air gas constant
  
  esat<-611*exp(17.27*Ta/(237.3+Ta)) # [Pa] Vapor pressure saturation
  ea<-esat*rH/100                        #Vapour pressure (Pa) 
  RH<-ea/esat # Relative Humidity
  Laten<-1000*(2501.3 - 2.361*(Ta)) # Latent heat vaporization/condensaition [J/kg]
  cp<-1005 + ((Ta +23.15)^2)/3364 # specific heat air  [J/kg K]
  gam<-cp*100*Pres/(0.622*Laten) # [Pa/C] psycrometric constant
  del<-(4098*esat)/((237.3+Ta)^2) # Pa/C
  Twb<-AWS$Ta - ( esat - ea )/( gam + del) # [C] Wet bulb temperature
  
  ## Reference elevation
  Zref<- -((Ta+15)/2+273.15)*(Rd/g)*log(Pres/P_Ref) # [m]
  Zref<-Zref/1000 # elevation, [m] => [km]
  
  ## Thresholds for discrimination precipitation
  dT<-0.215 - 0.099*RH + 1.018*RH^2
  dS<-2.374 - 1.634*RH
  T0<- -5.87 - 0.1042*Zref + 0.0885*Zref^2 + 16.06*RH - 9.614*RH^2
  Tmin<-T0
  Tmax<-T0
  
  idx<-dT/dS > log(2)
  Tmin[idx]<-T0[idx] - 
    dS[idx]*log(exp(dT[idx]/dS[idx]) - 
                  2*exp(-dT[idx]/dS[idx]))
  Tmax[idx]<-2*T0[idx] - Tmin[idx]
  rm(idx)
  
  ## Calculation of solid fraction of precipitation
  solid_fraction<-1/(1 + exp((Twb - T0)/dS)) # sleet
  
  Pr_sno<-vector()
  Pr_liq<-vector()
  for(i in 1:length(Twb)){
    # Discrimination of rain, sleet, and snow
    if(Twb[i] > Tmin[i] && Twb[i] < Tmax[i]){
      # sleet
      Pr_sno[i]<-AWS[i,'Rain_TB']*solid_fraction[i];
      Pr_liq[i]<-AWS[i,'Rain_TB']*(1-solid_fraction[i]);
    }
    if(Twb[i] <= Tmin[i]){
      # snow
      Pr_sno[i]<-AWS[i,'Rain_TB']
      Pr_liq[i]<- 0
    }
    if(Twb[i]>= Tmax[i]){
      # rain
      Pr_sno[i]<-0
      Pr_liq[i]<-AWS[i,'Rain_TB']
    }
  }
  
  return(cbind.data.frame(Pr_liq,Pr_sno))
}