#===============================================================================
# CORRECT FOR WIND-UNDERCATCH (MASUDA ET AL., 2019)
#===============================================================================
## Masuda et al., 2019 ESS (http://dx.doi.org/10.1029/2019EA000659)
## Tab. 2:
##  Unshielded RT-1 (Tipping Bucket) correction coefficient 
##  for CR (catch ratio) and m (correction parameter) 
##  (Yokoyama et al., 2003)
cc_s<-0.213  #correction coefficient snow 
cc_r<-0.0454 #correction coefficient rain

Ztb<-2     #tipping bucket height [m]
Zws<-1.8   #wind speed sensor height [m]
Z0<-0.1    #relative roughness [-]
## Wind speed at heigth of rain gauge
U<-AWS$ws*(log(Ztb)-log(Z0))/(log(Zws)-log(Z0))  #in R: log() = ln()
# ##check
# idx<-1:200
# plot(U[idx],type='l',col='darkorange')
# lines(AWS[idx,'ws'])

# ## Type of precipitation: solid or liquid 
prec_part$liq_rel<-0
prec_part[prec_part$TB > 0,'liq_rel']<-prec_part[prec_part$TB > 0,'liq'] / 
  prec_part[prec_part$TB > 0,'TB']
prec_part$sno_rel<-0
prec_part[prec_part$TB > 0,'sno_rel']<-prec_part[prec_part$TB > 0,'sno'] / 
  prec_part[prec_part$TB > 0,'TB']

## Catch ratio (partition between liquid and snow)
m<-vector()
for(i in 1:nrow(prec_part)){
  ifelse(prec_part[i,'liq_rel'] >= 0.5,
         m[i]<-cc_r,
         m[i]<-cc_s
  )
}
CR<-1/(1+m*U)

Pcor<-AWS$Rain_TB/CR #adjusted precipitation [mm/hr]

df<-cbind.data.frame(AWS$Date,AWS$Rain_TB,Pcor,CR)
colnames(df)<-c('Date','TB','corr','zCR')
df<-aggregate(df[,-1],list(Date=cut(as.POSIXct(df$Date),'day')),
              sum,na.rm=TRUE)
df[,-c(1,4)]<-cumsum(df[,-c(1,4)])
df$zCR<-df$zCR/24*1000