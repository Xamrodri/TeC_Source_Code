
%==================================================
%
%                       Precipitation Separation Scheme
%
%       This scheme is to discrimintate rain, sleet, and snow by using daily
%       average data. 
%
%       Inputs are daily precipitation, daily mean surface air temperature,
%       relative humidity, air pressure, and station elevation.
%
%       There are two output files.
%       (1) precipitation type. There are 4 types of precipitation:
%               type = 0 : no precipitation;
%               type = 1 : rain;
%               type = 2 : sleet (as a mixture of rain and snow);
%               type = 3 : snow.
%       (2) solid fraction of precipitation. It ranges from 0 to 1. When
%       there is no precipitation, the value is recorded as -1.
%
%       Reference:
%           Ding, B., Yang, K., Qin, J., Wang, L., Chen, Y., He, X., 2014. The 
%           dependence of precipitation types on surface elevation and meteorological 
%           conditions and its parameterization. J. Hydrol. 513, 154¨C163. 
%           http://dx.doi.org/10.1016/j.jhydrol.2014.03.038.
%
%       Institute of Tibetan Plateau Research, Chinese Academy of Sciences
%       E-mail: dingbh@itpcas.ac.cn
%       Web: http://dam.itpcas.ac.cn/
%
%==================================================
%
%       Author: Baohong DING
%       08/04/2014
%
%==================================================

function [type,solid_fraction,Snow_acc,Tw] = Precip_Type(precip,Ta,rh,pres,ele)

%% Unit conversion

Ta = Ta- 273.15;                  % Convert Temperature from Kelvin to Celsius for this function
rh = rh / 100;                     % relative humidity, [%] => [none]
ele = ele / 1000;                  % elevation, [m] => [km]
pres = pres / 100;

%% Wet bulb temperature

Tw = Wetbulb_Temp( Ta, rh, pres );            % wet bulb temperature [degree C]

%% Thresholds for discrimination precipitation

dT = 0.215 - 0.099.*rh + 1.018.*rh .*rh;
dS = 2.374 - 1.634.*rh;  
T0 = -5.87 - 0.1042.*ele + 0.0885.*ele .*ele + 16.06.*rh - 9.614.*rh .*rh;

T_thre1 = T0;
T_thre2 = T0;
    
num =  dT ./ dS > 0.6931 ;
T_thre1(num) = T0(num) - dS(num) .* log( exp( dT(num)./dS(num) ) - 2 * exp ( - dT(num)./dS(num) ) );
T_thre2(num) = 2 * T0(num) - T_thre1(num);

%% Discrimination of rain, sleet, and snow

num_rain = ( Tw - T_thre2 >= 0 );
num_sleet = ( Tw - T_thre1 > 0 & Tw - T_thre2 < 0 );
num_snow = ( Tw - T_thre1 <= 0 );

type(num_rain == 1) = 1;                                 % rain
type(num_sleet == 1) = 2;                               % sleet
type(num_snow == 1) = 3;                              % snow

num_NoPrecip = ( precip == 0 );
type(num_NoPrecip == 1) = 0;                         % no precipitation

%% Calculation of solid fraction of precipitation

solid_fraction = 1 ./ ( 1 + exp( ( Tw - T0 ) ./ dS) );      % sleet
solid_fraction( type==1 ) = 0;                               % rain
solid_fraction( type==3 ) = 1;                               % snow
solid_fraction( type==0 ) = -1;                              % no precipitation

%% Output of precipitation type & solid fraction of precipitation

Snow_acc(type <= 2) = 0;
Snow_acc(type == 3) = precip(type == 3);
