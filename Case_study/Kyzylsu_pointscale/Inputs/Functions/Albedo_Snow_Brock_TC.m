%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Albedo_Snow_Properties     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Asno,AMTa_out]=Albedo_Snow_Brock_TC(SWEtm1,AMTatm1,Psn_Th,Psn_sum,Ta_day,Br_Param,Cdeb,SWE,Aice,Deb_Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT
%%% MTa = maximum air temperature for the day divided so hourly
%%% MTatm1 = maximum air temperature for the day divided so hourly at time
%%% before
%%% Psn_max = maximum snow precipitation over previous 24 hours
%%% Psn_th = precipitation threshold (this is applied to an hourly value
%%% but found over the previous 24 hours)
%%% OUTPUT
% tau_sno [] %% Relative Age of snow
%snow_alb.dir_vis
%snow_alb.dif_vis
%snow_alb.dir_nir
%snow_alb.dif_nir
%%AMTa_out = accumulated maximum air temperature as output in tau_sno 
%%%%%%%%%%%%%%

%   Snow albedo function as described in Brock et al. (2000) and
%   implemented as in T and C

%   Snow albedo function as described in Brock et al. (2000)
%%%%
MTa=max(Ta_day)/24;
MTa(MTa<0)=0; %Only accumulate positive maximum temperature
%%%%%%%%

%%%%%%%%%%
SWE_m=SWE/1000; %Convert SWE to snow depth m w.e.

%Calculate maximum accumulated Ta
if (Psn_sum >= Psn_Th) || SWEtm1==0 %Snow greater than threshold or there was no snow on the timestep before
    ATa=MTa; %Reset to only be accumulated max T for this hour as new snow
    AMTa_out=ATa;
else %%(Psn_max < Psn_Th) Pr_sno<threshold
    ATa=MTa + AMTatm1; %Accumulated max T this timestep
    AMTa_out=ATa; %to go back as output for next timestep
end

%Determine underlying surface albedo
if Cdeb==0
    a_u=Aice; %clean ice
else
    a_u=Deb_Par.alb; %debris-covered ice
end

%Deep snow albedo
a_ds = Br_Param.a - Br_Param.b*log10(ATa);

%Shallow snow albedo
a_ss = a_u + Br_Param.c*exp(Br_Param.d*ATa);

%Final snow albedo
Asno = (1-exp(-SWE_m/Br_Param.e))*a_ds + exp(-SWE_m/Br_Param.e)*a_ss;

if Asno>0.85
    Asno=0.85;
elseif isnan(Asno)
    Asno=a_u; %So if=NaN (only occurs under conditions of Csno=1, SWE=0 and ATa=0, so its cold but the snow has melted, therefore give the underlying albedo)
end
   
return