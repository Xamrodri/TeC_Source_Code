%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% AUTOMATIC COMPUTATION OF RADIATION PARTIOTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Automatic_Radiation_Partition_I(Date,Lat,Lon,Zbas,DeltaGMT,Pr,Tdew,Rsw,GRAPH_VAR)
%%% PROGRAM  INPUTS
%%%% load Res_Location_Name.mat
%%%D N Pr Tdew DeltaGMT Lon Lat Zbas
%%% GRAPH_VAR 0/1 switcher to plot outputs 
%%% Rsw --> Measured Solar radiation
%load([cur_dir,'\Res_',Location_Name,'.mat']);
D=Date; N=zeros(size(D,1),1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
D=reshape(D,length(D),1);
[Anno,Mese,Giorno,Ora,Min,Sec]=datevec(D);
Ora= Ora + (Min/30>1);
Datam=[Anno,Mese,Giorno,Ora]; %% Datam %% [Yr, MO, DA, HR]
%%%%%%%%%%%%%%%%%%%%%%%%%%
dt= 1;
%%% Parameters LWP0 beta omega_A1 omega_A2
%%%%  Clear Sky beam
%beta = 0.0;%%0.045;%%0.095;%%% beta   Angstrom turbidity parameters [0-0.5]
%beta_A= 0.04;
uo  = 0.35; %% 0.28;  %%% [cm] ozone amount in vertical column [0.22-0.3]
un=  0.0002 ; % [cm] total nitrogen dioxide amount
alpha_A= 1.3;%% Angstrom turbidity parameters
%%% Clear Sky diffuse
omega_A1     = 0.92;%%0.920;%%%omega_A1  %% aerosol single-scattering albedo band 1 [0.84-0.94]
omega_A2     = 0.84;%%0.833;%%omega_A2 %% aerosol single-scattering albedo band 2  [0.84-0.94]
rho_g = 0.15; %% [] spatial average regional albedo
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Pc.beta_A=0.05; %% Worldwide average with alpha_A =1.3;
%if abs(Lat)<15;
if (Lat)<15 & Lat>-15
    Pc.beta_A=[0.054    0.054   0.054    0.054   0.054   0.054    0.054    0.054    0.054 0.054   0.054    0.054]; %%
    Pc.LWP_R=[70 70  70  70  70   70   70  70   70   70   70   70]; %%
elseif Lat<-15
    Pc.beta_A=[ 0.0850    0.0638    0.0454    0.0340    0.0352 0.0375    0.0500    0.0596    0.0789    0.0908    0.1049    0.1002  ]; %%
    Pc.LWP_R=[53.2845   57.1041   63.6880   75.5263   92.6457  109.3802 91.6259   68.6851   65.4887   60.8993   53.8308   52.2543 ]; %%
else
    Pc.beta_A=[  0.0375    0.0500    0.0596    0.0789    0.0908    0.1049    0.1002    0.0850    0.0638    0.0454    0.0340    0.0352]; %%
    Pc.LWP_R=[91.6259   68.6851   65.4887   60.8993   53.8308   52.2543   53.2845   57.1041   63.6880   75.5263   92.6457  109.3802]; %%
end
%%%%
NT=length(Datam(:,1));
NI = min(8760,NT); Ina = isnan(Rsw);
Rsw_I = Rsw(1:NI); 
if sum(Ina(1:NI))>= 2190
    Rsw_I = Rsw(NI+1:NI*2);
    if sum(Ina(NI+1:NI*2))>= 2190
        Rsw_I = Rsw(NI*2+1:NI*3);
        if sum(Ina(NI*2+1:NI*3))>= 2190
            Rsw_I = Rsw(NI*3+1:NI*4);
        end
    end
end
T_BF=[1.5 1.25 1 0.75 0.5 0.25 0 -0.25 -0.5 -0.67];
T_AF=[-0.5 -0.25 0 0.25 0.5 0.75 1  1.25  1.5  1.67];
for IS=1:length(T_BF)
    t_bef=T_BF(IS);
    t_aft=T_AF(IS); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i=0;
    SB=zeros(1,NI); SD=zeros(1,NI);SAB1=zeros(1,NI);SAB2=zeros(1,NI);SAD1=zeros(1,NI);SAD2=zeros(1,NI);PARB=zeros(1,NI);PARD=zeros(1,NI);
    for i=1:NI
        LWP0=Pc.LWP_R(Mese(i));
        beta_A=Pc.beta_A(Mese(i));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [SB(i),SD(i),SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i)]=ComputeRadiationForcingsN(Datam(i,:),DeltaGMT,Lon,Lat,...
            LWP0,N(i),Zbas,Tdew(i),beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g,t_bef,t_aft);
    end
    Rsws=SB+SD;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%
    n=length(Rsw_I); fr=24/dt;
    %%%%%%%%%%%%%%%%%%%%
    m=floor(n/fr); 
    Rswp=reshape(Rsw_I(1:m*fr),fr,m);
    Rswsp=reshape(Rsws(1:m*fr),fr,m);
    Xm=nanmean(Rswp');
    Xs=mean(Rswsp');
    rr = mean(Xs)./mean(Xm); Xs=Xs./rr; 
    vdiff=Xs-Xm; vdiff(vdiff<0)=2*vdiff(vdiff<0); 
    DIFF_VEC(IS) = sum((abs(vdiff)));
    if GRAPH_VAR==1
        t=1:fr;
        figure(600)
        subplot(5,2,IS);
        set(gca,'FontSize',10);
        plot(t,Xm,'r','LineWidth', 2);
        hold on ; grid on;
        plot(t,Xs,'b','LineWidth', 2);
        title('Radiation Total Daily cycle')
        legend('OBS','SIM')
        xlabel('Hour'); ylabel('[W/m^2]')
        mR=max(max(Xm,Xs))+50;
        axis([0 24 0 mR])
        %%%%%%%%%%%%
    end
end
%%% between 3 and 7 the most probable ;;
DIFF_VEC(3:7)=DIFF_VEC(3:7)-120; 
[v,p]=min(DIFF_VEC);
IS=p; clear v p
t_bef=T_BF(IS);
t_aft=T_AF(IS);
if GRAPH_VAR==1
    figure(620)
    set(gca,'FontSize',10);
    plot(T_AF,DIFF_VEC,'r','LineWidth', 2);
    xlabel('T_aft'); ylabel('Optimization'); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=0;
SB=zeros(1,NT); SD=zeros(1,NT);SAB1=zeros(1,NT);SAB2=zeros(1,NT);SAD1=zeros(1,NT);SAD2=zeros(1,NT);PARB=zeros(1,NT);PARD=zeros(1,NT);
if GRAPH_VAR==1; bau = waitbar(0,'Waiting Please...'); end 
for i=1:NT
    if GRAPH_VAR==1; waitbar(i/NT,bau); end 
    LWP0=Pc.LWP_R(Mese(i));
    beta_A=Pc.beta_A(Mese(i));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [SB(i),SD(i),SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i)]=ComputeRadiationForcingsN(Datam(i,:),DeltaGMT,Lon,Lat,...
        LWP0,N(i),Zbas,Tdew(i),beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g,t_bef,t_aft);
end
if GRAPH_VAR==1; close(bau) ; end 
%%%%%%%%
Rsws=SB+SD;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%save('nnn.mat','Location_Name','Rsw','Rsws','Pr','D','t_bef','t_aft','LWP0',...
%    'beta_A','alpha_A','omega_A1','omega_A2','uo','un','rho_g','Pc','dt','GRAPH_VAR','cur_dir'); 
%%%
%clear 
%%%
%load('nnn.mat'); 
[Nsim]=Estimate_CloudCover(Rsw,Rsws,Pr,D); 
%%%% 
%load([cur_dir,'\Res_',Location_Name,'.mat']);
D=Date; N=Nsim;
N(isnan(Rsw))=0.15;
%N(isnan(N))=nanmean(N);
%delete('nnn.mat'); 
%%%%%%%%%%%%%%%%%%%
%%%
i=0;
SB=zeros(1,NT); SD=zeros(1,NT);SAB1=zeros(1,NT);SAB2=zeros(1,NT);SAD1=zeros(1,NT);SAD2=zeros(1,NT);PARB=zeros(1,NT);PARD=zeros(1,NT);
if GRAPH_VAR==1; bau = waitbar(0,'Waiting Please...'); end 
for i=1:NT
   if GRAPH_VAR==1; waitbar(i/NT,bau); end 
    LWP0=Pc.LWP_R(Mese(i));
    beta_A=Pc.beta_A(Mese(i));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [SB(i),SD(i),SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i)]=ComputeRadiationForcingsN(Datam(i,:),DeltaGMT,Lon,Lat,...
        LWP0,N(i),Zbas,Tdew(i),beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g,t_bef,t_aft);
end
if GRAPH_VAR==1; close(bau); end 
Rsws=SB+SD;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsw(isnan(Rsw))=Rsws(isnan(Rsw));
if GRAPH_VAR==1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('Rsw SIM:')
    mean(Rsws)
    disp('RSW OBS:')
    mean(Rsw)
%%%%%% 
    n=length(Rsw); fr=24/dt;
    %%%%%%%%%%%%%%%%%%%%
    m=floor(n/fr);
    Rswp=reshape(Rsw(1:m*fr),fr,m);
    Rswsp=reshape(Rsws(1:m*fr),fr,m);
    Xm=nanmean(Rswp');
    Xs=mean(Rswsp');
    t=1:fr;
    figure(100)
    subplot(1,1,1);
    set(gca,'FontSize',10);
    plot(t,Xm,'r','LineWidth', 2);
    hold on ; grid on;
    plot(t,Xs,'b','LineWidth', 2);
    title('Radiation Total Daily cycle')
    legend('OBS','SIM')
    xlabel('Hour'); ylabel('[W/m^2]')
    mR=max(max(Xm,Xs))+50;
    axis([0 24 0 mR])
    %%%%%%%%%%%%  
end
[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD]=Ratio_Evaluator(Rsws,Rsw,SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD); 
end
function[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD]=Ratio_Evaluator(Rsws,Rsw,SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsw=reshape(Rsw,1,length(Rsw));  Rsw(isnan(Rsw))=Rsws(isnan(Rsw));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rd1=SAD1./(SD+SB); rd1(isnan(rd1))=0; rd1(rd1>1)=1;
rd2=SAD2./(SD+SB); rd2(isnan(rd2))=0; rd2(rd2>1)=1;
rb1=SAB1./(SD+SB); rb1(isnan(rb1))=0; rb1(rb1>1)=1;
rb2=SAB2./(SD+SB); rb2(isnan(rb2))=0; rb2(rb2>1)=1;
rparb=PARB./(SD+SB); rparb(isnan(rparb))=0; rparb(rparb>1)=1;
rpard=PARD./(SD+SB); rpard(isnan(rpard))=0; rpard(rpard>1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SAD1 = rd1.*Rsw; 
SAB1 = rb1.*Rsw; 
SAD2 = rd2.*Rsw; 
SAB2 = rb2.*Rsw; 
PARB = rparb.*Rsw; 
PARD = rpard.*Rsw; 
%%%%%%%%%%%%%%%%%%%%%
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Estimate Cloud_Cover       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Nsim]=Estimate_CloudCover(Rsw,Rsws,Pr,D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rsw=reshape(Rsw,length(Rsw),1); 
Rsws=reshape(Rsws,length(Rsws),1); 
Pr=reshape(Pr,length(Pr),1); 
D=reshape(D,length(D),1); 
%%%%%%%%%%%%%%%%%%%
Nsim = NaN*ones(length(Rsw),1);
A = NaN*ones(length(Rsw),1);
%%%%%%%%%%%%%%%%%%
RswN=Rsw(Rsws>150);
Rsw0=Rsws(Rsws>150);
%%%%%%%%%%%%%%%%%%%%%%%
A(Rsws>150)=RswN./Rsw0;
A(A>1)=1; A(A<0)=0;
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a= 0.75; b= 3.4;
%%% --->  A= 1-aN^b
%%%%%%%%%%%%%%%%%%%
Nsim = ((1/a)*(1-(A))).^(1/b);
Nsim(A==1)=0;
Nsim(Nsim>1)=1;
%%%%%%%%%%%%%%%%%%%%%%%%
if isnan(Nsim(1))
    Nsim(1) = nanmean(Nsim);
end
if isnan(Nsim(end))
    Nsim(end) = nanmean(Nsim);
end
x=1:length(Nsim);
xnt= x(not(isnan(Nsim)));
xn=x(isnan(Nsim));
Nn=interp1(xnt,Nsim(not(isnan(Nsim))),xn,'linear');
Nsim(xn)=Nn;
Nsim(Pr>0)=1;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute radiation forcings %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[EB,ED,EB1,EB2,ED1,ED2,PARB,PARD]=ComputeRadiationForcingsN(Datam,DeltaGMT,Lon,Lat,...
    LWP0,N,Zbas,Tdew,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g,t_bef,t_aft)
%%%OTUPUT
%%%  EB %% sky beam irradiance without terrain effects
%%%% ED %% sky total diffuse irradiance  without terrain effects
%%EB1 [W/m^2] sky beam irradiance VIS band[0.29 um - 0.70um ]
%%EB2  [W/m^2]   sky beam irradiance NIR band [0.70 um - 4.0 um ]
%%ED1 [W/m^2]    sky total diffuse flux at the ground  VIS band [0.29 um -0.70um ]
%%%ED2   [W/m^2]  sky  total diffuse flux at the ground NIR band  [0.70 um-4.0 um ]
%%%%%%%%%%
%%% INPUT
%%% Datam %% [Yr, MO, DA, HR]
%%% DeltaGMT [°]
%%% Lon [°]
%%% Lat [°]
%%% LWP0 [g/m^2]
%%% N [0-1] cloudiness
%Zbas [m a.s.l.]  watershed elevation
%Tdew [°C] dew point temperature
%beta_A   0.05 [0-1] Angstrom turbidity parameters
%alpha_A  1.3 [0.5-2.5]Angstrom turbidity parameters
%omega_A1  %% 0.92 [0.74-0.94]aerosol single-scattering albedo band 1
%omega_A2 %% 0.84 [0.74-0.94] aerosol single-scattering albedo band 2
%%% uo 0.35 [0.22-0.35]  %% [cm] ozone amount in vertical column
%%% un 0.0002 [0.0001-0.046] %% [cm] total nitrogen dioxide amount
%%%  rho_g 0.15 [0.05-0.3] %% [] spatial average regional albedo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants
So    = 1366.1;  % [W/m^2] Solar constant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[delta_S,h_S,zeta_S,T_sunrise,T_sunset,L_day,E0,jDay,Delta_TSL] = SetSunVariables(Datam,DeltaGMT,Lon,Lat,t_bef,t_aft);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction SetSunVariables   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sop = So*E0; %% Actual Solar constant [W/m^2]
%%%
%%% Partition Energy two bands Gueymard (2004; 2008)
So1 = Sop*0.4651;  %%  [W/m^2]Extraterrestrial  Radiation VIS band [0.29 um - 0.70 um ]
So2 = Sop*0.5195; %% [W/m^2] %%Extraterrestrial Radiation NIR band [0.70 um - 4.0 um ]
%LWP = exp(N*log(LWP0)) - 1; %%% Liquid Water Path whit cloudiness  [g/m^2]
LWP = LWP0*N; %% Liquid Water Path whit cloudiness  [g/m^2]
if(LWP < 1.1)
    Ntmp = 0; %% Cloudiness
else
    Ntmp = N; % Cloudiness
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(Ntmp,0) % Completely clear sky
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Subfunction SetClearSkyRadiation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [EB1,EB2,ED1,ED2,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] = SetClearSkyRadiation(h_S,Zbas,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
    %%%%%%%%%%%%%%%%%%%%%%
elseif Ntmp > 0  %% Overcast Condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %     Subfunction SetClearSkyRadiation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Eb1,Eb2,Ed1,Ed2,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] = SetClearSkyRadiation(h_S,Zbas,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Subfunction SetCloudySkyRadiation      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [EB1,EB2,ED1,ED2] = SetCloudySkyRadiation(LWP,h_S,Eb1,Eb2,Edp1,Edp2,N,rho_s1,rho_s2,rho_g); %
    %%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%Correction Flat Surface
EB1=EB1*sin(h_S);
EB2=EB2*sin(h_S);
EB = EB1 + EB2;%% [W/m^2] % beam irradiance
ED = ED1 + ED2; %% [W/m^2] %total diffuse flux at the ground
%%%% PAR Radiation estimation
PARB = EB1*Mb;
PARD = Mg*(EB1+ED1)- PARB;
%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction SetSunVariables   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta_S,h_S,zeta_S,T_sunrise,T_sunset,L_day,E0,jDay,Delta_TSL] = SetSunVariables(Datam,DeltaGMT,Lon,Lat,t_bef,t_aft)
%%% INPUT
%%% Datam %% [Yr, MO, DA, HR]
%%% DeltaGMT [°]
%%% Lon [°]
%%% Lat [°]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT
%%delta_S,  Solar declination
%%tau_S,  Hour angle of the sun
%%h_S, [rad] solar altitude
%%zeta_S, Sun's azimuth
%%T_sunrise, [h]  sunrise time,
%%T_sunset,  [h]sunset time,
%%L_day, [h] total day length
%%r_ES,[] ratio of the actual Earth-Sun to the mean Earth-Sun  distance
%%jDay,    Julian Day
%%Delta_TSL [h] Time difference between standard and local meridian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the julian day of the current time
days = [31 28 31 30 31 30 31 31 30 31 30 31];
nowYR=Datam(1); nowMO=Datam(2); nowDA=Datam(3); nowHR =Datam(4);
if(nowMO==1)
    jDay = nowDA;
elseif(nowMO==2)
    jDay = days(1) + nowDA;
else
    jDay = sum(days(1:(nowMO-1))) + nowDA;
    if(mod(nowYR,4)==0)
        if(mod(nowYR,400)==0)
            jDay = jDay + 1;
        elseif(mod(nowYR,100)~=0)
            jDay = jDay + 1;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute solar declination
delta_S = 23.45*pi/180*cos(2*pi/365*(172 - jDay));
%%% Compute time difference between standard and local meridian
if(Lon<0)
    Delta_TSL = -1/15*(15*abs(DeltaGMT) - abs(Lon));
else
    Delta_TSL = 1/15*(15*abs(DeltaGMT) - abs(Lon));
end
%t_bef=0.5;%
%t_aft=0.5;%
t= nowHR-t_bef:0.0166666:nowHR+t_aft;
for i=1:length(t)
    %%%  Compute hour angle of the sun
    if(t(i) < (12 + Delta_TSL))
        tau_S(i) = 15*pi/180*(t(i) + 12 - Delta_TSL);
    else
        tau_S(i) = 15*pi/180*(t(i) - 12 - Delta_TSL);
    end
end
%%%% Compute solar altitude
Lat_rad = Lat*pi/180;
sinh_S = sin(Lat_rad)*sin(delta_S) + cos(Lat_rad)*cos(delta_S)*cos(tau_S);
h_S = asin(sinh_S);
h_S=mean(h_S);
%%%% Compute Sun's azimuth
zeta_S = atan(-sin(tau_S)./(tan(delta_S)*cos(Lat_rad) - sin(Lat_rad)*cos(tau_S)));
%%%%%%%%%%%%%%%%%%%
for i=1:length(t)
    if (tau_S(i) >0 && tau_S(i) <= pi)
        if (zeta_S(i) > 0.)
            zeta_S(i) =  zeta_S(i) + pi;
        else
            zeta_S(i) = zeta_S(i) + (2.*pi);
        end
    elseif (tau_S(i) >=pi && tau_S(i) <= 2*pi)
        if (zeta_S(i) < 0.)
            zeta_S(i) =  zeta_S(i) + pi;
        end
    end
end
zeta_S=mean(zeta_S);
%%%% Compute sunrise time, sunset time, and total day length
T_sunrise = 180/(15*pi)*(2*pi - acos(-tan(delta_S)*tan(Lat_rad))) - 12;
T_sunset  = 180/(15*pi)*acos(-tan(delta_S)*tan(Lat_rad)) + 12;
L_day     = 360/(15*pi)*acos(-tan(delta_S)*tan(Lat_rad));
% Compute the ratio of the actual Earth-Sun to the mean Earth-Sun  distance
GA = 2*pi*(jDay-1)/365; %% daily angle
E0=1.00011+0.034221*cos(GA)+0.00128*sin(GA)+0.000719*cos(2*GA)+ ...
    0.000077*sin(2*GA); %%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Subfunction SetClearSkyRadiation      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Eb1,Eb2,Ed1,Ed2,Edp1,Edp2,rho_s1,rho_s2,Mb,Mg] = SetClearSkyRadiation(h_S,Zbas,Tdew,So1,So2,beta_A,alpha_A,omega_A1,omega_A2,uo,un,rho_g)
%%%% INPUT
%%h_S,[rad]  solar altitude
%Zbas [m a.s.l.]  watershed elevation
%Tdew [°C] dew point temperature
%So1 [W/m^2]  Extraterrestrial  Radiation VIS band [0.29 um - 0.70 um ]
%So2  % [W/m^2] %%Extraterrestrial Radiation NIR band [0.70 um - 4.0 um ]
%beta_A   0.05 [0-1] Angstrom turbidity parameters
%alpha_A  1.3 [0.5-2.5]Angstrom turbidity parameters
%omega_A1  %% 0.92 [0.74-0.94]aerosol single-scattering albedo band 1
%omega_A2 %% 0.84 [0.74-0.94] aerosol single-scattering albedo band 2
%%% uo 0.35 [0.22-0.35]  %% [cm] ozone amount in vertical column
%%% un 0.0002 [0.0001-0.046] %% [cm] total nitrogen dioxide amount
%%%  rho_g 0.15 [0.05-0.3] %% [] spatial average regional albedo
%%% OUTPUT
%%Eb1 [W/m^2]  clear sky beam irradiance VIS band [0.29 um - 0.70 um ]
%%Eb2  [W/m^2]  clear sky beam irradiance NIR band [0.70 um - 4.0 um ]
%%Ed1 [W/m^2]  clear sky total diffuse flux at the ground  VIS band [0.29 um - 0.70um ]
%%%Ed2   [W/m^2]  clear sky total diffuse flux at the ground NIR band [0.70 um - 4.0um ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_SD = 180/pi*h_S; %% [°]  solar altitude
Z=90-h_SD; %% sun zenit angle [°]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(h_S > 0.0) %% if there is sunshine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIRECT-BEAM TRANSMITTANCES %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Generally p/po = exp(-g*d_elevat/(Rd*Tm));  Rd =287.05; %% [J/kgK]  g= 9.81; %% [m/s^2] Tm = 15°C
    p=1013.25*exp(-Zbas/8434.5); %%% [mbar] Pressure after correction for differences in pressure between basin and seal level
    w=  exp(0.07*Tdew - 0.075); %% precipitable water [cm]
    %%% Gueymard 2003
    mR=(sin(h_S)+(0.48353*Z^(0.095846))/(96.741-Z)^1.1754)^-1; %%% Rayleigh scattering and uniformly mixd gas  Air mass
    mO=(sin(h_S) +(1.0651*Z^0.6379)/((101.8-Z)^2.2694))^-1; %% ozone absorption Air mass
    %% not necessary %% mn=(sin(h_S) +(1.1212*Z^1.6132)/((111.55-Z)^3.2629))^-1; %%
    mW=(sin(h_S) +(0.10648*Z^0.11423)/((93.781-Z)^1.9203))^-1;%% water vapor Air mass
    mA=(sin(h_S) +(0.16851*Z^0.18198)/((95.318-Z)^1.9542))^-1; %% aerosol extinction Air mass
    mRp=(p/1013.25)*mR;%% Rayleigh scattering and uniformely mixd gas  Air mass Correction pressure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute ozone transmittances
    f1= uo*(10.979-8.5421*uo)/(1+2.0115*uo+40.189*uo^2);
    f2= uo*(-0.027589-0.005138*uo)/(1-2.4857*uo+13.942*uo^2);
    f3= uo*(10.995-5.5001*uo)/(1+1.6784*uo + 42.406*uo^2);
    TO1 = (1+f1*mO+f2*mO^2)/(1+f3*mO);
    TO2 = 1.0;
    % Compute nitrogen dioxide transmittances
    g1 = (0.17499 +41.654*un -2146.4*un^2)/(1+22295.0*un^2);
    g2= un*(-1.2134+ 59.324*un)/(1+ 8847.8*un^2);
    g3= (0.17499 +61.658*un + 9196.4*un^2)/(1+74109.0*un^2);
    TN1 = min(1,(1+g1*mW+g2*mW^2)/(1+g3*mW));
    TN2 = 1.0;
    % Compute Rayleigh scattering transmittances
    TR1    = (1+ 1.8169*mRp -0.033454*mRp^2)/(1+ 2.063*mRp +0.31978*mRp^2);
    TR2    = (1 -0.010394*mRp)/(1- 0.00011042*mRp^2);
    % Compute the uniformly mixd gas transmittances
    TG1    = (1+ 0.95885*mRp -0.012871*mRp^2)/(1+ 0.96321*mRp +0.015455*mRp^2);
    TG2    = (1+ 0.27284*mRp -0.00063699*mRp^2)/(1+0.30306*mRp);
    % Compute water vapor transmmitances
    h1=w*(0.065445+0.00029901*w)/(1+1.2728*w);
    h2=w*(0.065687+0.0013218*w)/(1+1.2008*w);
    c1=w*(19.566-1.6506*w+1.0672*w^2)/(1+ 5.4248*w+1.6005*w^2);
    c2=w*(0.50158-0.14732*w+0.047584*w^2)/(1+ 1.1811*w+1.0699*w^2);
    c3=w*(21.286-0.39232*w+1.2692*w^2)/(1+ 4.8318*w+1.412*w^2);
    c4=w*(0.70992-0.23155*w+0.096514*w^2)/(1+ 0.44907*w+0.75425*w^2);
    TW1    = (1+h1*mW)/(1+h2*mW);
    TW2    = (1+c1*mW+c2*mW^2)/(1+c3*mW+c4*mW^2);
    % Compute average aerosol transmittances
    alpha1 = alpha_A; %%Angstrom turbidity parameters exponent
    alpha2 = alpha_A; %%Angstrom turbidity parameters exponent
    beta1 = beta_A*0.7^(alpha1-alpha2); %Angstrom turbidity parameters
    beta2 = beta_A ; %Angstrom turbidity parameters
    uA     = log(1 + mA*beta2);%%% ---> Gueymard, (1989) Parameterization for effective wavelength computation
    d0    = 0.57664 - 0.024743*alpha1 ;
    d1    = (0.093942 -0.2269*alpha1 + 0.12848*alpha1^2)/(1+0.6418*alpha1);
    d2 = (-0.093819 + 0.36668*alpha1 - 0.12775*alpha1^2)/(1-0.11651*alpha1);
    d3 = alpha1*(0.15232-0.087214*alpha1+0.012664*alpha1^2)/(1-0.90454*alpha1+0.26167*alpha1^2);
    e0 =(1.183 -0.022989*alpha2 + 0.020829*alpha2^2)/(1+0.11133*alpha2);
    e1 =(-0.50003 -0.18329*alpha2 + 0.23835*alpha2^2)/(1+1.6756*alpha2);
    e2 =(-0.50001 +1.1414*alpha2 + 0.0083589*alpha2^2)/(1+11.168*alpha2);
    e3 =(-0.70003 -0.73587*alpha2 + 0.51509*alpha2^2)/(1+4.7665*alpha2);
    le1    = (d0 + d1*uA + d2*uA^2)/(1+d3*uA^2);  %% effective wavelength for band 1
    le2    = (e0 + e1*uA + e2*uA^2)/(1+e3*uA);  %% effective wavelength for band 2  %%% <--
    tauA1= beta1*le1^(-alpha1);
    tauA2= beta2*le2^(-alpha2);
    TA1    = exp(-mA*tauA1);
    TA2    = exp(-mA*tauA2);
    % Return the transmittances for each band
    T1(1)  = TO1;
    T1(2)  = TR1;
    T1(3)  = TG1;
    T1(4)  = TW1;
    T1(5)  = TA1;
    T1(6)  = TN1;
    %%%
    T2(1)  = TO2;
    T2(2)  = TR2;
    T2(3)  = TG2;
    T2(4)  = TW2;
    T2(5)  = TA2;
    T2(6)  = TN2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute clear sky beam irradiance
    Eb1 = So1*prod(T1);
    Eb2 = So2*prod(T2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   DIFFUSE TRANSMITTANCES   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mWp= 1.66; % water vapor optical mass of reference
    TW1p    = (1+h1*mWp)/(1+h2*mWp);
    TW2p    = (1+c1*mWp+c2*mWp^2)/(1+c3*mWp+c4*mWp^2);
    TN1p = min(1,(1+g1*mWp+g2*mWp^2)/(1+g3*mWp));
    TN2p = 1.0;
    clear g1 g2 g3 h1 h2
    %%%%%%%%%%
    TAs1 = exp(-mA*omega_A1*tauA1);  %% aerosol transmittances due to scattering
    TAs2 = exp(-mA*omega_A2*tauA2);  %% aerosol transmittances due to scattering
    %%%%%%%%%%%%%%%%%%%%%%%
    BR1   = 0.5*(0.89013-0.0049558*mR+0.000045721*mR^2);  %% forward scattering fractions for Rayleigh extinction
    BR2  =  0.5 ; %% forward scattering fractions for Rayleigh extinction
    Ba   = 1 - exp(-0.6931 - 1.8326*sin(h_S));%% fractions aerosol scattered fluxes
    %%%%%%%%%%
    g0= (3.715 +0.368*mA +0.036294*mA^2)/(1+0.0009391*mA^2);
    g1= (-0.164-0.72567*mA +0.20701*mA^2)/(1+0.0019012*mA^2);
    g2= (-0.052288+0.31902*mA +0.17871*mA^2)/(1+0.0069592*mA^2);
    h0= (3.4352 +0.65267*mA +0.00034328*mA^2)/(1+0.034388*mA^1.5);
    h1= (1.231 -1.63853*mA +0.20667*mA^2)/(1+0.1451*mA^1.5);
    h2= (0.8889 -0.55063*mA +0.50152*mA^2)/(1+0.14865*mA^1.5);
    F1=(g0+g1*tauA1)/(1+g2*tauA1);
    F2=(h0+h1*tauA2)/(1+h2*tauA2);
    %%%% Compute the clear  sky albedo %%%%
    rho_s1=(0.13363 +0.00077358*alpha1 ...
        +beta1*(0.37567+0.22946*alpha1)/(1-0.10832*alpha1))/(1+beta1*(0.84057+0.68683*alpha1)/(1-0.08158*alpha1));
    rho_s2=(0.010191 +0.00085547*alpha2 ...
        +beta2*(0.14618+0.062758*alpha2)/(1-0.19402*alpha2))/(1+beta2*(0.58101+0.17426*alpha2)/(1-0.17586*alpha2));
    % Compute incidente irradiance a perfectly absorbing ground
    Edp1= TO1*TG1*TN1p*TW1p*(BR1*(1-TR1)*(TA1^0.25) +Ba*F1*TR1*(1-TAs1^0.25))*So1*sin(h_S);
    Edp2= TO2*TG2*TN2p*TW2p*(BR2*(1-TR2)*(TA2^0.25) +Ba*F2*TR2*(1-TAs2^0.25))*So2*sin(h_S);
    %% Backscattered diffuse component
    Edd1= rho_g*rho_s1*(Eb1*sin(h_S) +Edp1)/(1-rho_g*rho_s1);
    Edd2= rho_g*rho_s2*(Eb2*sin(h_S) +Edp2)/(1-rho_g*rho_s2);
    %%%% %%% Compute diffuse irradiances
    Ed1= Edp1 +Edd1;
    Ed2= Edp2 +Edd2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Compute PAR Ratios
    m15 = min(mR,15);
    t0 = (0.90227 + 0.29*m15 + 0.22928*(m15^2)-0.0046842*(m15^3))/(1 +0.35474*m15 +0.19721*(m15^2));
    t1 = (-0.10591 + 0.15416*m15 - 0.048486*(m15^2)+0.0045932*(m15^3))/(1 -0.29044*m15 +0.026267*(m15^2));
    t2 = (0.47291 - 0.44639*m15 + 0.1414*(m15^2)-0.014978*(m15^3))/(1 -0.37798*m15 +0.052154*(m15^2));
    t3 = (0.077407+ 0.18897*m15 - 0.072869*(m15^2)+ 0.0068684*(m15^3))/(1 -0.25237*m15 +0.020566*(m15^2));
    v0 = (0.82725+ 0.86015*m15 + 0.00713*(m15^2)+ 0.00020289*(m15^3))/(1 +0.90358*m15 +0.015481*(m15^2));
    v1 = (-0.089088+ 0.089226*m15 - 0.021442*(m15^2)+ 0.0017054*(m15^3))/(1 -0.28573*m15 +0.024153*(m15^2));
    v2 = (-0.05342- 0.0034387*m15 + 0.0050661*(m15^2)- 0.00062569*(m15^3))/(1 -0.32663*m15 +0.029382*(m15^2));
    v3 = (-0.17797+ 0.13134*m15 - 0.030129*(m15^2)+ 0.0023343*(m15^3))/(1 -0.28211*m15 +0.023712*(m15^2));
    beta_e=beta1*(le1^(1.3-alpha1));
    Mb=(t0 +t1*beta_e +t2*(beta_e^2))/(1+t3*(beta_e^2)); %% PAR ratio for beam
    Mg=(v0 +v1*beta_e +v2*(beta_e^2))/(1+v3*(beta_e^2)); %% PAR ratio for global
else %%% not sunshine
    Eb1 = 0;
    Eb2 = 0;
    Ed1 = 0;
    Ed2 = 0;
    Edp1=0;
    Edp2=0;
    rho_s1=NaN;
    rho_s2=NaN;
    Mb=0;
    Mg=0;
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Subfunction SetCloudySkyRadiation      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EB1,EB2,ED1,ED2] = SetCloudySkyRadiation(LWP,h_S,Eb1,Eb2,Edp1,Edp2,N,rho_s1,rho_s2,rho_g)
%%%% INPUT
%%h_S,[rad]  solar altitude
%LWP % Liquid Water Path whit cloudiness  [g/m^2]
%%%Eb1 [W/m^2]  clear sky beam irradiance VIS band [0.29 um - 0.70 um ]
%%Eb2  [W/m^2]  clear sky beam irradiance NIR band [0.70 um - 4.0 um ]
%%Edp1 [W/m^2]  clear sky incidente irradiance  diffuse flux at the ground  VIS band [0.29 um - 0.70um ]
%%%Edp2   [W/m^2] clear sky incidente irradiance diffuse flux at the ground NIR band [0.70 um -4.0um ]
%%% N [0-1] cloudiness
%%% rho_s1 rho_s2   clear sky albedo
%%% rho_g ground albedo
%%% OUTPUT
%%EB1 [W/m^2]  cloud sky beam irradiance VIS band [0.29 um - 0.70 um ]
%%EB2  [W/m^2]  cloud sky beam irradiance NIR band [0.70 um - 4.0 um ]
%%ED1 [W/m^2]   cloud sky total diffuse flux at the ground  VIS band [0.29 um - 0.70um ]
%%%ED2   [W/m^2]  cloud sky  total diffuse flux at the ground NIR band [0.70 um - 4.0um ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(h_S > 0) % if there is sunshine
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Partition in Two principals band
    K = [0.4651 0.5195 0.5195 0.5195];
    %% 4 band approach of Slingo 1989
    %% Slingo (1989) considered four spectral  bands, one in UV/VIS and three in NIR wavelength
    %% intervals: [0.25um ÷ 0.69um], [0.69um ÷ 1.19um], [1.19um ÷ 2.38um], [2.38?m ÷ 4.0um]
    %% Parameterization of Slingo (1989)
    k = [0.460 0.326 0.181 0.033]; %% respective fraction of solar irradiance at the top of the atmosphere
    a = [2.817 2.682 2.264 1.281]*1e-2;  %%[m^2/g]
    b = [1.305 1.346 1.454 1.641]; %% [um m^2 /g]
    c = [-5.62e-8 -6.94e-6 4.64e-4 2.01e-1]; %%[]
    d = [1.63e-7 2.35e-5 1.24e-3 7.56e-3]; %%[1/um]
    e = [0.829 0.794 0.754 0.826];  %[]
    f = [2.482 4.226 6.560 4.353]*1e-3; %%[1/um]
    %%% Estimate effective radius of drop-size distribution
    tn1=(10^(0.2633 + 1.7095*log(log10(LWP)))); %% optical thickness band [0.29 um - 0.70 um ]
    tn2=(10^(0.3492 + 1.6518*log(log10(LWP)))); %% optical thickness band [0.70 um - 4.0um ]
    re1 = 1.5*LWP/tn1;  %% [um] effective radius of cloud-droplet size distribution  band [0.29 um - 0.70 um ]
    re2 = 1.5*LWP/tn2;  %% [um] effective radius of cloud-droplet size distribution [0.70 um - 4.0um ]
    re  = [re2 re1 re1 re1]; %%[um] effective radius of cloud-droplet size distribution 4 bands
    re(re < 4.2) = 4.2;  %% range re adjustament
    re(re > 16.6) = 16.6;  %% range re adjustament
    %%%% Compute cloudy sky direct beam irradiance and cloudy sky diffuse beam irradiance
    tau         = LWP.*(a + b./re);  %%  cloud optical depth 4 bands
    omega_tilde = 1 - (c + d.*re);  %% single  scatter  albedo
    g           = e + f.*re; %%  asymmetry parameter
    %%%%%%
    beta0   = (3/7).*(1 - g); %% fraction of the scattered diffuse radiation
    betah   = 0.5 - (3.*sin(h_S).*g)./(4.*(1 + g));  %% fraction of the scattered direct radiation
    f2      = g.^2;
    U1      = 7/4; %% reciprocals of the effective cosines for the diffuse upward flux
    U2      = (7/4).*(1 - (1 - omega_tilde)./(7.*omega_tilde.*beta0));  % reciprocals of the effective cosines for the diffuse  downward flux
    alpha1  = U1.*(1 - omega_tilde.*(1 - beta0));
    alpha2  = U2.*omega_tilde.*beta0;
    alpha3  = (1 - f2).*omega_tilde.*betah;
    alpha4  = (1 - f2).*omega_tilde.*(1 - betah);
    eps     = sqrt(alpha1.^2 - alpha2.^2);
    M       = alpha2./(alpha1 + eps);
    E       = exp(-eps.*tau);
    gamma1  = ((1 - omega_tilde.*f2).*alpha3 - sin(h_S).*(alpha1.*alpha3 + alpha2.*alpha4))./((1 - omega_tilde.*f2).^2 - eps.^2.*sin(h_S).^2);
    gamma2  = (-(1 - omega_tilde.*f2).*alpha4 - sin(h_S).*(alpha1.*alpha4 + alpha2.*alpha3))./((1 - omega_tilde.*f2).^2 - eps.^2.*sin(h_S).^2);
    %%%%%%%%%
    TDB         = exp(-(1 - omega_tilde.*f2).*tau./sin(h_S)); %% cloud transmissivity for the direct beam flux in 4 bands j
    TDB(TDB < 0.0) = 0.0;
    RDIF    = M.*(1 - E.^2)./(1 - E.^2.*M.^2);  %%% diffuse reflectivity for diffuse incident radiation
    RDIF(RDIF < 0.0) = 0.0;
    TDIF    = E.*(1 - M.^2)./(1 - E.^2.*M.^2);  %%%  diffuse transmissivity for diffuse incident radiation
    TDIF(TDIF < 0.0) = 0.0;
    RDIR    = (-gamma2.*RDIF - gamma1.*TDB.*TDIF + gamma1); %%% diffuse reflectivity for direct incident radiation
    TDIR    = (-gamma2.*TDIF - gamma1.*TDB.*RDIF + gamma2.*TDB);  %%% diffuse transmissivity for direct incident radiation
    %%%%%%%%%%%%
    for j=1:4  %%% check TDB + TDIR <= 1
        if  (TDB(j) + TDIR(j)) > 1
            SDB_DIR = TDB(j) + TDIR(j);
            TDIR(j) = TDIR(j)/SDB_DIR;
            TDB(j)  = TDB(j)/SDB_DIR;
        end
    end
    %%%%%%% cloud sky beam irradiance
    EB1     = Eb1*((1 - N) + TDB(1)*N)*k(1)/K(1);
    EB2_all = Eb2.*(TDB(2:4).*k(2:4))./K(2:4);
    EB2     = (1 - N)*Eb2 + N.*sum(EB2_all);
    %%%%%% cloud sky  total diffuse flux
    EDp1     = (1-N)*Edp1 + N*(TDIR(1).*Eb1 + TDIF(1).*Edp1).*k(1)./K(1);  %%
    EDp2_all = (TDIR(2:4).*Eb2 + TDIF(2:4).*Edp2).*k(2:4)./K(2:4);
    EDp2     = (1 - N).*Edp2 + N.*sum(EDp2_all);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Sky albedo for cloudy sky
    rho_csb1 = (1-N)*rho_s1 + N*RDIR(1)*k(1)/K(1);
    rho_csd1 = (1-N)*rho_s1 + N*RDIF(1)*k(1)/K(1);
    rho_csb2 = (1-N)*rho_s2 + N*sum(RDIR(2:4).*k(2:4)./K(2:4));
    rho_csd2 = (1-N)*rho_s2 + N*sum(RDIF(2:4).*k(2:4)./K(2:4));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Backscattered diffuse component  after sky attenuation
    EDd1= rho_g*rho_csb1*(EB1*sin(h_S))/(1-rho_g*rho_csb1) +  rho_g*rho_csd1*(EDp1)/(1-rho_g*rho_csd1);
    EDd2= rho_g*rho_csb2*(EB2*sin(h_S))/(1-rho_g*rho_csb2) +  rho_g*rho_csd2*(EDp2)/(1-rho_g*rho_csd2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Diffuse radiation
    ED1= EDp1 +EDd1;
    ED2= EDp2 +EDd2;
else % not sunshine
    EB1 = 0.0;
    EB2 = 0.0;
    ED1 = 0.0;
    ED2 = 0.0;
end
end
