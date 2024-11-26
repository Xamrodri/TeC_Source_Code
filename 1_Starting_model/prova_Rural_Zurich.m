%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% WORKING LAUNCH PAD HBM  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% INPUT MANAGER %%%%%%%%%%%%%%
current_directory = cd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NN=6096;% %%% time Step
%%%%%%%%% Time step 
dt=3600; %%[s] %%% 
dth=1; %%[h]
%%%%%%%%%%%%
ms = 16;  %% Number of soil layers 
cc = 1;   %% Crown area
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% METEO INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id_location = 'ZURICH_SMA';
load('M:\19_ISTA\1_TC\3_Model_Source\TeC_Source_Code\Inputs\Data_Run_Zurich_Fluntern.mat')
Date=D; clear D 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x1=1;
x2=6096; %Depends on the time step, defined by the data added.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Date=Date(x1:x2); %Date
Pr=Pr(x1:x2); %Precipitation
Ta=Ta(x1:x2); %Air temperature
Ws=Ws(x1:x2); %Wind speed
Ws(Ws<=0)=0.01; %Wind speed equal to 0.01 if it is lower than zero
ea=ea(x1:x2); %Vapor pressure
SAD1=SAD1(x1:x2); %First band diffuse radiation
SAD2=SAD2(x1:x2); %Second band diffuse radiation
SAB1=SAB1(x1:x2); %First band direct radiation 
Pre=Pre(x1:x2); %Atmospheric pressure
SAB2=SAB2(x1:x2); %First band direct radiation
N=N(x1:x2); %Cloud cover or longwave incoming radiation 
Tdew=Tdew(x1:x2);
esat=esat(x1:x2); %Vapor pressure at saturation
PARB=PARB(x1:x2); %PAR Radiation direct
PARD = PARD(x1:x2); %PAR radiation diffuse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
t_bef= -0.67; %Integration interval for solar variables
t_aft= 1.67; %Integration interval for solar variables
%%%%%%%%%%%%%%%%%%%%
Ds=esat-ea; %% [Pa] Vapor Pressure Deficit
Ds(Ds<0)=0; %%All negative values are set to zero

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Carbon data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('M:\19_ISTA\1_TC\3_Model_Source\TeC_Source_Code\Inputs\Ca_Data.mat'); %CO2 atmospheric concentration (raw data)
d1 = find(abs(Date_CO2-Date(1))<1/36);d2 = find(abs(Date_CO2-Date(end))<1/36);
Ca=Ca(d1:d2); %CO2 atmospheric concentration
clear d1 d2 Date_CO2 
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%Time management
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%dt,Pr(i),Ta(i),Ws(i),ea(i),Pre(i),Rdir(i),Rdif(i),N(i),z,Tdew(i),esat(i),.
[YE,MO,DA,HO,MI,SE] = datevec(Date); %Date separated in Year(YE)-Month(MO)-Day(DA)-
                                     %Hour(HO)-Minute(MI)-Second(SE)
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%DIRECTORIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAM_IC = strcat(current_directory,'\MOD_PARAM_',id_location);
%PARAM_IC = strcat(current_directory,'/MOD_PARAM_',id_location);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Directory = uigetdir('Window','Insert Directory Noname Package') ;
Directory='M:\19_ISTA\1_TC\3_Model_Source\TeC_Source_Code\T&C_Code';
cd(Directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Calling and setting the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAIN_FRAME;
cd(current_directory);