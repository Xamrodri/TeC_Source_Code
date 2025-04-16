%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

PARAMETER FILE TEMPLATE FOR T&C 

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DIRECTORY
cur_dir=cd;
%cd(Directory)

%% Basic Grid Cell Information
% Not relevant for plot scale 

fpr=1;
%SvF=1; %% Sky View Factor
SN=0; % Stream Identifier
%Slo_top=0;  % Slope [fraction dy/dx]
Slo_pot=zeros(1,ms); % Slope of hydraulic head [fraction dy/dx]
Asur = 1./cos(atan(Slo_top)); % Real Area/Projected Area [m^2/m^2]
Ared = 1; % Reduction due to soil rock content 
aR =0; % anisotropy ratio %Kh=Ks*aR;
% cellsize=1; %%[m^2];
aTop = 1000*cellsize^2./cellsize; %  Ratio betweeen Area/ContourLenght [mm]

%% Rainfall disaggregation information 
a_dis = NaN ;
pow_dis = NaN;

%% SOIL PARAMETERS
%==========================================================================
%
%==========================================================================
Kbot = [0.0 0.0 0.0 0.0 0.0 0.0 5];     % Conductivity at the bedrock layer [mm/h] 
Krock = [NaN NaN NaN NaN NaN NaN 0.15]; % Conductivity of Fractured Rock [mm/h] 

Kbot= Kbot(ksv(ij)); % Conductivity of the bedrock [mm/h] 
Krock=Krock(ksv(ij)); % Hydraulic conductivity fractured rock [mm/h]

% Soil layers
Zs= [0    10    20    50   100   150   200   300   400    700   1000]; % Depth of top of the soil layer [mm],  ms+1
Zdes = 10; % Depth of evaporation layer [mm]
Zinf=  10; % Depth of infiltration layer (=first layer) [mm]
Zbio = 250; % Depth of the active Biogeochemistry zone [mm]

if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end

[EvL_Zs]=Evaporation_layers(Zs,Zdes); % Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); % Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio); % Infiltration Depth Layer fraction
dz= diff(Zs); % Thickness of the Layers [mm]
Dz=zeros(1,ms); % Soil layer thickness [mm]

for i = 1:ms
    if i>1
        Dz(i)= (dz(i)+ dz(i-1))/2; % Delta Depth Between Middle Layer [mm]
    else
        Dz(i)=dz(1)/2; % Delta Depth Between First Middle Layer and soil surface [mm]
    end
end

%% SOIL PARAMETERS 
% Pcla= 0.2;
% Psan= 0.4;
% Porg= 0.025;
Color_Class = 0;
%%%%%%%%%%%%%%%%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%%%%%%%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%
SPAR=2; %%% SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%%%%%%%%
Osat=Osat*ones(1,ms);
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]

%%%%%%%%%%%%%%%% Define soil layer depth by specifying 
%%%%%%%%%%%%%%%% impermeable layer at given layer
Soil_th=Zs(end)*(SOIL_TH(ij)/100); %% convert from relative to absolute depth
[~, ix]= min(abs(Zs-Soil_th));
if ix==1
    ix=2;
end
Ks_Zs(ix-1:end)=10^4; 
clear ix


%% Soil parameter function

% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]

% Function
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp
Ofc(Ks_Zs<0.2)=Osat(Ks_Zs<0.2);
Oice = 0;


%% SNOW AND ICE PARAMETERS
%==========================================================================
%
%==========================================================================

TminS=-1.1; %% Threshold temperature snow
TmaxS= 2.5; %% Threshold temperature snow
ros_max1 = 580; %520; %600; %%% [kg/m^3]
ros_max2 = 300; %320; %450; %%% [kg/m^3]
Th_Pr_sno = 8; %%% [mm/h] Threshold Intensity of snow to consider a New SnowFall

%========================= ICE Parameter =============================

Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
WatFreez_Th = -8; %% [ï¿½C] Threshold for freezing water
dz_ice = 0.45; %% [mm / h] Water Freezing Layer progression without snow-layer  

%======================= Albedo for snow/ice ========================

% Ameas=rand(NN,1); %Making dummy variable here, but add measured if you have it, it will only be used if one of the switches is on. 
% Any missing data leave as NaN and modelled Albedo will be used.
%Aice_meas_on_hourly = NaN(NN,1);
%Asno_meas_on_hourly = NaN(NN,1);

%idxa=isnan(Ameas)==1;
%Aice_meas_on_hourly(idxa)=0;  %Use modelled ice albedo if not measured
%Aice_meas_on_hourly(~idxa)=alpha; %Switch; When = 1 use measured ice albedo when available (Change to 0 if you don't want to use measured albedo at all)
%Asno_meas_on_hourly(idxa)=0;  %Use modelled snow albedo if not measured
%Asno_meas_on_hourly(~idxa)=alpha; %Switch; When = 1 use measured snow albedo when available (Change to 0 if you don't want to use measured albedo at all)

Aice = 0.28; %% Default ice albedo (needed for Restating_parameters_temporary until in loop)

if exist('Afirn','var')
     Aice=Afirn(ij); %Use landsat distributed measured albedo
end 

%% Debris Cover Glacier
%==========================================================================
% Debris Cover Glacier
% MAKE SURE DEBRIS THICKNESS 0 FOR CLEAN GLACIER!!
%==========================================================================

albs   =     [0.153  0.115     0.13     0.13   0.13     0.13 ]; % Ask Achille??
lans   =     [1.65   0.985     1.45     0.94   0.94     0.94 ];
zoms   =     [0.38   0.081     0.15     0.016  0.016    0.016];

%Deb_Par.alb= albs(s);
%Deb_Par.e_sur =  0.94;
%Deb_Par.lan = lans(s);
%Deb_Par.rho = 1496;  % [kg/m^3]
%Deb_Par.cs = 948;   % [J/kg K]
%Deb_Par.zom = zoms(s);

dbThick=DEB_MAP(ij);%% [mm]

%%  ROOT PARAMETER 

%categories  [fir    larch     grass    shrub  BLever   BLdec]  
ExEM   =     [0.0    0.0       0.0      0.0    0.0      0.0  ]; %Fraction of EM mycorrhizal on the toal mycorrhizal [-]
% Selection based on II
ExEM = ExEM(II);

CASE_ROOT= 1;  % Type of Root Profile
%categories  [fir    larch     grass    shrub  BLever   BLdec]  
ZR95_H   =   [800    800       0        0      1000     800]; % Root depth 95 percentile, high vegetation [mm]
ZR95_L   =   [0      0         200      600    0        0  ]; % Root depth 95 percentile, low vegetation [mm]
ZR50_H   =   [NaN    NaN       NaN      NaN    NaN      NaN]; % Root depth 50 percentile, high vegetation [mm]
ZR50_L   =   [NaN    NaN       NaN      NaN    NaN      NaN]; % Root depth 50 percentile, low vegetation [mm]
ZRmax_H  =   [NaN    NaN       NaN      NaN    NaN      NaN]; % Maximum root depth, high vegetation [mm]
ZRmax_L  =   [NaN    NaN       NaN      NaN    NaN      NaN]; % Maximum root depth, low vegetation [mm]

% Selection based on II
ZR95_H =ZR95_H(II); ZR50_H =ZR50_H(II); ZRmax_H =ZRmax_H(II);
ZR95_L =ZR95_L(II); ZR50_L =ZR50_L(II); ZRmax_L =ZRmax_L(II);

% If element ij is rock, then:
if ksv(ij) == 7 
    ZR95_H = [0];
    ZR95_L = [0]; 
    ZR50_H = [NaN];
    ZR50_L = [NaN];
    ZRmax_H = [NaN];
    ZRmax_L = [NaN];
end

%%  VEGETATION SECTION
%==========================================================================
%{
Classes in T&C:
        1 Fir (evergr.)
        2 Larch (decid.)
        3 Grass C3
        4 Shrub (decid.)
        5 Broadleaf evergreen
        6 Broadleaf deciduous
        7 Rock
%}  
%==========================================================================

%% INTERCEPTION PARAMETERS 

In_max_urb= 5;    % Maximum interception capacity in urban [mm]
In_max_rock= 0.1; % Maximum interception capacity in rocks [mm]
Kct=0.75;         % Foliage cover decay factor for throughfall [-]

%% Interception Parameter
gcI=3.7;       % Interception parameter [1/mm]
KcI=0.06;      % Interception drainage rate coefficient [mm] - Mahfouf and Jacquemin 1989
Sp_SN_In= 5.9; % Specific interception of rainfall for unit leaf area. Average of high vegetation [mm/LAI]

%categories  [fir    larch     grass    shrub  BLever   BLdec]  
Sp_LAI_H_In= [0.1     0.1      0.2      0.2    0.2      0.2 ]; % Specific interception of rainfall for unit leaf area, high vegetation [mm/LAI]
Sp_LAI_L_In= [0.2     0.2      0.2      0.1    0.2      0.2 ]; % Specific interception of rainfall for unit leaf area, low vegetation [mm/LAI]

% Selection based on II
Sp_LAI_H_In =Sp_LAI_H_In(II);
Sp_LAI_L_In =Sp_LAI_L_In(II);

%% Leaf Dimension

%categories  [fir    larch     grass    shrub  BLever   BLdec]  
d_leaf_H =   [0.25   0.8       2        2      5        4 ]; % Leaf characteristic dimension, high vegetation [cm]
d_leaf_L =   [2      2         0.8      3      2        2];  % Leaf characteristic dimension, low vegetation [cm]

% Selection based on II
d_leaf_H = d_leaf_H(II);
d_leaf_L =d_leaf_L(II);

%% Veg Biochemical parameter

%categories  [fir    larch     grass    shrub  BLever   BLdec]  
KnitH      = [0.35   0.2       NaN      NaN    0.35     0.30 ]; % Canopy Nitrogen Decay coefficient, high vegetation [-]
KnitL      = [NaN    NaN       0.15     0.25   NaN      NaN  ]; % Canopy Nitrogen Decay coefficient, low vegetation [-]

mSl_H      = [0      0         NaN      NaN    0        0    ]; % Linear coefficient of increasing specific leaf area with LAI, high vegetation [m2 PFT /gC] 
mSl_L      = [NaN    NaN       0        0      NaN      NaN  ]; % Linear coefficient of increasing specific leaf area with LAI, low vegetation  [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree

% Selection based on II
KnitH =KnitH(II); 
KnitL =KnitL(II);
mSl_H =mSl_H(II); 
mSl_L =mSl_L(II);

%%  Photosynthesis Parameter

% High vegetation
%categories  [fir    larch     grass    shrub  BLever   BLdec] 
FI_H   =     [0.081  0.081     NaN      NaN    0.081    0.081 ]; % Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H   =     [800    700       NaN      NaN    800      1000  ]; % Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis [Pa] 
a1_H   =     [6      6         NaN      NaN    7        6     ]; % WUE parameter [-]  
go_H   =     [0.01   0.01      NaN      NaN    0.01     0.01  ]; % minimum Stomatal Conductance [mol / s m^2] 
CT_H   =     [3      3         NaN      NaN    3        3     ]; % Photosyntesis pathway - Typology for Plants --> C3 or C4 
DSE_H  =     [0.649  0.66      NaN      NaN    0.649    0.649 ]; % Activation Energy - Plant Dependent [kJ/mol] 
Ha_H   =     [72     94        NaN      NaN    72       76    ]; % Entropy factor - Plant Dependent [kJ / mol K]  
gmes_H =     [Inf    Inf       NaN      NaN    Inf      Inf   ]; % Mesophyll conductance [mol CO2 / s m^2 ];  
rjv_H  =     [2      1.5       NaN      NaN    2        2     ]; % Ratio Jmax - Vmax  [umol electrons / umolCO2 ]

% Selection based on II
FI_H=FI_H(II); Do_H=Do_H(II); a1_H=a1_H(II); go_H=go_H(II);
CT_H=CT_H(II); DSE_H=DSE_H(II); Ha_H=Ha_H(II); gmes_H=gmes_H(II);
rjv_H=rjv_H(II);

% Low vegetation
%categories  [fir    larch     grass    shrub  BLever   BLdec] 
FI_L   =     [NaN    NaN       0.081    0.081  NaN      NaN  ]; % Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L   =     [NaN    NaN       1000     1000   NaN      NaN  ]; % [Pa] 
a1_L   =     [NaN    NaN       8        7      NaN      NaN  ]; % [-] WUE parameter 
go_L   =     [NaN    NaN       0.01     0.01   NaN      NaN  ]; % [mol / s m^2] minimum Stomatal Conductance
CT_L   =     [NaN    NaN       3        3      NaN      NaN  ]; %--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L  =     [NaN    NaN       0.649    0.649  NaN      NaN  ]; % [kJ/mol] Activation Energy - Plant Dependent
Ha_L   =     [NaN    NaN       62       72     NaN      NaN  ]; % [kJ / mol K]  entropy factor - Plant Dependent
gmes_L =     [NaN    NaN       Inf      Inf    NaN      NaN  ]; % [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L  =     [NaN    NaN       1.9      2.2    NaN      NaN  ]; % Ratio Jmax - Vmax  [umol electrons / umolCO2 ]

% Selection based on II
FI_L=FI_L(II); Do_L=Do_L(II); a1_L=a1_L(II); go_L=go_L(II);
CT_L=CT_L(II); DSE_L=DSE_L(II); Ha_L=Ha_L(II); gmes_L=gmes_L(II);
rjv_L=rjv_L(II);

%categories  [fir    larch     grass    shrub  BLever   BLdec] 
Vmax_H  =    [45     64        0        0      32       48   ]; % Maximum Rubisco Capacity [umol CO2 /m2 s]  
Vmax_L  =    [0      0         50       46     0        0    ]; % Maximum Rubisco Capacity [umol CO2 /m2 s]

% Selection based on II
Vmax_H  =   Vmax_H(II); Vmax_L =Vmax_L(II);


%% Hydraulic Parameters
%categories    [fir    larch   grass    shrub  BLever   BLdec] 
Psi_sto_00_H = [-0.8   -0.8    NaN      NaN    -1.0     -0.8 ]; % Water Potential at 2% loss conductivity [MPa] 
Psi_sto_50_H = [-2.5   -2.5    NaN      NaN    -2.8     -2.5 ]; % Water Potential at 50% loss conductivity [MPa]  

% Leaf
PsiL00_H     = [-1     -1      NaN      NaN    -1.2     -1.0 ]; % Water Potential at 2% loss conductivity [MPa] 
PsiL50_H     = [-3.2   -3.2    NaN      NaN    -4.0     -3.0 ]; % Water Potential at 50% loss conductivity [MPa]  
Kleaf_max_H  = [10     10      NaN      NaN    10       10   ]; % Leaf maximum hydraulic conductivity [mmolH20 m^2 leaf s /MPa]
Cl_H         = [1200   1200    NaN      NaN    1200     1200 ]; % Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]

% Xylem
Axyl_H       = [15     15      NaN      NaN    15       15   ]; % [cm^2 stem /m^2 PFT]
Kx_max_H     = [80000  80000   NaN      NaN    80000    80000]; % Xylem Conductivity specific for water. 5550-555550 [mmolH20 /m s MPa]  
PsiX50_H     = [-5     -5      NaN      NaN    -6       -4.5 ]; % Water Potential at 50% loss conductivity [MPa]
Cx_H         = [150    150     NaN      NaN    150      150  ]; % [kg / m^3 sapwood MPa]

% Stomata
Psi_sto_00_L = [NaN    NaN    -0.5      -1     NaN      NaN  ]; % Water Potential at 2% loss conductivity  [MPa]  
Psi_sto_50_L = [NaN    NaN    -2.8      -3     NaN      NaN  ]; % Water Potential at 50% loss conductivity [MPa]  

% Leaf
PsiL00_L     = [NaN    NaN    -1        -2.5   NaN      NaN  ]; % Water Potential at 2% loss conductivity [MPa]  
PsiL50_L     = [NaN    NaN    -3.5      -4.5   NaN      NaN  ]; % Water Potential at 50% loss conductivity [MPa]  
Kleaf_max_L  = [NaN    NaN    5         5      NaN      NaN  ]; % [mmolH20 m^2 leaf s /MPa]
Cl_L         = [NaN    NaN    1200      1200   NaN      NaN  ]; % Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]

% Xylem
Axyl_L       = [NaN    NaN    0         0      NaN      NaN  ]; % Xylem area over PFT area [cm^2 stem /m^2 PFT]
Kx_max_L     = [NaN    NaN    80000     80000  NaN      NaN  ]; % Xylem Conductivity specific for water;5550-555550 [mmolH20 /m s MPa]  
PsiX50_L     = [NaN    NaN    -4.5      -9     NaN      NaN  ]; % Water Potential at 50% loss conductivity[MPa] 
Cx_L         = [NaN    NaN    150       150    NaN      NaN  ]; % Steam capacitance low vegetation [kg / m^3 sapwood MPa]

% Selection based on II
Psi_sto_50_H =Psi_sto_50_H(II);  Psi_sto_00_H =Psi_sto_00_H(II);
PsiL00_H = PsiL00_H(II); PsiL50_H=PsiL50_H(II);  Kleaf_max_H=Kleaf_max_H(II);
Cl_H=Cl_H(II); Axyl_H=Axyl_H(II); Kx_max_H=Kx_max_H(II); PsiX50_H=PsiX50_H(II); Cx_H=Cx_H(II);
Psi_sto_50_L =Psi_sto_50_L(II);  Psi_sto_00_L =Psi_sto_00_L(II);
PsiL00_L = PsiL00_L(II); PsiL50_L=PsiL50_L(II);  Kleaf_max_L=Kleaf_max_L(II);
Cl_L=Cl_L(II); Axyl_L=Axyl_L(II); Kx_max_L=Kx_max_L(II); PsiX50_L=PsiX50_L(II); Cx_L=Cx_L(II);

%% Root Parameters (Function)
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

%% Growth Parameters
% High vegetation
%categories    [fir    larch   grass    shrub  BLever   BLdec] 
PsiG50_H   =   [-0.5   -0.8    NaN      NaN    NaN      NaN  ]; % Water potential at 50% impairment of growth and allocation control [MPa]
PsiG99_H   =   [-2.5   -2.5    NaN      NaN    NaN      NaN  ]; % Water potential at 90% impairment of growth and allocation control [MPa]
gcoef_H    =   [3.5    3.5     NaN      NaN    NaN      NaN  ]; % Parameter for maximum growth in perfect conditions, related to Env. controls of growth [gC/m2 day]

% Low vegetation
PsiG50_L   =   [NaN    NaN     -2.8     -3     NaN      NaN  ]; % Water potential at 50% impairment of growth and allocation control [MPa]
PsiG99_L   =   [NaN    NaN     -4       -4.5   NaN      NaN  ]; % Water potential at 90% impairment of growth and allocation control [MPa]
gcoef_L    =   [NaN    NaN     3.5      3.5    NaN      NaN  ]; % Parameter for maximum growth in perfect conditions, related to Env. controls of growth[gC/m2 day]

% Selection based on II
PsiG50_H=PsiG50_H(II); PsiG99_H=PsiG99_H(II); gcoef_H=gcoef_H(II);
PsiG50_L=PsiG50_L(II); PsiG99_L=PsiG99_L(II); gcoef_L=gcoef_L(II);

%% Vegetation Optical Parameter
%categories    [fir    larch   grass    shrub  BLever   BLdec] 
OPT_PROP_H  =  [2 3 0 0 5 7];   % =PFT_Class for "Veg_Optical_Parameter"-function
OPT_PROP_L  =  [0 0 13 2 0 0];

% Selection based on II
OPT_PROP_H = OPT_PROP_H(II);
OPT_PROP_L = OPT_PROP_L(II);

for i=1:cc
    %%%%%%%% Vegetation Optical Parameter
    [PFT_opt_H(i)]=Veg_Optical_Parameter(OPT_PROP_H(i));
    [PFT_opt_L(i)]=Veg_Optical_Parameter(OPT_PROP_L(i));
end

%categories    [fir    larch   grass    shrub  BLever   BLdec] 
OM_H  =        [1      1       NaN      NaN    1        1    ]; % Within canopy clumping factor [-]
OM_L  =        [NaN    NaN     1        1      NaN      NaN  ]; % Within canopy clumping factor [-]

% Selection based on II
OM_H=OM_H(II); OM_L=OM_L(II);

%% Specific leaf area of litter
Sllit = 2 ; % Litter Specific Leaf area [m2 Litter / kg DM]

%% High Vegetation
%BLever: broadleaf evergreen vegetation dec;
%BLdec: broadleaf deciduous vegetation dec

%categories   [fir   larch      grass  shrub  BLever       BLdec]  
aSE_H    =    [0        1       NaN    NaN    0            1    ]; % Allocation to reserve carbohydrate Values: 1 for Seasonal Plant and 0 for Evergreen
Sl_H     =    [0.010    0.025   NaN    NaN    0.016        0.020]; % Specific leaf area of  biomass [m^2 /gC]. Values: 0.05 -0.005
Nl_H     =    [42       26      NaN    NaN    40           28   ]; % Leaf Nitrogen Concentration [kgC/kgN ] 
r_H      =    [0.058    0.055   NaN    NaN    0.045        0.035]; % respiration rate at 10° [gC/gN d ]. Values: [0.066 -0.011]
gR_H     =    [0.25     0.25    NaN    NaN    0.25         0.25 ]; % Growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
Tcold_H  =    [-40      -3      NaN    NaN    3            3.5  ]; % Cold Leaf Shed [°C]
age_cr_H =    [950      180     NaN    NaN    365          180  ]; % Critical Leaf Age [day]
Trr_H    =    [0.25     3       NaN    NaN    0.5          3.5  ]; % Translocation rate [gC /m^2 d]
LtR_H    =    [0.8      0.8     NaN    NaN    1.0          0.9  ]; % Leaf to Root ratio maximum
eps_ac_H =    [0.25     1       NaN    NaN    0.5          1.0  ]; % Allocation to reserve parameter [0-1]
fab_H    =    [0.74     0.8     NaN    NaN    0.74         0.74 ]; % fraction above-ground sapwood and reserve
ff_r_H   =    [0.1      0.1     NaN    NaN    0.1          0.1  ]; % Reference allocation to Fruit and reproduction
Wm_H     =    [0        0       NaN    NaN    0            0    ]; % wood turnover coefficient [1/d]

dd_max_H = 1./[150      200     NaN    NaN   200           100  ]; % Death maximum for drought [1/d] 
dc_C_H   = 1./[5        10      NaN    NaN   365           10   ]; % Factor of increasing mortality for cold
drn_H    = 1./[900      1100    NaN    NaN   550           800  ]; % Turnover root  [1/d]
dsn_H    = 1./[1100     750     NaN    NaN   800           700  ]; % Normal transfer rate sapwood [1/d]
Mf_H     = 1./[80       50      NaN    NaN   80            80   ]; % Fruit maturation turnover [1/d]
Klf_H    = 1./[40       30      NaN    NaN   30            28   ]; % Dead Leaves fall turnover [1/d]

fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
% [Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
% [Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
% [Stoich_H(3)]=Veg_Stoichiometric_Parameter(Nl_H(3));
% [Stoich_H(4)]=Veg_Stoichiometric_Parameter(Nl_H(4));
% [Stoich_H(5)]=Veg_Stoichiometric_Parameter(Nl_H(5));
% [Stoich_H(6)]=Veg_Stoichiometric_Parameter(Nl_H(6));

%% Phenology 
%categories   [fir     larch    grass  shrub  BLever    BLdec]  
Bfac_lo_H =   [0.99    0.99     NaN    NaN    0.95      0.95 ]; % Leaf Onset Water Stress
Bfac_ls_H =   [NaN     NaN      NaN    NaN    NaN       NaN  ]; % Not-used 
Tlo_H     =   [4.5     3.5      NaN    NaN    5.5       2.8  ]; % Mean Temperature for Leaf onset
Tls_H     =   [NaN     NaN      NaN    NaN    NaN       NaN  ]; % Not-used 
PAR_th_H  =   [NaN     NaN      NaN    NaN    NaN       NaN  ]; % Light Phenology Threshold 
dmg_H     =   [30      30       NaN    NaN    45        30   ]; % Day of Max Growth
LAI_min_H =   [0.001   0.01     NaN    NaN    0.001     0.01 ];
mjDay_H   =   [220     250      NaN    NaN    250       250  ]; % Maximum Julian day for leaf onset
LDay_min_H =  [12.8    12.7     NaN    NaN    12.1      11.7 ]; % Minimum Day duration for leaf onset
LDay_cr_H =   [11.8    11.6     NaN    NaN    11.8      12.0 ]; % Threshold for senescence day light [h]

% Selection based on II
Sl_H =Sl_H(II); Nl_H=Nl_H(II);
r_H=r_H(II); gR_H=gR_H(II); aSE_H=aSE_H(II); dd_max_H=dd_max_H(II);
dc_C_H=dc_C_H(II); Tcold_H=Tcold_H(II); drn_H=drn_H(II);
dsn_H=dsn_H(II);  age_cr_H=age_cr_H(II);
Bfac_lo_H=Bfac_lo_H(II); Bfac_ls_H=Bfac_ls_H(II);
Tlo_H = Tlo_H(II);  Tls_H=Tls_H(II);
dmg_H = dmg_H(II); LAI_min_H=LAI_min_H(II);
Trr_H = Trr_H(II);  mjDay_H=mjDay_H(II);
LDay_min_H= LDay_min_H(II); LtR_H =LtR_H(II);
Mf_H= Mf_H(II);  Wm_H= Wm_H(II);  eps_ac_H = eps_ac_H(II);
LDay_cr_H = LDay_cr_H(II);  Klf_H = Klf_H(II);
fab_H = fab_H(II); fbe_H = fbe_H(II); ff_r_H = ff_r_H(II);


for i=1:cc
    [Stoich_H(i)]=Veg_Stoichiometric_Parameter(Nl_H(i));
    [ParEx_H(i)]=Exudation_Parameter(0);
    [Mpar_H(i)]=Vegetation_Management_Parameter;
end

%% Low Vegetation 
%categories   [fir     larch    grass  shrub    BLever    BLdec]  
aSE_L   =     [NaN     NaN      2      1        NaN       NaN]; % Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_L    =     [NaN     NaN      0.016  0.018    NaN       NaN]; % Specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_L    =     [NaN     NaN      23     35       NaN       NaN]; % Leaf Nitrogen Concentration [kgC/kgN ]
r_L     =     [NaN     NaN      0.025  0.025    NaN       NaN]; % Respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_L    =     [NaN     NaN      0.25   0.25     NaN       NaN]; % Growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
Tcold_L =     [NaN     NaN      3      2        NaN       NaN]; % Cold Leaf Shed [°C] 


dd_max_L = 1./[NaN     NaN      45     100      NaN       NaN]; % Death maximum for drought [1/d]
dc_C_L   = 1./[NaN     NaN      32     15       NaN       NaN]; % Factor of increasing mortality for cold
drn_L    = 1./[NaN     NaN      950    1600     NaN       NaN]; % turnover root  [1/d]
dsn_L    = 1./[NaN     NaN      365    1800     NaN       NaN]; % normal transfer rate sapwood [1/d]
Mf_L     = 1./[NaN     NaN      50     50       NaN       NaN]; % Fruit maturation turnover [1/d]
Klf_L    = 1./[NaN     NaN      20     40       NaN       NaN]; % Dead Leaves fall turnover [1/d]

age_cr_L =    [NaN     NaN      180    180      NaN       NaN]; % [day] Critical Leaf Age
Trr_L    =    [NaN     NaN      2      0.6      NaN       NaN]; % Translocation rate [gC /m^2 d]
LtR_L    =    [NaN     NaN      0.7    1        NaN       NaN]; % Leaf to Root ratio maximum
Wm_L     =    [NaN     NaN      0      0        NaN       NaN] ;% wood turnover coefficient [1/d]
eps_ac_L =    [NaN     NaN      1      1        NaN       NaN]; % Allocation to reserve parameter [0-1]
fab_L    =    [NaN     NaN      0      0.75     NaN       NaN]; % fraction above-ground sapwood and reserve
fbe_L    = 1-fab_L;                          % fraction below-ground sapwood and reserve
ff_r_L   =    [NaN     NaN      0.1    0.1      NaN       NaN]; % Reference allocation to Fruit and reproduction

% Phenology 
Bfac_lo_L  =  [NaN     NaN      0.99   0.99     NaN       NaN]; % Leaf Onset Water Stress
Bfac_ls_L  =  [NaN     NaN      0.15   NaN      NaN       NaN] ;% 
Tlo_L      =  [NaN     NaN      2.5    2.5      NaN       NaN]; % Mean Temperature for Leaf onset
Tls_L      =  [NaN     NaN      NaN    NaN      NaN       NaN]; % Not-used 
PAR_th_L   =  [NaN     NaN      NaN    NaN      NaN       NaN]; % Light Phenology Threshold 
dmg_L      =  [NaN     NaN      20     25       NaN       NaN]; % Day of Max Growth
LAI_min_L  =  [NaN     NaN      0.05   0.001    NaN       NaN];
mjDay_L    =  [NaN     NaN      250    180      NaN       NaN]; % Maximum Julian day for leaf onset
LDay_min_L =  [NaN     NaN      12     12.2     NaN       NaN]; % Minimum Day duration for leaf onset
LDay_cr_L  =  [NaN     NaN      12     12       NaN       NaN]; % Threshold for senescence day light [h]

% Selection of parameters
Sl_L =Sl_L(II); Nl_L=Nl_L(II);
r_L=r_L(II); gR_L=gR_L(II); aSE_L=aSE_L(II); dd_max_L=dd_max_L(II);
dc_C_L=dc_C_L(II); Tcold_L=Tcold_L(II); drn_L=drn_L(II);
dsn_L=dsn_L(II);  age_cr_L=age_cr_L(II);
Bfac_lo_L=Bfac_lo_L(II); Bfac_ls_L=Bfac_ls_L(II);
Tlo_L = Tlo_L(II);  Tls_L=Tls_L(II);
dmg_L = dmg_L(II); LAI_min_L=LAI_min_L(II);
Trr_L = Trr_L(II);  mjDay_L=mjDay_L(II);
LDay_min_L= LDay_min_L(II); LtR_L =LtR_L(II);
Mf_L= Mf_L(II);  Wm_L= Wm_L(II);  eps_ac_L = eps_ac_L(II);
LDay_cr_L = LDay_cr_L(II);  Klf_L = Klf_L(II);
fab_L = fab_L(II); fbe_L = fbe_L(II); ff_r_L = ff_r_L(II);

for i=1:cc
    [Stoich_L(i)]=Veg_Stoichiometric_Parameter(Nl_L(i));
    [ParEx_L(i)]=Exudation_Parameter(0);
    [Mpar_L(i)]=Vegetation_Management_Parameter;
end

%% INITIAL CONDITIONS
%==========================================================================
%
%==========================================================================

L_day=zeros(NNd,1);
for j=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')

if OPT_SoilBiogeochemistry == 1
    %%%%%
end
%%%

%% HYDROLOGY
%==========================================================================
% Initial conditions Hydrology/non-Vegetation
%==========================================================================
SWE(1)=SWEtm1(ij); %% [mm]
SND(1)=SNDtm1(ij); %% [m]
Ts(1)=Ta(1)+2;
Tdamp(1)=Ta(1);
Tdp(1,:)= Ta(1)*ones(1,ms);

%%% Snow_alb = soil_alb initial
snow_alb.dir_vis = 0.6;
snow_alb.dif_vis = 0.6;
snow_alb.dir_nir = 0.6;
snow_alb.dif_nir = 0.6;

In_L(1,:)=In_Ltm1(ij); In_H(1,:)=In_Htm1(ij);
In_urb(1)=In_urbtm1(ij); In_rock(1)= In_rocktm1(ij);
In_Litter(1)=In_Littertm1(ij);
SP_wc(1)=SP_wctm1(ij) ; %%[mm]
In_SWE(1)= In_SWEtm1(ij);
ros(1)= rostm1(ij);
t_sls(1)= t_slstm1(ij);
e_sno(1) = e_snotm1(ij);
tau_sno(1) = tau_snotm1(ij);
EK(1)=EKtm1(ij);
WAT(1) = WATtm1(ij);% 
ICE(1) = ICEtm1(ij);% 
IP_wc(1)= IP_wctm1(ij);
ICE_D(1)= ICE_Dtm1(ij);%   ; 
FROCK(1)=FROCKtm1(ij);
Ws_under(1)=Ws_undertm1(ij); 
%%%%%%%%%%%%%% Volume [mm]
O(1,:)= Ofc;
%%%%%%%%%%%%%%%%%%%
V(1,:) = (O(1,:) -Ohy).*dz.*0; % TRYING WITH 0 initial water content in soil
cd(cur_dir)
%%%%%%%%%%%%%%%%%

%% VEGETATION
%==========================================================================
% Initial conditions Vegetation 
%==========================================================================

Ci_sunL(1,:) = [Ca(1)]; % [umolCO2/mol]
Ci_sunH(1,:) = [Ca(1)]; % [umolCO2/mol]
Ci_shdL(1,:) = [Ca(1)]; % [umolCO2/mol]
Ci_shdH(1,:) = [Ca(1)]; % [umolCO2/mol]
%%%%%%%%%%%%%%%%%%
%%% B1 Leaves - Grass  %%% B2 Sapwood  %%% B3 Fine Root  %%% B4 Carbohydrate Reserve
%%% B5 Fruit and Flower %%% B6 Heartwood - Dead Sapwood %%% B7 Leaves - Grass -- Standing Dead
%%%%%%%%%%%%%%%%%%
LAI_H(1,:) = LAI_Htm1(ij); B_H(1,:,:)= B_Htm1(ij,:,:); Rrootl_H(1,:)= Rrootl_Htm1(ij); 
PHE_S_H(1,:)= PHE_S_Htm1(ij); dflo_H(1,:)=dflo_Htm1(ij); AgeL_H(1,:)=AgeL_Htm1(ij);
e_rel_H(1,:)=e_rel_Htm1(ij); hc_H(1,:) =hc_Htm1(ij); SAI_H(1,:) = SAI_Htm1(ij);
%%%%%%%%%%%%%%%%%%
LAI_L(1,:) = LAI_Ltm1(ij); B_L(1,:,:)= B_Ltm1(ij,:,:); Rrootl_L(1,:)= Rrootl_Ltm1(ij);
PHE_S_L(1,:)=PHE_S_Ltm1(ij); dflo_L(1,:)=dflo_Ltm1(ij); AgeL_L(1,:)=AgeL_Ltm1(ij);
e_rel_L(1,:)=e_rel_Ltm1(ij); hc_L(1,:) =hc_Ltm1(ij); SAI_L(1,:) =SAI_Ltm1(ij);
%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%%%%%%%%%
RexmyI(1,:)= [0 0 0];
%%%%%%%%%%%%%%%%%%
Nreserve_H(1,:)= Nreserve_Htm1(ij); Preserve_H(1,:)= Preserve_Htm1(ij); Kreserve_H(1,:)= Kreserve_Htm1(ij);
FNC_H(1,:)=FNC_Htm1(ij); NupI_H(1,:,:)= NupI_Htm1(ij,:,:); Nreserve_L(1,:)= Nreserve_Ltm1(ij);
Preserve_L(1,:)=Preserve_Ltm1(ij); Kreserve_L(1,:)=Kreserve_Ltm1(ij); FNC_L(1,:)=FNC_Ltm1(ij); NupI_L(1,:,:)= NupI_Ltm1(ij,:,:);
%%%%%%%%%%%%%%%%%%
TdpI_H(1,:)=TdpI_Htm1(ij); TdpI_L(1,:)=TdpI_Ltm1(ij);