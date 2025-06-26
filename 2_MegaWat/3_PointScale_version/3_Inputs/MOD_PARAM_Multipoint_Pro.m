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
%         Depth1 Depth2 Depth3  Depth4 Depth5 Depth6 Depth7 Depth8 Depth9  Depth10
Kbot  =  [0.0    0.0    0.0     0.0    0.0    0.0    5      0.0    0.0     0.0    ]; % Conductivity at the bedrock layer [mm/h] 
Krock =  [NaN    NaN    NaN     NaN    NaN    NaN    0.15   NaN    NaN     NaN    ]; % Conductivity of Fractured Rock [mm/h] 

Kbot = Kbot(ksv(ij)); % Conductivity of the bedrock [mm/h] 
Krock =Krock(ksv(ij)); % Hydraulic conductivity fractured rock [mm/h]

% Soil layers
%           Depth1 Depth2 Depth3  Depth4 Depth5 Depth6 Depth7 Depth8 Depth9  Depth10
Zs= [0      10     20      50     100    150    200    300    400     700    1500]; % Depth of top of the soil layer [mm],  ms+1
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

% Function Soil_parametersII
% Output variables [Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII
[Ofc,~,~,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
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

%categories  [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A  grass_B     shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
albs   =     [0.153       0.13      0.13       0.13     0.13        0.13      0.13       0.13   0.13         0.13       0.13       NaN  ]; % Ask Achille??
lans   =     [1.65        1.45      1.45       1.45     1.45        1.45      1.45       1.45   0.94         0.94       0.94       NaN  ];
zoms   =     [0.38        0.15      0.15       0.15     0.15        0.15      0.15       0.15   0.016        0.016      0.016      NaN  ];

%Deb_Par.alb= albs(s);
%Deb_Par.e_sur =  0.94;
%Deb_Par.lan = lans(s);
%Deb_Par.rho = 1496;  % [kg/m^3]
%Deb_Par.cs = 948;   % [J/kg K]
%Deb_Par.zom = zoms(s);

dbThick=DEB_MAP(ij);%% [mm]



%%  VEGETATION SECTION
%==========================================================================
%{
Classes in T&C:
        1 Fir (evergr.) high elevation (>700m)
            A) Based on MOD_PARAM_Renon (by Simone Fatichi)
            B) Onclude all the Pines/Firs/Spruce in the mountainous regions (>700 m)
        2 CROPS 
            2.1) Crops_WW: Winter wheat
            2.2) Crops_WB: Winter barley
            2.3) Crops_S: Sunflower
            2.4) Crops_R: Rapeseed


        3 Grass C3_A: Parameters for Festuca nigrescens
             A) Based on ValidAlinya

        4 Grass CR_B: Parameters for Nardetum alpigenum
             B) Based on Monte_Bondone

        5 Shrub (decid.), sclerophyllous vegetation, mediterranean shrublands
            A) In Apennines
                A.1) High altitude: Dominiated by Pino-Juniperetalia. 
                     Pinus and Juniperus genera.
                A.2) Subalpine belt: Dominated by Vaccinium
            B) Based on Noe because of Juniperus.
               Parameters belongs to Juniperus phoenicea
        6 Broadleaf evergreen
        
        7 Broadleaf deciduous low elevation (<700m)
            A) Based on NGreece from Simone
            

        8 Broadleaf deciduous high elevation (>700m)
            A) Based on Collelongo from Simone
            B) beech (Fagus sylvatica)
        9 NoVEG
%}  
%==========================================================================

%%  ROOT PARAMETER 

%categories  [fir_high    Crops_WB  Crops_WB   Crops_S  Crops_R    grass_A  grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
ExEM   =     [1.0         0.0       0.0        0.0      0.0        0.0      0.0        0.0    0.0          0.0        0.0        NaN  ]; % Fraction of ectomychorrizae per area  [-]
% Selection based on II
ExEM = ExEM(II);

CASE_ROOT= 1;  % Type of Root Profile
%categories  [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A  grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
ZR95_H   =   [800         0.0       0.0        0.0      0.0         0.0      0.0        0.0    1000         1200       1500       NaN]; % Root depth 95 percentile, high vegetation [mm]
ZR95_L   =   [0           950       900        900      500         250      250        1000   0            250        0          NaN]; % Root depth 95 percentile, low vegetation [mm]
ZR50_H   =   [NaN         NaN       NaN        NaN      NaN         NaN      NaN        NaN    NaN          NaN        NaN        NaN]; % Root depth 50 percentile, high vegetation [mm]
ZR50_L   =   [NaN         NaN       NaN        NaN      NaN         NaN      NaN        NaN    NaN          NaN        NaN        NaN]; % Root depth 50 percentile, low vegetation [mm]
ZRmax_H  =   [NaN         NaN       NaN        NaN      NaN         NaN      NaN        NaN    NaN          NaN        NaN        NaN]; % Maximum root depth, high vegetation [mm]
ZRmax_L  =   [NaN         NaN       NaN        NaN      NaN         NaN      NaN        NaN    NaN          NaN        NaN        NaN]; % Maximum root depth, low vegetation [mm]

% Selection based on II
ZR95_H =ZR95_H(II); ZR50_H =ZR50_H(II); ZRmax_H =ZRmax_H(II);
ZR95_L =ZR95_L(II); ZR50_L =ZR50_L(II); ZRmax_L =ZRmax_L(II);

% If element ij is rock, then:
%if ksv(ij) == 7 
%    ZR95_H = 0;
%    ZR95_L = 0; 
%    ZR50_H = NaN;
%    ZR50_L = NaN;
%    ZRmax_H = NaN;
%    ZRmax_L = NaN;
%end

%% INTERCEPTION PARAMETERS 

In_max_urb= 5;    % Maximum interception capacity in urban [mm]
In_max_rock= 0.1; % Maximum interception capacity in rocks [mm]
Kct=0.75;         % Foliage cover decay factor for throughfall [-]

%% Interception Parameter
gcI=3.7;       % Interception parameter [1/mm]
KcI=0.06;      % Interception drainage rate coefficient [mm] - Mahfouf and Jacquemin 1989
Sp_SN_In= 5.9; % Specific interception of rainfall for unit leaf area. Average of high vegetation [mm/LAI]

%categories  [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A  grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
Sp_LAI_H_In= [0.1         0.2       0.2        0.2      0.2         0.2      0.2        0.2    0.2          0.15       0.15       NaN]; % Specific interception of rainfall for unit leaf area, high vegetation [mm/LAI]
Sp_LAI_L_In= [0.2         0.2       0.2        0.2      0.2         0.2      0.2        0.2    0.2          0.10       0.2        NaN]; % Specific interception of rainfall for unit leaf area, low vegetation [mm/LAI]

% Selection based on II
Sp_LAI_H_In =Sp_LAI_H_In(II);
Sp_LAI_L_In =Sp_LAI_L_In(II);

%% Leaf Dimension

%categories  [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub  BLever_high   BLdec_low  BLdec_high NoVeg]    
d_leaf_H =   [0.25        NaN       NaN        NaN      NaN         3.5       3.5        0.0    5             2.0        4          NaN  ]; % Leaf characteristic dimension, high vegetation [cm]
d_leaf_L =   [0.8         5.0       5.0        5.0      5.0         0.8       0.8        0.7    2             0.8        2          NaN  ]; % Leaf characteristic dimension, low vegetation [cm]

% Selection based on II
d_leaf_H = d_leaf_H(II);
d_leaf_L =d_leaf_L(II);

%% Veg Biochemical parameter

%categories  [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
KnitH      = [0.35        0.2       0.2        0.2      0.2         0.2       0.2        0.0    0.35         0.40       0.25       NaN  ]; % Canopy Nitrogen Decay coefficient, high vegetation [-]
KnitL      = [0.5         0.15      0.15       0.15     0.15        0.15      0.15       0.15   NaN          0.15       NaN        NaN  ]; % Canopy Nitrogen Decay coefficient, low vegetation [-]

mSl_H      = [0           0.0       0.0        0.0      0.0         0.0       0.0        0.0    0            0.0        0.0        NaN  ]; % Linear coefficient of increasing specific leaf area with LAI, high vegetation [m2 PFT /gC] 
mSl_L      = [0           0.0       0.0        0.0      0.0         0.0       0.0        0.0    NaN          0.0        NaN        NaN  ]; % Linear coefficient of increasing specific leaf area with LAI, low vegetation  [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree

% Selection based on II
KnitH =KnitH(II); 
KnitL =KnitL(II);
mSl_H =mSl_H(II); 
mSl_L =mSl_L(II);

%%  Photosynthesis Parameter

% High vegetation
%categories  [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R    grass_A  grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
FI_H   =     [0.081       0.081     0.081      0.081    0.081      0.081    0.081      0.081  0.081        0.081      0.081      NaN  ]; % Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H   =     [1000        1000      1000       1000     1000       1000     1000       1000   800          1000       1000       NaN  ]; % Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis [Pa] 
a1_H   =     [7           7         7          7        7          7        7          7      7            7          6          NaN  ]; % WUE parameter [-]  
go_H   =     [0.01        0.01      0.01       0.01     0.01       0.01     0.01       0.01   0.01         0.01       0.01       NaN  ]; % minimum Stomatal Conductance [mol / s m^2] 
CT_H   =     [3           4         4          4        4          3        3          3      3            3          3          NaN  ]; % Photosyntesis pathway - Typology for Plants --> C3 or C4 
DSE_H  =     [0.649       0.649     0.649      0.649    0.649      0.649    0.649      0.649  0.649        0.649      0.649      NaN  ]; % Activation Energy - Plant Dependent [kJ/mol] 
Ha_H   =     [78          72        72         72       72         72       72         72     72           76         76         NaN  ]; % Entropy factor - Plant Dependent [kJ / mol K]  
gmes_H =     [Inf         Inf       Inf        Inf      Inf        Inf      Inf        Inf    Inf          Inf        Inf        NaN  ]; % Mesophyll conductance [mol CO2 / s m^2 ];  
rjv_H  =     [1.9         1.97      1.97       1.97     1.97       1.97     1.97       2.1    2.0          2.1        2.4        NaN  ]; % Ratio Jmax - Vmax  [umol electrons / umolCO2 ]

% Selection based on II
FI_H=FI_H(II); Do_H=Do_H(II); a1_H=a1_H(II); go_H=go_H(II);
CT_H=CT_H(II); DSE_H=DSE_H(II); Ha_H=Ha_H(II); gmes_H=gmes_H(II);
rjv_H=rjv_H(II);

% Low vegetation
%categories  [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub  BLever_high   BLdec_low   BLdec_high NoVeg]  
FI_L   =     [0.081       0.081     0.081      0.081    0.081       0.081     0.081      0.081  NaN           0.081       NaN        NaN  ]; % Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L   =     [2000        1000      1000       1000     1000        1000      1000       1000   NaN           1000        NaN        NaN  ]; % [Pa] 
a1_L   =     [7           8         7          7        7.5         7         6          5      NaN           8           NaN        NaN  ]; % [-] WUE parameter 
go_L   =     [0.01        0.01      0.01       0.01     0.01        0.01      0.01       0.01   NaN           0.01        NaN        NaN  ]; % [mol / s m^2] minimum Stomatal Conductance
CT_L   =     [3           3         3          3        3           3         3          3      NaN           3           NaN        NaN  ]; %--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L  =     [0.66        0.649     0.649      0.649    0.649       0.649     0.656      0.649  NaN           0.649       NaN        NaN  ]; % [kJ/mol] Activation Energy - Plant Dependent
Ha_L   =     [48          72        72         72       72          72        55         72     NaN           56          NaN        NaN  ]; % [kJ / mol K]  entropy factor - Plant Dependent
gmes_L =     [NaN         Inf       Inf        Inf      Inf         Inf       Inf        Inf    NaN           Inf         NaN        NaN  ]; % [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L  =     [2.6         2.4       2.4        2.4      2.4         1.9       2.4        2.0    NaN           2.2         NaN        NaN  ]; % Ratio Jmax - Vmax  [umol electrons / umolCO2 ]

% Selection based on II
FI_L=FI_L(II); Do_L=Do_L(II); a1_L=a1_L(II); go_L=go_L(II);
CT_L=CT_L(II); DSE_L=DSE_L(II); Ha_L=Ha_L(II); gmes_L=gmes_L(II);
rjv_L=rjv_L(II);

%categories  [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A  grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
Vmax_H  =    [45          0.0       0.0        0.0      0.0         0.0      0.0        0.0    32           44         62         NaN  ]; % Maximum Rubisco Capacity [umol CO2 /m2 s]  
Vmax_L  =    [0           105       105        105      115          64       62         52     0           85         0          NaN  ]; % Maximum Rubisco Capacity [umol CO2 /m2 s]

% Selection based on II
Vmax_H  =   Vmax_H(II); Vmax_L =Vmax_L(II);


%% Hydraulic Parameters
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
Psi_sto_00_H = [-0.5        -0.5      -0.5       -0.5     -0.5        -0.5      -0.5       -1.8   -1.0         -0.5       -0.8       NaN  ]; % Water Potential at 2% loss conductivity [MPa] 
Psi_sto_50_H = [-2.5        -2.0      -2.0       -2.0     -2.0        -2.0      -2.0       -3.5   -2.8         -3.2       -2.5       NaN  ]; % Water Potential at 50% loss conductivity [MPa]  

% Leaf
PsiL00_H     = [-1          -2.5      -2.5       -2.5     -2.7        -2.7      -2.7      -1.3    -1.2         -1.5       -1.2       NaN  ]; % Water Potential at 2% loss conductivity [MPa] 
PsiL50_H     = [-3.2        -3.5      -3.5       -3.5     -5.6        -5.6      -5.6      -4.5    -4.0         -5.2       -3.5       NaN  ]; % Water Potential at 50% loss conductivity [MPa]  
Kleaf_max_H  = [10           5         5          5        5           5         5.0       20      10           10         10        NaN  ]; % Leaf maximum hydraulic conductivity [mmolH20 m^2 leaf s /MPa]
Cl_H         = [1200        1200       1200       1200     1200        1200      1200      1200    1200         1200       1200      NaN  ]; % Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]

% Xylem
Axyl_H       = [15          15         15         15       15          15        15        6.0     15           15         15        NaN  ]; % [cm^2 stem /m^2 PFT]
Kx_max_H     = [80000       80000      80000      80000    80000       80000     80000     80000   80000        80000      80000     NaN  ]; % Xylem Conductivity specific for water. 5550-555550 [mmolH20 /m s MPa]  
PsiX50_H     = [-5          -5        -3.5       -3.5     -3.5        -3.5      -3.5      -6.5    -6           -7         -5.5       NaN  ]; % Water Potential at 50% loss conductivity [MPa]
Cx_H         = [150         150        150        150      150         150       150       80      150          150        150       NaN  ]; % [kg / m^3 sapwood MPa]

% Stomata
Psi_sto_00_L = [-0.8        -0.7      -0.7       -0.7     -0.5        -0.5      -0.5       -0.7    NaN         -0.8        NaN       NaN  ]; % Water Potential at 2% loss conductivity  [MPa]  
Psi_sto_50_L = [-3.0        -2.9      -2.9       -1.1     -2.5        -2.5      -3.0       -5.5    NaN         -3.0        NaN       NaN  ]; % Water Potential at 50% loss conductivity [MPa]  

% Leaf
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A    grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
PsiL00_L     = [-1.1        -1.2      -1.2       -1.2     -1.2        -0.9      -0.5       -1.4      NaN          -1.1       NaN        NaN  ]; % Water Potential at 2% loss conductivity [MPa]  
PsiL50_L     = [-4.0        -3.5      -3.5       -1.5     -3.5        -3.0      -3.0       -7.0      NaN          -4.0       NaN        NaN  ]; % Water Potential at 50% loss conductivity [MPa]  
Kleaf_max_L  = [ 5.0         5.0       5.0        5.0      5.0         5.0       5.0        10.0     NaN           5         NaN        NaN  ]; % [mmolH20 m^2 leaf s /MPa]
Cl_L         = [1200         1200      1200       1200     1200        1200      1200       1200     NaN           1200      NaN        NaN  ]; % Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]

% Xylem
Axyl_L       = [0.0         0.0        0.0        0.0      0.0         0.0        0.0        6.0     NaN           0.0       NaN        NaN  ]; % Xylem area over PFT area [cm^2 stem /m^2 PFT]
Kx_max_L     = [80000       80000      80000      80000    80000       80000      80000      80000   NaN           80000     NaN        NaN  ]; % Xylem Conductivity specific for water;5550-555550 [mmolH20 /m s MPa]  
PsiX50_L     = [-4.5       -4.5       -4.5       -4.5     -4.5        -4.5       -4.5       -9.0     NaN          -4.5       NaN        NaN  ]; % Water Potential at 50% loss conductivity[MPa] 
Cx_L         = [150         150        150        150      150         150        150        80      NaN           150       NaN        NaN  ]; % Steam capacitance low vegetation [kg / m^3 sapwood MPa]

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
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R    grass_A   grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
PsiG50_H   =   [-0.5       -0.45      -0.45      -0.45    -0.45     -0.45     -0.45      -1.4    NaN          -0.5       -0.8        NaN  ]; % Water potential at 50% impairment of growth and allocation control [MPa]
PsiG99_H   =   [-2.5       -1.2       -1.2       -1.2     -1.2      -1.2      -1.2       -1.2    NaN          -3.2       -2.5        NaN  ]; % Water potential at 90% impairment of growth and allocation control [MPa]
gcoef_H    =   [ 3.5        3.5        3.5        3.5      3.5       3.5       3.5        3.5    NaN           4.5        4.5        NaN  ]; % Parameter for maximum growth in perfect conditions, related to Env. controls of growth [gC/m2 day]

% Low vegetation
PsiG50_L   =   [-1.45      -1.2       -1.2       -1.2     -1.2      -2.5       -3.0       -1.4   NaN          -1.45       NaN        NaN  ]; % Water potential at 50% impairment of growth and allocation control [MPa]
PsiG99_L   =   [-4.0       -3.5       -3.5       -1.5     -3.5      -3.0       -4.0       -5.5   NaN          -4.0        NaN        NaN  ]; % Water potential at 90% impairment of growth and allocation control [MPa]
gcoef_L    =   [3.5         3.5        3.5        3.5      3.5       3.5        3.5        3.5   NaN           3.5        NaN        NaN  ]; % Parameter for maximum growth in perfect conditions, related to Env. controls of growth[gC/m2 day]

% Selection based on II
PsiG50_H=PsiG50_H(II); PsiG99_H=PsiG99_H(II); gcoef_H=gcoef_H(II);
PsiG50_L=PsiG50_L(II); PsiG99_L=PsiG99_L(II); gcoef_L=gcoef_L(II);

%% Vegetation Optical Parameter
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
OPT_PROP_H  =  [2           0.0       0.0        0.0      0.0         0.0       0.0        0.0    5            7          7          NaN  ];   % =PFT_Class for "Veg_Optical_Parameter"-function
OPT_PROP_L  =  [0           16        16         16       16          13        13         1      0            13         0          NaN  ];

% Selection based on II
OPT_PROP_H = OPT_PROP_H(II);
OPT_PROP_L = OPT_PROP_L(II);

%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R    grass_A   grass_B    shrub  BLever_high   BLdec_low  BLdec_high NoVeg]  
OM_H  =        [1           1         1          1        1          1         1          1      1             1          1          NaN  ]; % Within canopy clumping factor [-]
OM_L  =        [1           1         1          1        1          1         1          1      NaN           1          NaN        NaN  ]; % Within canopy clumping factor [-]

% Selection based on II
OM_H=OM_H(II); OM_L=OM_L(II);

%% Specific leaf area of litter
Sllit = 2 ; % Litter Specific Leaf area [m2 Litter / kg DM]

%% High Vegetation
%BLever: broadleaf evergreen vegetation dec;
%BLdec: broadleaf deciduous vegetation dec

%categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A  grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
aSE_H    =    [0           1         1          1        1           1        1          1        0            1          1          NaN  ]; % Allocation to reserve carbohydrate Values: 1 for Seasonal Plant and 0 for Evergreen
Sl_H     =    [0.014       0.016     0.016      0.016    0.016       0.016    0.016      0.015    0.016        0.016      0.017      NaN  ]; % Specific leaf area of  biomass [m^2 /gC]. Values: 0.05 -0.005
Nl_H     =    [42          30        30         30       30          30       30         30       40           28         28         NaN  ]; % Leaf Nitrogen Concentration [kgC/kgN ] 
r_H      =    [0.062       0.030     0.030      0.030    0.030       0.030    0.030      0.030    0.045        0.032      0.032      NaN  ]; % respiration rate at 10° [gC/gN d ]. Values: [0.066 -0.011]
gR_H     =    [0.25        0.25      0.25       0.25     0.25        0.25     0.25       0.25     0.25         0.25       0.25       NaN  ]; % Growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
Tcold_H  =    [-20         7         7          7        7           7        7          7        3            6.5        0.0        NaN  ]; % Cold Leaf Shed [°C]
age_cr_H =    [900         150       150        150      150         150      150        150      365          180        150        NaN  ]; % Critical Leaf Age [day]
Trr_H    =    [0.25        3.5       3.5        3.5      3.5         3.5      3.5        3.5      0.5          5.5        5          NaN  ]; % Translocation rate [gC /m^2 d]
LtR_H    =    [0.8         0.8       0.8        0.8      0.8         1.0      1.0        1.0      1.0          0.6        0.5        NaN  ]; % Leaf to Root ratio maximum
eps_ac_H =    [0.20        1.0       1.0        1.0      1.0         1.0      1.0        1.0      0.5          1          1.0        NaN  ]; % Allocation to reserve parameter [0-1]
fab_H    =    [0.74        0.74      0.74       0.74     0.74        0.74     0.74       0.74     0.74         0.74       0.74       NaN  ]; % fraction above-ground sapwood and reserve
ff_r_H   =    [0.1         0.1       0.1        0.1      0.1         0.1      0.1        0.1      0.1          0.1        0.1        NaN  ]; % Reference allocation to Fruit and reproduction
Wm_H     =    [0           0         0          0        0           1/16425  1/16425    0.0      0.0          0          0.0        NaN  ]; % wood turnover coefficient [1/d]

%categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A  grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
dd_max_H = 1./[150         365       365        365      365         365      45         365    200          100        100        NaN  ]; % Death maximum for drought [1/d] 
dc_C_H   = 1./[4.9         182       182        182      182         182      182        182    365          10.2       11         NaN  ]; % Factor of increasing mortality for cold
drn_H    = 1./[900         1095      1095       1095     1095        1095     1095       1095   550          800        1200       NaN  ]; % Turnover root  [1/d]
dsn_H    = 1./[1000        365       365        365      365         365      365        365    800          500        800        NaN  ]; % Normal transfer rate sapwood [1/d]
Mf_H     = 1./[80          50        50         50       50          50       50         0      80           50         80         NaN  ]; % Fruit maturation turnover [1/d]
Klf_H    = 1./[40          15        15         15       15          15       15         15     30           28         28         NaN  ]; % Dead Leaves fall turnover [1/d]

%% check Mf_H in shrub for 1/0

fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
% [Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
% [Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
% [Stoich_H(3)]=Veg_Stoichiometric_Parameter(Nl_H(3));
% [Stoich_H(4)]=Veg_Stoichiometric_Parameter(Nl_H(4));
% [Stoich_H(5)]=Veg_Stoichiometric_Parameter(Nl_H(5));
% [Stoich_H(6)]=Veg_Stoichiometric_Parameter(Nl_H(6));

%% Phenology 
%categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A  grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
Bfac_lo_H =   [0.99        0.95      0.95       0.95     0.95        0.95     0.95       0.95   0.95         0.99       0.95       NaN  ]; % Leaf Onset Water Stress
Bfac_ls_H =   [NaN         NaN       NaN        NaN      NaN         NaN      NaN        NaN    NaN          NaN        NaN        NaN  ]; % Not-used 
Tlo_H     =   [4.5         12.9      12.9       12.9     12.9        12.9     12.9       12.9   5.5          8.8        2.0        NaN  ]; % Mean Temperature for Leaf onset
Tls_H     =   [NaN         NaN       NaN        NaN      NaN         NaN      NaN        NaN    NaN          NaN        NaN        NaN  ]; % Not-used 
PAR_th_H  =   [NaN         NaN       NaN        NaN      NaN         NaN      NaN        NaN    NaN          NaN        NaN        NaN  ]; % Light Phenology Threshold 
dmg_H     =   [30          35        35         35       35          35       35         35     45           30         30         NaN  ]; % Day of Max Growth
LAI_min_H =   [0.001       0.01      0.01       0.01     0.01        0.01     0.01       0.01   0.001        0.05       0.01       NaN  ];
mjDay_H   =   [220         180       180        180      180         180      180        180    250          180        250        NaN  ]; % Maximum Julian day for leaf onset
LDay_min_H =  [14.05       12.58     12.58      12.58    12.58       11.0     11.0       11.0   12.1         11.3       12.8       NaN  ]; % Minimum Day duration for leaf onset
LDay_cr_H =   [9.45        12.3      12.3       12.3     12.3        12.3     12.3       12.3   11.8         11.3       10.8       NaN  ]; % Threshold for senescence day light [h]

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



%% Low Vegetation 
%--------------------------------------------------------------------------
% MOD_PARAM_Aurade contains aSE_L params as 5 for all crops. The model
% crashes under these paramaters. Potential values are 0, 1 or 2. Check
% this with Simone. 
% Note that wheat and barley are grasses
% Note that sunflower and rapeseed are not decidious plants, nor evergreen,
% not grasses.
%--------------------------------------------------------------------------
%categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
aSE_L   =     [2           1         1          1        1           2         2          0        NaN          2          NaN        NaN  ]; % Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen -- 2 Grass species -
Sl_L    =     [0.028       0.032     0.038      0.030    0.042       0.022     0.023      0.016    NaN          0.026      NaN        NaN  ]; % Specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_L    =     [23          16        16         20       16          23        23         40       NaN          23         NaN        NaN  ]; % Leaf Nitrogen Concentration [kgC/kgN ]
r_L     =     [0.055       0.025     0.025      0.025    0.025       0.060     0.060      0.036    NaN          0.045      NaN        NaN  ]; % Respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_L    =     [0.25        0.25      0.25       0.25     0.25        0.25      0.25       0.25     NaN          0.025      NaN        NaN  ]; % Growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
Tcold_L =     [0.0         0.0       0.0        0.0      0.0        -2.0       1.0        1.0      NaN          0          NaN        NaN  ]; % Cold Leaf Shed [°C] 

%categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high   BLdec_low  BLdec_high NoVeg]  
dd_max_L = 1./[360         50        50         50       50          20        45         365       NaN          30         NaN        NaN  ]; % Death maximum for drought [1/d]
dc_C_L   = 1./[52          52        52         52       52          52        52         182       NaN          60.83      NaN        NaN  ]; % Factor of increasing mortality for cold
drn_L    = 1./[570         365       365        365      365         550       450        900       NaN          750        NaN        NaN  ]; % turnover root  [1/d]
dsn_L    = 1./[365         365       365        365      1E10        365       365        600       NaN          365        NaN        NaN  ]; % normal transfer rate sapwood [1/d]
Mf_L     = 1./[NaN         365       365        365      1E10        50        50         50        NaN          50         NaN        NaN  ]; % Fruit maturation turnover [1/d]
Klf_L    = 1./[NaN         20        20         20       20          40        50         50        NaN          30         NaN        NaN  ]; % Dead Leaves fall turnover [1/d]

age_cr_L =    [180         110       85         90       80          180       180        730       NaN          250        NaN        NaN  ]; % [day] Critical Leaf Age
Trr_L    =    [4.0         8.5       6.5        3.5      6.5         2.0       3.5        0.4       NaN          1.0        NaN        NaN  ]; % Translocation rate [gC /m^2 d]
LtR_L    =    [0.8         1.0       1.1        1.0      1.1         0.45      0.35       0.5       NaN          0.5        NaN        NaN  ]; % Leaf to Root ratio maximum
Wm_L     =    [NaN         0.0       0.0        0.0      0.0         0.0       0.0        0.0       NaN          0.0        NaN        NaN  ] ;% wood turnover coefficient [1/d]
eps_ac_L =    [0.5         0.2       0.2        0.2      0.2         1.0       0.2        0.6       NaN          1.0        NaN        NaN  ]; % Allocation to reserve parameter [0-1]
fab_L    =    [NaN         1.0       1.0        1.0      1.0         0.0       0.0        0.75      NaN          0.0        NaN        NaN  ]; % fraction above-ground sapwood and reserve
fbe_L    = 1-fab_L;                           % fraction below-ground sapwood and reserve
ff_r_L   =    [NaN         0.5       0.5        0.2      0.5         0.1       0.1        0.1       NaN          0.1        NaN        NaN  ]; % Reference allocation to Fruit and reproduction

% Phenology 
%categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high   BLdec_low  BLdec_high NoVeg]  
Bfac_lo_L  =  [0.99        0.99      0.99        0.99    0.99        0.99       0.99       0.99      NaN         0.99       NaN        NaN  ]; % Leaf Onset Water Stress
Bfac_ls_L  =  [NaN         NaN       NaN         NaN     NaN         0.15       NaN        NaN       NaN         NaN        NaN        NaN  ] ;% 
Tlo_L      =  [6.0         6.0       5.0         12.0    5.5         1.0       -1.0        11.0      NaN         2.0        NaN        NaN  ]; % Mean Temperature for Leaf onset
Tls_L      =  [NaN         NaN       NaN         NaN     NaN         NaN        NaN        NaN       NaN         NaN        NaN        NaN  ]; % Not-used 
PAR_th_L   =  [NaN         NaN       NaN         NaN     NaN         NaN        NaN        NaN       NaN         NaN        NaN        NaN  ]; % Light Phenology Threshold 
dmg_L      =  [25          70        45          40      50          15         20         15        NaN         20         NaN        NaN  ]; % Day of Max Growth
LAI_min_L  =  [0.1         0.01      0.01        0.01    0.01        0.05       0.1        0.001     NaN         0.1        NaN        NaN  ];
mjDay_L    =  [220         200       200         367     200         250        250        210       NaN         366        NaN        NaN  ]; % Maximum Julian day for leaf onset
LDay_min_L =  [NaN         11        11          10.5    11.0        12         12.2       12.1      NaN         9.0        NaN        NaN  ]; % Minimum Day duration for leaf onset
LDay_cr_L  =  [12.5        11.2      11.2        10.5    11.4        10.2       11.2       10.2      NaN         9.0        NaN        NaN  ]; % Threshold for senescence day light [h]

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


%% INITIAL CONDITIONS
%==========================================================================
%
%==========================================================================

L_day=zeros(NNd,1);
for j=2:24:NN
    [~,~,~,~,~,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('L_day')

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

%% VEGETATION - Carbon
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
LAI_H(1,:) = LAI_Htm1(ij,:);  Rrootl_H(1,:)= Rrootl_Htm1(ij,:); 
PHE_S_H(1,:)= PHE_S_Htm1(ij,:); dflo_H(1,:)=dflo_Htm1(ij,:); AgeL_H(1,:)=AgeL_Htm1(ij,:);
e_rel_H(1,:)=e_rel_Htm1(ij,:); hc_H(1,:) =hc_Htm1(ij,:); SAI_H(1,:) = SAI_Htm1(ij,:);
B_H(1,:,:)= B_Htm1(ij,:,:);
%%%%%%%%%%%%%%%%%%
LAI_L(1,:) = LAI_Ltm1(ij,:); B_L(1,:,:)= B_Ltm1(ij,:,:); Rrootl_L(1,:)= Rrootl_Ltm1(ij,:);
PHE_S_L(1,:)=PHE_S_Ltm1(ij,:); dflo_L(1,:)=dflo_Ltm1(ij,:); AgeL_L(1,:)=AgeL_Ltm1(ij,:);
e_rel_L(1,:)=e_rel_Ltm1(ij,:); hc_L(1,:) =hc_Ltm1(ij,:); SAI_L(1,:) =SAI_Ltm1(ij,:);
%%%%%%%%%%%%%%%%%%
BLit(1)=0.0;  %% [kg DM /m2 PFT] Litter Biomass
%%%%%%%%%%%%%%%%%%
RexmyI(1,:)= [0 0 0];
%%%%%%%%%%%%%%%%%%
Nreserve_H(1,:)= Nreserve_Htm1(ij,:,:); 
Preserve_H(1,:)= Preserve_Htm1(ij,:,:); 
Kreserve_H(1,:)= Kreserve_Htm1(ij,:,:);
FNC_H(1,:)=FNC_Htm1(ij,:,:); 
NupI_H(1,:,:)= NupI_Htm1(ij,:,:); 
Nreserve_L(1,:)= Nreserve_Ltm1(ij,:,:);
Preserve_L(1,:)=Preserve_Ltm1(ij,:,:); 
Kreserve_L(1,:)=Kreserve_Ltm1(ij,:,:);
FNC_L(1,:)=FNC_Ltm1(ij);
NupI_L(1,:,:)= NupI_Ltm1(ij,:,:);
%%%%%%%%%%%%%%%%%%
TdpI_H(1,:)=TdpI_Htm1(ij,:,:); 
TdpI_L(1,:)=TdpI_Ltm1(ij,:,:);