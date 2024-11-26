function[aR,Zs,...parfor
    EvL_Zs,Inf_Zs,Bio_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Ks_Zs,Dz,...
    ms,Kbot,Krock,zatm,...
    Ccrown,Cbare,Crock,Curb,Cwat,...
    Color_Class,OM_H,OM_L,PFT_opt_H,PFT_opt_L,d_leaf_H,d_leaf_L,...
    SPAR,Phy,Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,fpr,...
    VegH_Param_Dyn,VegL_Param_Dyn,...
    Stoich_H,aSE_H,Stoich_L,aSE_L,fab_H,fbe_H,fab_L,fbe_L,...
    ZR95_H,ZR95_L,In_max_urb,In_max_rock,K_usle,...
    Urb_Par,Deb_Par,Zs_deb,... 
    Sllit,Kct,ExEM,ParEx_H,Mpar_H,ParEx_L,Mpar_L]=PARAMETERS_SOIL(ANSWER,Psan,Pcla,Porg,dbThick, md_max,zatm_surface,Afirn,Soil_th,s)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PARAMETERS SOILS AND VEGETATION  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs
%%% ANSWER;
%%% 1 Fir (everg.)
%%% 2 Larch (decid.)
%%% 3 Grass C3
%%% 4 Shrub (decid.)
%%% 5 broadleaf evergreen vegetation dec.
%%% 6 broadleaf deciduous vegetation dec.
%%% 7 Rock

%%% --> II (define high and low vegetation
%%% Fir (H) / Larch (H) / Grass C3 (L) / Shrub Winter Dec. (L) 
%%% zatm_surface - instrument height above shrub, grass meadow and rock

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TO ADAPT IN AAA_CHA_TC...:   
%%%  "CONDITION INITIAL 1 STEP" (SOIL MOISTURE, HIGH/LOW VEGETATION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ms=10; %%% Number of soil layers (has to correspond to "ms_max" in the launcher)
zatm =[18 18 zatm_surface zatm_surface 18 18 zatm_surface] ; %% Reference Height
zatm = zatm(ANSWER);
fpr = 1;
%%%%
aR =100;
%Kh=Ks*aR;
Kbot = [0.0 0.0 0.0 0.0 0.0 0.0 0.0];   %% [mm/h] Conductivity at the bedrock layer
Krock = [NaN NaN NaN NaN NaN NaN 0.15]; %% [mm/h] Conductivity of Fractured Rock
Kbot= Kbot(ANSWER);
Krock=Krock(ANSWER);

%%%%%%%%%%%%%%%%%% Soil layer depths [mm]
%%%%%%%%%%%%%%%%%%  "Zs" & "dz" have to correspond to "vi" in the launcher (->SOIL MOISTURE)
Zs= [0    10    20    50   100   150   200   300   400    700   1000]; %% ms+1
% Zs= [0    10    20    50   100   150   200   300   400   500   600   800  1000  1250  1500]; %% ms+1
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end
Zdes = 10; %%% Evaporation depth
Zinf = 10; %%% Infiltration depth
Zbio = 250;
[EvL_Zs]=Evaporation_layers(Zs,Zdes); %%% Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); %%% Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio);
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
Dz=zeros(1,ms);
for ii = 1:ms
    if ii>1
        Dz(ii)= (dz(ii)+ dz(ii-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(ii)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end

%%%%%%%%%% SOIL INPUT
Color_Class = 0;
%%%%
switch ANSWER
    case 1
        %%%% LAND COVER PARTITION  Fir - evergreen
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [1 0 0 0 0 0]>0;  
    case 2
        %%%% LAND COVER PARTITION  Larch - deciduous
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];  
        cc=length(Ccrown);%% Crown area
        II = [0 1 0 0 0 0]>0;  
    case 3
        %%%% LAND COVER PARTITION  Grass C3
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 1 0 0 0]>0;  
    case 4
        %%%% LAND COVER PARTITION  Shrub dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.1; Ccrown = [0.9];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 1 0 0]>0;  
    case 5
        %%%% LAND COVER PARTITION  broadleaf evergreen vegetation dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 0 1 0]>0;  
    case 6
        %%%% LAND COVER PARTITION  broadleaf deciduous vegetation dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 0 0 1]>0;  
    case 7
        %%%% LAND COVER PARTITION  Rock
        Cwat = 0; Curb = 0.0 ; Crock = 1.0;
        Cbare = 0.0; Ccrown = [0.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 1 0 0]>0;  
    otherwise
        disp('INDEX FOR SOIL VEGETATION PARAMETER INCONSISTENT')
        return
end
%%%%%%%%%%%%%%%%%%%
SPAR=2; %%% SOIL PARAMETER TYPE (1-VanGenuchten 2-Saxton-Rawls)
%%%%
[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters(Psan,Pcla,Porg);
%%%%%
rsd=rsd*ones(1,ms);
lan_dry=lan_dry*ones(1,ms);
lan_s =lan_s*ones(1,ms);
cv_s = cv_s*ones(1,ms);
%%%
%nVG=L+1;
%alpVG = 1/(-101.9368*Pe); %%[1/mm]%;
p=3+2/L;
m=2/(p-1); nVG= 1/(1-m);
alpVG=(((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p^2)/(147.8+8.1*p+0.092*p^2)))^-1; %%[1/mm]%;
%%%%%%
Osat=Osat*ones(1,ms);
%Ohy = Ohy*ones(1,ms) ; %% [-]
L=L*ones(1,ms);
Pe = Pe*ones(1,ms);
O33 = O33*ones(1,ms);
alpVG= alpVG*ones(1,ms); %% [1/mm]
nVG= nVG*ones(1,ms); %% [-]
Ks_Zs= Ks*ones(1,ms); %%[mm/h]

%%%%%%%%%%%%%%%% Define soil layer depth by specifying 
%%%%%%%%%%%%%%%% impermeable layer at given layer
Soil_th=Zs(end)*(Soil_th/100); %% convert from relative to absolute depth
[~, ix]= min(abs(Zs-Soil_th));
if ix==1
    ix=2;
end
Ks_Zs(ix-1:end)=1e-12; 
clear ix

%%%%%%%%%%%%%%%% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
clear Oss_Lp Owp_Lp
% force field capacity (Ofc) to equal Osat for impermeable layers
Ofc(Ks_Zs<0.2)=Osat(Ks_Zs<0.2);

Oice = 0;
%%%%%%%%%%%%%%%%%%%
%================= OTHER PARAMETER ================================
%Interception
In_max_urb=5;
In_max_rock=0.1; %% [mm]

%================== SNOW PARAMETER =================================

TminS=-1.1;%% Threshold temperature snow
TmaxS= 2.5;%% Threshold temperature snow
ros_max1=580; %600; %%% [kg/m^3]
ros_max2=300; %450; %%% [kg/m^3]
Th_Pr_sno = 8.0; %%% [mm/day] Threshold Intensity of snow to consider a New SnowFall

% =================== ICE Parameter ================================

Ice_wc_sp =0.01; %% [-] Specific Maximum water content ice
ros_Ice_thr = 500 ; %% [kg/m^3] Density Thrshold to transform snow into ice
WatFreez_Th = -8; %% [?C] Threshold for freezing lake water
dz_ice = 0.45; %% [mm / h] Water Freezing Layer progression without snow-layer

%====================== Ice Albedo ==============================

Aice = Afirn; %% [-] Ice albedo Use standard albedo

%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExEM = [0.0 0.0 0.0 0.0 0.0 0.0]; %% Fraction of ectomychorrizae per area
ExEM = ExEM(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS VEGETATION
%%% cc -- number of crown area
%%% Root Depth
CASE_ROOT=1;  %%% Type of Root Profile
%%%
ZR95_H = [800  800    0    0 1000 800]; %% [mm]
ZR95_L = [0    0  200  600 0 0 ]; %% [mm]
ZR50_H = [NaN  NaN  NaN  NaN  NaN  NaN];
ZR50_L = [NaN  NaN  NaN  NaN NaN  NaN];
ZRmax_H = [NaN  NaN  NaN  NaN NaN  NaN];
ZRmax_L = [NaN  NaN  NaN  NaN NaN  NaN];
ZR95_H =ZR95_H(II); ZR50_H =ZR50_H(II); ZRmax_H =ZRmax_H(II);
ZR95_L =ZR95_L(II); ZR50_L =ZR50_L(II); ZRmax_L =ZRmax_L(II);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% NO ROOTS IN LANDCOVER tYPE "ROCK"
if ANSWER == 7 
    ZR95_H = [0]; %% [mm]
    ZR95_L = [0]; %% [mm]
    ZR50_H = [NaN];
    ZR50_L = [NaN];
    ZRmax_H = [NaN];
    ZRmax_L = [NaN];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
% Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]
%%%%%%%%% Interception Parameter
Sp_LAI_H_In= [0.1         0.1         0.2         0.2  0.2         0.2]; %%[mm/LAI]
Sp_LAI_L_In= [0.2         0.2         0.2         0.1  0.2         0.2 ]; %%[mm/LAI]
Sp_LAI_H_In =Sp_LAI_H_In(II);
Sp_LAI_L_In =Sp_LAI_L_In(II);
%%%%%%%%%%% Leaf Dimension
%%%%%%%%%%% Leaf Dimension
d_leaf_H= [0.25         0.8           2           2  5   4 ]; %%[cm]
d_leaf_L= [2           2         0.8           3  2  2];  %% [cm]
d_leaf_H =d_leaf_H(II);
d_leaf_L =d_leaf_L(II);
%%%%%%%% Biochemical parameter
%% Veg Biochemical parameter
KnitH=[0.35           0.2           NaN           NaN  0.35 0.30 ]; %%% Canopy Nitrogen Decay
KnitL=[NaN           NaN          0.15          0.25  NaN           NaN];
%%%%%
mSl_H = [0    0  NaN  NaN 0    0];%% [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = [NaN  NaN    0    0 NaN  NaN ]; % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree
%%%%%
KnitH =KnitH(II);  mSl_H =mSl_H(II);
KnitL =KnitL(II);  mSl_L =mSl_L(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  Photosynthesis Parameter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_H=[0.081         0.081           NaN           NaN 0.081         0.081] ;% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_H=[800  700  NaN  NaN 800  1000]; %%[Pa] 
a1_H=[6    6  NaN  NaN 7 6 ];  %%% [-] WUE parameter 
go_H=[0.01          0.01           NaN           NaN 0.01          0.01 ];% [mol / s m^2] minimum Stomatal Conductance
CT_H=[3    3  NaN  NaN 3    3 ]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_H =[0.649          0.66           NaN           NaN 0.649   0.649  ];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_H =[72   94  NaN  NaN 72 76]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_H=[Inf Inf   NaN    NaN Inf Inf  ]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_H=[2           1.5           NaN           NaN 2 2 ] ; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_H=FI_H(II); Do_H=Do_H(II); a1_H=a1_H(II); go_H=go_H(II);
CT_H=CT_H(II); DSE_H=DSE_H(II); Ha_H=Ha_H(II); gmes_H=gmes_H(II);
rjv_H=rjv_H(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FI_L=[NaN           NaN         0.081         0.081 NaN           NaN  ];% Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L=[NaN   NaN  1000  1000 NaN           NaN  ]; %%[Pa] 
a1_L=[NaN  NaN    8    7 NaN           NaN  ];  %%% [-] WUE parameter 
go_L=[NaN           NaN          0.01          0.01 NaN           NaN  ];% [mol / s m^2] minimum Stomatal Conductance
CT_L=[NaN  NaN    3    3 NaN NaN]; %%--> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
DSE_L =[NaN           NaN         0.649         0.649 NaN           NaN  ];  %% [kJ/mol] Activation Energy - Plant Dependent
Ha_L =[NaN  NaN   62   72 NaN           NaN  ]; %% [kJ / mol K]  entropy factor - Plant Dependent
gmes_L=[NaN    NaN    Inf Inf NaN           NaN  ]; %% [mol CO2 / s m^2 ];  mesophyll conductance
rjv_L=[NaN           NaN           1.9           2.2 NaN           NaN  ]; %%% Ratio Jmax - Vmax  [umol electrons / umolCO2 ]
%%%%%%%
FI_L=FI_L(II); Do_L=Do_L(II); a1_L=a1_L(II); go_L=go_L(II);
CT_L=CT_L(II); DSE_L=DSE_L(II); Ha_L=Ha_L(II); gmes_L=gmes_L(II);
rjv_L=rjv_L(II);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydraulic Parameters
Psi_sto_00_H = [-0.8          -0.8           NaN           NaN  -1.0 -0.8]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_H = [-2.5          -2.5           NaN           NaN -2.8 -2.5] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_H = [-1   -1  NaN  NaN -1.2 -1.0]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_H = [-3.2          -3.2           NaN           NaN -4.0 -3.0] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_H = [10   10  NaN  NaN 10   10 ] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_H  = [1200  1200   NaN   NaN 1200  1200  ];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_H = [15   15  NaN  NaN 15   15 ] ; %% [cm^2 stem /m^2 PFT]
Kx_max_H = [80000  80000    NaN    NaN 80000  80000 ];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_H = [-5   -5  NaN  NaN -6   -4.5]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_H= [150  150  NaN  NaN 150  150 ]; %%% [kg / m^3 sapwood MPa]
%%------------------------
%%% Stomata
Psi_sto_00_L= [NaN           NaN          -0.5            -1 NaN           NaN  ]; %% [MPa]  Water Potential at 2% loss conductivity
Psi_sto_50_L = [NaN           NaN          -2.8            -3 NaN           NaN  ] ;%% [MPa]  Water Potential at 50% loss conductivity
%%% Leaf
PsiL00_L = [NaN           NaN            -1          -2.5 NaN           NaN  ]; %% [MPa]  Water Potential at 2% loss conductivity
PsiL50_L = [NaN           NaN          -3.5          -4.5 NaN           NaN  ] ;%%[MPa]  Water Potential at 50% loss conductivity
Kleaf_max_L = [NaN  NaN    5    5 NaN           NaN  ] ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
Cl_L  = [NaN   NaN  1200  1200 NaN  NaN];  %%%  Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]%
%%% Xylem
Axyl_L = [NaN  NaN    0    0 NaN  NaN] ; %% [cm^2 stem /m^2 PFT]
Kx_max_L = [NaN    NaN  80000  80000 NaN  NaN];  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
PsiX50_L = [NaN           NaN          -4.5            -9 NaN  NaN]; %%[MPa]  Water Potential at 50% loss conductivity
Cx_L= [NaN  NaN  150  150 NaN  NaN]; %%% [kg / m^3 sapwood MPa]
%%%
Psi_sto_50_H =Psi_sto_50_H(II);  Psi_sto_00_H =Psi_sto_00_H(II);
PsiL00_H = PsiL00_H(II); PsiL50_H=PsiL50_H(II);  Kleaf_max_H=Kleaf_max_H(II);
Cl_H=Cl_H(II); Axyl_H=Axyl_H(II); Kx_max_H=Kx_max_H(II); PsiX50_H=PsiX50_H(II); Cx_H=Cx_H(II);
Psi_sto_50_L =Psi_sto_50_L(II);  Psi_sto_00_L =Psi_sto_00_L(II);
PsiL00_L = PsiL00_L(II); PsiL50_L=PsiL50_L(II);  Kleaf_max_L=Kleaf_max_L(II);
Cl_L=Cl_L(II); Axyl_L=Axyl_L(II); Kx_max_L=Kx_max_L(II); PsiX50_L=PsiX50_L(II); Cx_L=Cx_L(II);

%%%%%%%%%%%%%%%% Root Parameters
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

%% Growth Parameters
PsiG50_H= [-0.5          -0.8           NaN           NaN NaN  NaN];  %%[MPa]
PsiG99_H= [-2.5          -2.5           NaN           NaN NaN  NaN];  %%[MPa]
gcoef_H = [3.5           3.5           NaN           NaN NaN  NaN]; % [gC/m2 day]
%%------
PsiG50_L= [NaN           NaN          -2.8            -3 NaN  NaN];
PsiG99_L= [NaN           NaN            -4          -4.5 NaN  NaN];
gcoef_L = [NaN           NaN           3.5           3.5 NaN  NaN]; % [gC/m2 day]
%%%%%
PsiG50_H=PsiG50_H(II); PsiG99_H=PsiG99_H(II); gcoef_H=gcoef_H(II);
PsiG50_L=PsiG50_L(II); PsiG99_L=PsiG99_L(II); gcoef_L=gcoef_L(II);

OPT_PROP_H =[2 3 0 0 5 7];   % =PFT_Class for "Veg_Optical_Parameter"-function
OPT_PROP_L =[0 0 13 2 0 0];
OPT_PROP_H=OPT_PROP_H(II);
OPT_PROP_L=OPT_PROP_L(II);
for i=1:cc
    %%%%%%%% Vegetation Optical Parameter
    [PFT_opt_H(i)]=Veg_Optical_Parameter(OPT_PROP_H(i));
    [PFT_opt_L(i)]=Veg_Optical_Parameter(OPT_PROP_L(i));
end

OM_H=1;
OM_L=1;
Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VEGETATION PART %%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HIGH VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%
aSE_H=  [0    1  NaN  NaN 0 1 ]; %% Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_H =  [0.010         0.025           NaN           NaN 0.016 0.020] ; %  [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_H =  [42   26  NaN  NaN 40 28]; %[kgC/kgN ] Leaf Nitrogen Concentration
r_H =   [0.058         0.055           NaN           NaN 0.045         0.035 ];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_H=   [0.25          0.25           NaN           NaN 0.25          0.25  ];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_H= 1./[150  200  NaN  NaN 200 100];%%%  [1/d]  death maximum for drought
dc_C_H =  1./[5      10          NaN           NaN 365 10]; %%[Factor of increasing mortality for cold]
Tcold_H = [-40   -3  NaN  NaN  3 3.5]; %% [°C] Cold Leaf Shed
drn_H=  1./[900  1100   NaN   NaN 550 800]; %% turnover root  [1/d]
dsn_H= 1./[1100   750   NaN   NaN 800 700]; % normal transfer rate sapwood [1/d]
age_cr_H= [950  180  NaN  NaN 365 180]; %% [day] Critical Leaf Age
Trr_H = [0.25             3           NaN           NaN 0.5 3.5]; %% Translocation rate [gC /m^2 d]
LtR_H = [0.8           0.8           NaN           NaN 1.0 0.9]; %%% Leaf to Root ratio maximum
Mf_H= 1./[80   50  NaN  NaN 80 80]; %% fruit maturation turnover [1/d]
Wm_H= [0    0  NaN  NaN 0 0] ; % wood turnover coefficient [1/d]
eps_ac_H = [0.25             1           NaN           NaN 0.5 1.0]; %% Allocation to reserve parameter [0-1]
Klf_H = 1./[40   30  NaN  NaN 30 28]; %% Dead Leaves fall turnover [1/d]
fab_H = [0.74           0.8           NaN           NaN 0.74  0.74 ]; %% fraction above-ground sapwood and reserve
fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
ff_r_H= [0.1           0.1           NaN           NaN 0.1           0.1  ]; %% Reference allocation to Fruit and reproduction
%%%% Phenology 
Bfac_lo_H= [0.99          0.99           NaN           NaN 0.95         0.95   ]; %% Leaf Onset Water Stress
Bfac_ls_H= [NaN  NaN  NaN  NaN NaN  NaN] ;  % Not-used 
Tlo_H = [4.5           3.5           NaN           NaN 5.5 2.8]; %% Mean Temperature for Leaf onset
Tls_H = [NaN  NaN  NaN  NaN NaN  NaN] ; %%% Not-used 
PAR_th_H= [NaN  NaN  NaN  NaN NaN  NaN]; %% Light Phenology Threshold 
dmg_H= [30   30  NaN  NaN 45 30]; %%%  Day of Max Growth
LAI_min_H = [0.001          0.01           NaN           NaN 0.001          0.01 ];
mjDay_H = [220  250  NaN  NaN 250 250]; %% Maximum Julian day for leaf onset
LDay_min_H =[12.8           12.7           NaN           NaN 12.1 11.7]; %% Minimum Day duration for leaf onset
LDay_cr_H = [11.8          11.6           NaN           NaN 11.8 12.0]; %%%  Threshold for senescence day light [h]
%%%
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LOW VEGETATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
aSE_L=  [NaN  NaN    2    1     NaN  NaN  ]; % Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_L =  [NaN  NaN   0.016 0.018 NaN  NaN  ]; % [m^2 gC] specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_L =  [NaN  NaN   23     35   NaN  NaN  ]; % [kgC/kgN ] Leaf Nitrogen Concentration
r_L =  [NaN           NaN         0.025         0.025 NaN  NaN  ];  %% respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_L=  [NaN           NaN          0.25          0.25 NaN  NaN  ];%  growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
dd_max_L= 1./[NaN  NaN   45  100 NaN  NaN  ];%%%  [1/d]  death maximum for drought
dc_C_L =  1./[NaN           NaN      32     15 NaN  NaN  ]; %%[Factor of increasing mortality for cold]
Tcold_L = [NaN  NaN    3    2 NaN  NaN  ]; %% [°C] Cold Leaf Shed
drn_L=  1./[NaN   NaN   950  1600 NaN  NaN  ]; %% turnover root  [1/d]
dsn_L= 1./[NaN   NaN   365  1800 NaN  NaN  ]; % normal transfer rate sapwood [1/d]
age_cr_L= [NaN  NaN  180  180 NaN  NaN  ]; %% [day] Critical Leaf Age
Trr_L = [NaN           NaN             2           0.6 NaN  NaN  ]; %% Translocation rate [gC /m^2 d]
LtR_L = [NaN           NaN           0.7             1 NaN  NaN  ]; %%% Leaf to Root ratio maximum
Mf_L= 1./[NaN  NaN   50   50 NaN  NaN  ]; %% fruit maturation turnover [1/d]
Wm_L= [NaN  NaN    0    0 NaN  NaN  ] ; % wood turnover coefficient [1/d]
eps_ac_L = [NaN  NaN    1    1 NaN  NaN  ]; %% Allocation to reserve parameter [0-1]
Klf_L = 1./[NaN  NaN   20   40 NaN  NaN  ]; %% Dead Leaves fall turnover [1/d]
fab_L = [NaN           NaN             0          0.75 NaN  NaN  ]; %% fraction above-ground sapwood and reserve
fbe_L = 1-fab_L; %% fraction below-ground sapwood and reserve
ff_r_L= [NaN           NaN           0.1           0.1 NaN  NaN  ]; %% Reference allocation to Fruit and reproduction
%%%% Phenology 
Bfac_lo_L= [NaN           NaN          0.99          0.99 NaN  NaN  ]; %% Leaf Onset Water Stress
Bfac_ls_L= [NaN  NaN  0.15 NaN NaN  NaN  ] ; % 
Tlo_L = [NaN  NaN   2.5    2.5 NaN  NaN  ]; %% Mean Temperature for Leaf onset
Tls_L = [NaN  NaN  NaN  NaN NaN  NaN  ]; %% Not-used 
PAR_th_L= [NaN  NaN  NaN  NaN NaN  NaN  ]; %% Light Phenology Threshold 
dmg_L= [NaN  NaN   30  25 NaN  NaN  ]; %%%  Day of Max Growth
LAI_min_L = [NaN           NaN          0.05         0.001 NaN  NaN  ];
mjDay_L = [NaN  NaN  250  180 NaN  NaN  ]; %% Maximum Julian day for leaf onset
LDay_min_L =[NaN  NaN  12  12.2 NaN  NaN  ]; %% Minimum Day duration for leaf onset
LDay_cr_L = [NaN  NaN   12   12 NaN  NaN  ]; %%%  Threshold for senescence day light [h]
%%%
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
if ANSWER == 3
    Mpar_L(1).jDay_cut=[180:243];
    Mpar_L(1).LAI_cut=[-0.1]; %% LAI of grass after cut
end
%%%%
Vmax_H = [45  64   0   0 32 48]; % [umol CO2 /m2 s] - Maximum Rubisco Capacity 
Vmax_L = [0   0  50  46 0 0]; % [umol CO2 /m2 s]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vmax_H =Vmax_H(II); Vmax_L =Vmax_L(II);
%%%%%
%%%%%%%%%%%%%%%%%
Restating_parameters;
if dbThick>5
    nst=md_max;
    k=(dbThick/5)^(1/(nst-1));
    Zs_deb = [0 5*k.^(1:nst-1)]; %% [mm]
    clear nst k
else
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Restate debris parameters here, otherwise missing as output??
albs = [0.153 0.115 0.13 0.13 0.13 0.13];
lans = [1.65 0.985 1.45 0.94 0.94 0.94];
zoms = [0.38 0.081 0.15 0.016 0.016 0.016];

Deb_Par.alb= albs(s);
Deb_Par.e_sur =  0.94;
Deb_Par.lan = lans(s);
Deb_Par.rho = 1496;  % [kg/m^3]
Deb_Par.cs = 948;   % [J/kg K]
Deb_Par.zom = zoms(s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%