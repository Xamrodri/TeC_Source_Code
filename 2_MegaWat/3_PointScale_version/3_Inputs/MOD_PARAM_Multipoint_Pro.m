%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

PARAMETER FILE TEMPLATE FOR T&C 

This code uses TT_par defined in the Starter_MultiPoint.m
TT_par is the table with the model parameters. It is started in
Starter_MultiPoint.m and sent to the function in run_Point_Por.m

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
albs = table2array(TT_par(strcmp(TT_par.Parameters,'albs'),II));
lans = table2array(TT_par(strcmp(TT_par.Parameters,'lans'),II));
zoms = table2array(TT_par(strcmp(TT_par.Parameters,'zoms'),II));


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

ExEM = table2array(TT_par(strcmp(TT_par.Parameters,'ExEM'),II));

CASE_ROOT= 1;  % Type of Root Profile

% Selection based on II
ZR95_H = table2array(TT_par(strcmp(TT_par.Parameters,'ZR95_H'),II));
ZR95_L = table2array(TT_par(strcmp(TT_par.Parameters,'ZR95_L'),II));
ZR50_H = table2array(TT_par(strcmp(TT_par.Parameters,'ZR50_H'),II));
ZR50_L = table2array(TT_par(strcmp(TT_par.Parameters,'ZR50_L'),II));
ZRmax_H = table2array(TT_par(strcmp(TT_par.Parameters,'ZRmax_H'),II));
ZRmax_L = table2array(TT_par(strcmp(TT_par.Parameters,'ZRmax_L'),II));


%% INTERCEPTION PARAMETERS 

In_max_urb= 5;    % Maximum interception capacity in urban [mm]
In_max_rock= 0.1; % Maximum interception capacity in rocks [mm]
Kct=0.75;         % Foliage cover decay factor for throughfall [-]

%% Interception Parameter
gcI=3.7;       % Interception parameter [1/mm]
KcI=0.06;      % Interception drainage rate coefficient [mm] - Mahfouf and Jacquemin 1989
Sp_SN_In= 5.9; % Specific interception of rainfall for unit leaf area. Average of high vegetation [mm/LAI]

Sp_LAI_H_In = table2array(TT_par(strcmp(TT_par.Parameters,'Sp_LAI_H_In'),II));
Sp_LAI_L_In = table2array(TT_par(strcmp(TT_par.Parameters,'Sp_LAI_L_In'),II));

%% Leaf Dimension
d_leaf_H = table2array(TT_par(strcmp(TT_par.Parameters,'d_leaf_H'),II));
d_leaf_L = table2array(TT_par(strcmp(TT_par.Parameters,'d_leaf_L'),II));

%% Veg Biochemical parameter
KnitH = table2array(TT_par(strcmp(TT_par.Parameters,'KnitH'),II));
KnitL = table2array(TT_par(strcmp(TT_par.Parameters,'KnitL'),II));

mSl_H = table2array(TT_par(strcmp(TT_par.Parameters,'mSl_H'),II));
mSl_L = table2array(TT_par(strcmp(TT_par.Parameters,'mSl_L'),II));

%%  Photosynthesis Parameter
% High vegetation
FI_H = table2array(TT_par(strcmp(TT_par.Parameters,'FI_H'),II));
Do_H = table2array(TT_par(strcmp(TT_par.Parameters,'Do_H'),II));
a1_H = table2array(TT_par(strcmp(TT_par.Parameters,'a1_H'),II));
go_H = table2array(TT_par(strcmp(TT_par.Parameters,'go_H'),II));
CT_H = table2array(TT_par(strcmp(TT_par.Parameters,'CT_H'),II));
DSE_H = table2array(TT_par(strcmp(TT_par.Parameters,'DSE_H'),II));
Ha_H = table2array(TT_par(strcmp(TT_par.Parameters,'Ha_H'),II));
gmes_H = table2array(TT_par(strcmp(TT_par.Parameters,'gmes_H'),II));
rjv_H = table2array(TT_par(strcmp(TT_par.Parameters,'rjv_H'),II));

% Low vegetation
FI_L = table2array(TT_par(strcmp(TT_par.Parameters,'FI_L'),II));
Do_L = table2array(TT_par(strcmp(TT_par.Parameters,'Do_L'),II));
a1_L = table2array(TT_par(strcmp(TT_par.Parameters,'a1_L'),II));
go_L = table2array(TT_par(strcmp(TT_par.Parameters,'go_L'),II));
CT_L = table2array(TT_par(strcmp(TT_par.Parameters,'CT_L'),II));
DSE_L = table2array(TT_par(strcmp(TT_par.Parameters,'DSE_L'),II));
Ha_L = table2array(TT_par(strcmp(TT_par.Parameters,'Ha_L'),II));
gmes_L = table2array(TT_par(strcmp(TT_par.Parameters,'gmes_L'),II));
rjv_L = table2array(TT_par(strcmp(TT_par.Parameters,'rjv_L'),II));

Vmax_H = table2array(TT_par(strcmp(TT_par.Parameters,'Vmax_H'),II));
Vmax_L = table2array(TT_par(strcmp(TT_par.Parameters,'Vmax_L'),II));

%% Hydraulic Parameters
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
Psi_sto_00_H = table2array(TT_par(strcmp(TT_par.Parameters,'Psi_sto_00_H'),II));
Psi_sto_50_H = table2array(TT_par(strcmp(TT_par.Parameters,'Psi_sto_50_H'),II));

% Leaf
PsiL00_H = table2array(TT_par(strcmp(TT_par.Parameters,'PsiL00_H'),II));
PsiL50_H = table2array(TT_par(strcmp(TT_par.Parameters,'PsiL50_H'),II));
Kleaf_max_H = table2array(TT_par(strcmp(TT_par.Parameters,'Kleaf_max_H'),II));
Cl_H = table2array(TT_par(strcmp(TT_par.Parameters,'Cl_H'),II));

% Xylem
Axyl_H= table2array(TT_par(strcmp(TT_par.Parameters,'Axyl_H'),II));
Kx_max_H = table2array(TT_par(strcmp(TT_par.Parameters,'Kx_max_H'),II));
PsiX50_H = table2array(TT_par(strcmp(TT_par.Parameters,'PsiX50_H'),II));
Cx_H = table2array(TT_par(strcmp(TT_par.Parameters,'Cx_H'),II));

% Stomata
Psi_sto_00_L = table2array(TT_par(strcmp(TT_par.Parameters,'Psi_sto_00_L'),II));
Psi_sto_50_L = table2array(TT_par(strcmp(TT_par.Parameters,'Psi_sto_50_L'),II));

% Leaf
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A    grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
PsiL00_L = table2array(TT_par(strcmp(TT_par.Parameters,'PsiL00_L'),II));
PsiL50_L = table2array(TT_par(strcmp(TT_par.Parameters,'PsiL50_L'),II));
Kleaf_max_L = table2array(TT_par(strcmp(TT_par.Parameters,'Kleaf_max_L'),II));
Cl_L = table2array(TT_par(strcmp(TT_par.Parameters,'Cl_L'),II));

% Xylem
Axyl_L = table2array(TT_par(strcmp(TT_par.Parameters,'Axyl_L'),II));
Kx_max_L = table2array(TT_par(strcmp(TT_par.Parameters,'Kx_max_L'),II));
PsiX50_L = table2array(TT_par(strcmp(TT_par.Parameters,'PsiX50_L'),II));
Cx_L = table2array(TT_par(strcmp(TT_par.Parameters,'Cx_L'),II));


%% Root Parameters (Function)
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

%% Growth Parameters
% High vegetation
%categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R    grass_A   grass_B    shrub  BLever_high  BLdec_low  BLdec_high NoVeg]  
PsiG50_H = table2array(TT_par(strcmp(TT_par.Parameters,'PsiG50_H'),II));
PsiG99_H = table2array(TT_par(strcmp(TT_par.Parameters,'PsiG99_H'),II));
gcoef_H = table2array(TT_par(strcmp(TT_par.Parameters,'gcoef_H'),II));

% Low vegetation
PsiG50_L = table2array(TT_par(strcmp(TT_par.Parameters,'PsiG50_L'),II));
PsiG99_L = table2array(TT_par(strcmp(TT_par.Parameters,'PsiG99_L'),II));
gcoef_L = table2array(TT_par(strcmp(TT_par.Parameters,'gcoef_L'),II));


%% Vegetation Optical Parameter
OPT_PROP_H = table2array(TT_par(strcmp(TT_par.Parameters,'OPT_PROP_H'),II));
OPT_PROP_L = table2array(TT_par(strcmp(TT_par.Parameters,'OPT_PROP_L'),II));

OM_H = table2array(TT_par(strcmp(TT_par.Parameters,'OM_H'),II));
OM_L = table2array(TT_par(strcmp(TT_par.Parameters,'OM_L'),II));


%% Specific leaf area of litter
Sllit = 2 ; % Litter Specific Leaf area [m2 Litter / kg DM]

%% High Vegetation
%BLever: broadleaf evergreen vegetation dec;
%BLdec: broadleaf deciduous vegetation dec

aSE_H = table2array(TT_par(strcmp(TT_par.Parameters,'aSE_H'),II));
Sl_H = table2array(TT_par(strcmp(TT_par.Parameters,'Sl_H'),II));
Nl_H = table2array(TT_par(strcmp(TT_par.Parameters,'Nl_H'),II));
r_H = table2array(TT_par(strcmp(TT_par.Parameters,'r_H'),II));
gR_H = table2array(TT_par(strcmp(TT_par.Parameters,'gR_H'),II));
Tcold_H = table2array(TT_par(strcmp(TT_par.Parameters,'Tcold_H'),II));
age_cr_H = table2array(TT_par(strcmp(TT_par.Parameters,'age_cr_H'),II));
Trr_H = table2array(TT_par(strcmp(TT_par.Parameters,'Trr_H'),II));
LtR_H = table2array(TT_par(strcmp(TT_par.Parameters,'LtR_H'),II));
eps_ac_H = table2array(TT_par(strcmp(TT_par.Parameters,'eps_ac_H'),II));
fab_H = table2array(TT_par(strcmp(TT_par.Parameters,'fab_H'),II));
ff_r_H = table2array(TT_par(strcmp(TT_par.Parameters,'ff_r_H'),II));
Wm_H = table2array(TT_par(strcmp(TT_par.Parameters,'Wm_H'),II));


dd_max_H = table2array(TT_par(strcmp(TT_par.Parameters,'dd_max_H'),II));
dc_C_H = table2array(TT_par(strcmp(TT_par.Parameters,'dc_C_H'),II));
drn_H = table2array(TT_par(strcmp(TT_par.Parameters,'drn_H'),II));
dsn_H = table2array(TT_par(strcmp(TT_par.Parameters,'dsn_H'),II));
Mf_H = table2array(TT_par(strcmp(TT_par.Parameters,'Mf_H'),II));
Klf_H = table2array(TT_par(strcmp(TT_par.Parameters,'Klf_H'),II));

%% check Mf_H in shrub for 1/0

fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
% [Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
% [Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
% [Stoich_H(3)]=Veg_Stoichiometric_Parameter(Nl_H(3));
% [Stoich_H(4)]=Veg_Stoichiometric_Parameter(Nl_H(4));
% [Stoich_H(5)]=Veg_Stoichiometric_Parameter(Nl_H(5));
% [Stoich_H(6)]=Veg_Stoichiometric_Parameter(Nl_H(6));

%% Phenology 
Bfac_lo_H = table2array(TT_par(strcmp(TT_par.Parameters,'Bfac_lo_H'),II));
Bfac_ls_H = table2array(TT_par(strcmp(TT_par.Parameters,'Bfac_ls_H'),II));
Tlo_H = table2array(TT_par(strcmp(TT_par.Parameters,'Tlo_H'),II));
Tls_H = table2array(TT_par(strcmp(TT_par.Parameters,'Tls_H'),II));
PAR_th_H = table2array(TT_par(strcmp(TT_par.Parameters,'PAR_th_H'),II));
dmg_H = table2array(TT_par(strcmp(TT_par.Parameters,'dmg_H'),II));
LAI_min_H = table2array(TT_par(strcmp(TT_par.Parameters,'LAI_min_H'),II));
mjDay_H  = table2array(TT_par(strcmp(TT_par.Parameters,'mjDay_H'),II));
LDay_min_H  = table2array(TT_par(strcmp(TT_par.Parameters,'LDay_min_H'),II));
LDay_cr_H = table2array(TT_par(strcmp(TT_par.Parameters,'LDay_cr_H'),II));


%% Low Vegetation 
%--------------------------------------------------------------------------
% MOD_PARAM_Aurade contains aSE_L params as 5 for all crops. The model
% crashes under these paramaters. Potential values are 0, 1 or 2. Check
% this with Simone. 
% Note that wheat and barley are grasses
% Note that sunflower and rapeseed are not decidious plants, nor evergreen,
% not grasses.
%--------------------------------------------------------------------------

aSE_L = table2array(TT_par(strcmp(TT_par.Parameters,'aSE_L'),II));
Sl_L = table2array(TT_par(strcmp(TT_par.Parameters,'Sl_L'),II));
Nl_L = table2array(TT_par(strcmp(TT_par.Parameters,'Nl_L'),II));
r_L = table2array(TT_par(strcmp(TT_par.Parameters,'r_L'),II));
gR_L = table2array(TT_par(strcmp(TT_par.Parameters,'gR_L'),II));
Tcold_L = table2array(TT_par(strcmp(TT_par.Parameters,'Tcold_L'),II));

dd_max_L = table2array(TT_par(strcmp(TT_par.Parameters,'dd_max_L'),II));
dc_C_L = table2array(TT_par(strcmp(TT_par.Parameters,'dc_C_L'),II));
drn_L = table2array(TT_par(strcmp(TT_par.Parameters,'drn_L'),II));
dsn_L = table2array(TT_par(strcmp(TT_par.Parameters,'dsn_L'),II));
Mf_L = table2array(TT_par(strcmp(TT_par.Parameters,'Mf_L'),II));
Klf_L = table2array(TT_par(strcmp(TT_par.Parameters,'Klf_L'),II));

age_cr_L = table2array(TT_par(strcmp(TT_par.Parameters,'age_cr_L'),II));
Trr_L = table2array(TT_par(strcmp(TT_par.Parameters,'Trr_L'),II));
LtR_L = table2array(TT_par(strcmp(TT_par.Parameters,'LtR_L'),II));
Wm_L = table2array(TT_par(strcmp(TT_par.Parameters,'Wm_L'),II));
eps_ac_L = table2array(TT_par(strcmp(TT_par.Parameters,'eps_ac_L'),II));
fab_L = table2array(TT_par(strcmp(TT_par.Parameters,'fab_L'),II));
ff_r_L = table2array(TT_par(strcmp(TT_par.Parameters,'ff_r_L'),II));
fbe_L    = 1-fab_L;  % fraction below-ground sapwood and reserve

% Phenology 
%categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high   BLdec_low  BLdec_high NoVeg]  
Bfac_lo_L = table2array(TT_par(strcmp(TT_par.Parameters,'Bfac_lo_L'),II));
Bfac_ls_L = table2array(TT_par(strcmp(TT_par.Parameters,'Bfac_ls_L'),II));
Tlo_L = table2array(TT_par(strcmp(TT_par.Parameters,'Tlo_L'),II));
Tls_L = table2array(TT_par(strcmp(TT_par.Parameters,'Tls_L'),II));
PAR_th_L = table2array(TT_par(strcmp(TT_par.Parameters,'PAR_th_L'),II));
dmg_L = table2array(TT_par(strcmp(TT_par.Parameters,'dmg_L'),II));
LAI_min_L = table2array(TT_par(strcmp(TT_par.Parameters,'LAI_min_L'),II));
mjDay_L = table2array(TT_par(strcmp(TT_par.Parameters,'mjDay_L'),II));
LDay_min_L = table2array(TT_par(strcmp(TT_par.Parameters,'LDay_min_L'),II));
LDay_cr_L = table2array(TT_par(strcmp(TT_par.Parameters,'LDay_cr_L'),II));


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