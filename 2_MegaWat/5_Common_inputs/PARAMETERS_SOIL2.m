%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% PARAMETERS SOILS AND VEGETATION  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% --> II (define high and low vegetation
%%% Fir (H) / Larch (H) / Grass C3 (L) / Shrub Winter Dec. (L) 
%%% zatm_surface - instrument height above shrub, grass meadow and rock

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function[aR,              Zs,         EvL_Zs,      Inf_Zs,     Bio_Zs, ...
         Zinf,            RfH_Zs,     RfL_Zs,      dz,         Ks_Zs, ...
         Dz,              ms,         Kbot,        Krock,      zatm,...
         Ccrown,          Cbare,      Crock,       Curb,       Cwat,...
         Color_Class,     OM_H,       OM_L,        PFT_opt_H,  PFT_opt_L, ...
         d_leaf_H,        d_leaf_L,   SPAR,Phy,    Soil_Param, Interc_Param, ...
         SnowIce_Param,   VegH_Param, VegL_Param,  fpr,        VegH_Param_Dyn, ...
         VegL_Param_Dyn,  Stoich_H,   aSE_H,       Stoich_L,   aSE_L, ...
         fab_H,           fbe_H,      fab_L,       fbe_L,      ZR95_H, ...
         ZR95_L,          In_max_urb, In_max_rock, K_usle,     Urb_Par, ...
         Deb_Par,         Zs_deb,     Sllit,       Kct,        ExEM, ...
         ParEx_H,         Mpar_H,     ParEx_L,     Mpar_L] ...
                     =PARAMETERS_SOIL( ...
         cell_class,      Psan,       Pcla,        Porg,       dbThick, ...
         md_max,          Afirn,      Soil_th,     POI,        TT_par)

fpr = 1;
aR =100;
%Kh=Ks*aR;

%% SOIL PARAMETERS
%==========================================================================
ms=10; % Number of soil layers (has to correspond to "ms_max" in the launcher)

%         Depth1 Depth2 Depth3  Depth4 Depth5 Depth6 Depth7 Depth8 Depth9  Depth10
Kbot  =  [0.0    0.0    0.0     0.0    0.0    0.0    5      0.0    0.0     0.0  ]; % Conductivity at the bedrock layer [mm/h] 
Krock =  [NaN    NaN    NaN     NaN    NaN    NaN    0.15   NaN    NaN     NaN  ]; % Conductivity of Fractured Rock [mm/h] 

Kbot = Kbot(cell_class); % Conductivity of the bedrock [mm/h] 
Krock =Krock(cell_class); % Hydraulic conductivity fractured rock [mm/h]

%% Soil layer depths [mm]
%--------------------------------------------------------------------------
% "Zs" & "dz" have to correspond to "vi" in the launcher (->SOIL MOISTURE)
% Soil layers
%--------------------------------------------------------------------------

%    Depth1  Depth2  Depth3   Depth4  Depth5  Depth6   Depth7   Depth8   Depth9  Depth10
Zs= [0       10      20       50      100     150      200      300      1000    1500     2500]; % Depth of top of the soil layer [mm],  ms+1
Zdes = 10; % Depth of evaporation layer [mm]
Zinf=  10; % Depth of infiltration layer (=first layer) [mm]
Zbio = 250; % Depth of the active Biogeochemistry zone [mm]

% Warning for soil
%--------------------------------------------------------------------------
if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    return
end

[EvL_Zs]=Evaporation_layers(Zs,Zdes); % Evaporation Layer fraction
[Inf_Zs]=Evaporation_layers(Zs,Zinf); % Infiltration Depth Layer fraction
[Bio_Zs]=Evaporation_layers(Zs,Zbio);
dz= diff(Zs); % Thickness of the Layers [mm]
Dz=zeros(1,ms);

for ii = 1:ms
    if ii>1
        Dz(ii)= (dz(ii)+ dz(ii-1))/2; %%% Delta Depth Between Middle Layer  [mm]
    else
        Dz(ii)=dz(1)/2; %%% Delta Depth Between First Middle Layer and soil surface [mm]
    end
end

%% VEGETATION
Color_Class = 0;

Cwat    = POI.Cwat{cell_class}; 
Curb    = POI.Curb{cell_class}; 
Crock   = POI.Crock{cell_class}; 
Cbare   = POI.Cbare{cell_class}; 
Ccrown  = POI.Ccrowns{cell_class};
cc      = length(POI.Ccrowns{cell_class});
II      = POI.Veg_type{cell_class}; 

% Elevation 
zatm = table2array(TT_par(strcmp(TT_par.Parameters,'zatm_surface'),II)); %% Reference Height
zatm = max(zatm);
%% SOIL
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

%%
ExEM = table2array(TT_par(strcmp(TT_par.Parameters,'ExEM'),II));

%% PARAMETERS VEGETATION
%==========================================================================
ExEM = table2array(TT_par(strcmp(TT_par.Parameters,'ExEM'),II));

CASE_ROOT= 1;  % Type of Root Profile

% Selection based on II
ZR95_H = table2array(TT_par(strcmp(TT_par.Parameters,'ZR95_H'),II));
ZR95_L = table2array(TT_par(strcmp(TT_par.Parameters,'ZR95_L'),II));
ZR50_H = table2array(TT_par(strcmp(TT_par.Parameters,'ZR50_H'),II));
ZR50_L = table2array(TT_par(strcmp(TT_par.Parameters,'ZR50_L'),II));
ZRmax_H = table2array(TT_par(strcmp(TT_par.Parameters,'ZRmax_H'),II));
ZRmax_L = table2array(TT_par(strcmp(TT_par.Parameters,'ZRmax_L'),II));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kct=0.75; %%% Factor Vegetation Cover --- for throughfall
% Interception Parameter
gcI=3.7; %%% [1/mm]
KcI=0.06; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Interception Parameter
Sp_SN_In= 5.9; %% [mm/LAI]

% Interception Parameter
%--------------------------------------------------------------------------
Sp_LAI_H_In = table2array(TT_par(strcmp(TT_par.Parameters,'Sp_LAI_H_In'),II));
Sp_LAI_L_In = table2array(TT_par(strcmp(TT_par.Parameters,'Sp_LAI_L_In'),II));

% Leaf Dimension
%--------------------------------------------------------------------------
d_leaf_H = table2array(TT_par(strcmp(TT_par.Parameters,'d_leaf_H'),II));
d_leaf_L = table2array(TT_par(strcmp(TT_par.Parameters,'d_leaf_L'),II));

%% Biochemical parameterS
%==========================================================================

% Veg Biochemical parameter
%--------------------------------------------------------------------------
KnitH = table2array(TT_par(strcmp(TT_par.Parameters,'KnitH'),II));
KnitL = table2array(TT_par(strcmp(TT_par.Parameters,'KnitL'),II));

mSl_H = table2array(TT_par(strcmp(TT_par.Parameters,'mSl_H'),II)); % [m2 PFT /gC]  Linear increase in Sla with LAI
mSl_L = table2array(TT_par(strcmp(TT_par.Parameters,'mSl_L'),II)); % [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree


%  Photosynthesis Parameter
%--------------------------------------------------------------------------
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


% Hydraulic Parameters
%--------------------------------------------------------------------------
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

% Root Parameters (Function)
[RfH_Zs,RfL_Zs]=Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L);

% Growth Parameters
%--------------------------------------------------------------------------
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


for i=1:cc
    %%%%%%%% Vegetation Optical Parameter
    [PFT_opt_H(i)]=Veg_Optical_Parameter(OPT_PROP_H(i));
    [PFT_opt_L(i)]=Veg_Optical_Parameter(OPT_PROP_L(i));
end

% Specific leaf area of litter
%--------------------------------------------------------------------------
Sllit = 2 ; % Litter Specific Leaf area [m2 Litter / kg DM]

%% VEGETATION PART
%==========================================================================

%% High Vegetation
%--------------------------------------------------------------------------
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


fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve

% Phenology 
%--------------------------------------------------------------------------
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


for i=1:cc
    [Stoich_H(i)]=Veg_Stoichiometric_Parameter(Nl_H(i));
    [ParEx_H(i)]=Exudation_Parameter(0);
    [Mpar_H(i)]=Vegetation_Management_Parameter;
end



%% LOW VEGETATION
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

fbe_L = 1-fab_L; %% fraction below-ground sapwood and reserve

% Phenology 
%--------------------------------------------------------------------------
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

for i=1:cc
    [Stoich_L(i)]=Veg_Stoichiometric_Parameter(Nl_L(i));
    [ParEx_L(i)]=Exudation_Parameter(0);
    [Mpar_L(i)]=Vegetation_Management_Parameter;
end

%if ANSWER == 3
%    Mpar_L(1).jDay_cut=[180:243];
%    Mpar_L(1).LAI_cut=[-0.1]; %% LAI of grass after cut
%end



%%
Restating_parameters;

%if dbThick>5
%    nst=md_max;
%    k=(dbThick/5)^(1/(nst-1));
%    Zs_deb = [0 5*k.^(1:nst-1)]; %% [mm]
%    clear nst k
%else
%end


%% Restate debris parameters here, otherwise missing as output??
albs = table2array(TT_par(strcmp(TT_par.Parameters,'albs'),II));
lans = table2array(TT_par(strcmp(TT_par.Parameters,'lans'),II));
zoms = table2array(TT_par(strcmp(TT_par.Parameters,'zoms'),II));

%Deb_Par.alb= albs(s);
%Deb_Par.e_sur =  0.94;
%Deb_Par.lan = lans(s);
%Deb_Par.rho = 1496;  % [kg/m^3]
%Deb_Par.cs = 948;   % [J/kg K]
%Deb_Par.zom = zoms(s);
