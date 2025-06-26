%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% CONDITION INITIAL 1 STEP  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[]=INIT_COND_v3(num_cell,m_cell,n_cell,...
    cc_max,ms_max,md_max,...
    MASKn,GLH,Ca,SNOWD,SNOWALB,out, II_vect, ij_vect)

%% Debugging
%II_vect = cell2mat(POI.II);
%ij_vect = POI.ij;
%CP = 8; 
%cc_max = 3;

%% CODE
%%%%%%%%% TEMPORAL INITIALIZATION
An_H_t=     zeros(num_cell,cc_max,24);
An_L_t=     zeros(num_cell,cc_max,24);
O_t=        zeros(num_cell,ms_max,24);
PAR_t=      zeros(num_cell,24);
Pr_sno_t=   zeros(num_cell,24);
Psi_l_H_t=  zeros(num_cell,cc_max,24);
Psi_l_L_t=  zeros(num_cell,cc_max,24);
Psi_x_H_t=  zeros(num_cell,cc_max,24);
Psi_x_L_t=  zeros(num_cell,cc_max,24);
Rdark_H_t=  zeros(num_cell,cc_max,24);
Rdark_L_t=  zeros(num_cell,cc_max,24);
Ta_t=       zeros(num_cell,24);
Tdp_H_t=    zeros(num_cell,cc_max,24);
Tdp_L_t=    zeros(num_cell,cc_max,24);
Tdp_t=      zeros(num_cell,ms_max,24);
V_t=        zeros(num_cell,ms_max,24);

%%%%%%%%% SOIL MOISTUREPER LAYER
% "vi" has to correspond to "Zs" and "dz" in "PARAMETERS_SOIL_....m" (->Soil layer depths)
vi= [0    10    20    50   100   150   200   300   400    700   1000];
  
CP= 8;  %%Number of Carbon Pool (defined in "VEGETATION_MODULE_PAR.m"


%% GENERAL VEGETATION / HYDROLOGY
%==========================================================================
%
%==========================================================================
alp_soiltm1=    0*MASKn; % Soil relative humidity [-]
b_soiltm1=      0*MASKn; % Soil resistance beta factor [-]
Bamtm1=         0*MASKn;
Bemtm1=         0*MASKn;
BLittm1=        zeros(num_cell,cc_max); % Total litter content on the surface [kg DM /m2]
Ccrown_t_tm1=   ones(num_cell,cc_max); % Time dynamic vegetated fraction for each Ccrown [-]
Cicetm1=        0*MASKn; % Presence [1] or absence of ice [0]
Cicewtm1=       0*MASKn; % Presence [1] or absence of frozen water [0]
CK1tm1=         0*MASKn; % Check on Mass Balance 1 [mm/h]
Csnotm1=        0*MASKn; % Presence [1] or absence of snow [0]
Csnowtm1=       0*MASKn; % Presence [1] or absence of snow above frozen water [0]
dQ_Stm1=        0*MASKn; % Residual from energy budget [W/m2]
DQ_Stm1=        0*MASKn; % Residual from energy budget [W/m2]
dQVEGtm1=       0*MASKn; % Residual from energy budget of snow free vegetation [W/m2]
DT_Stm1=        0*MASKn; % Residual temperature difference in the energy budget [�C]
dw_SNOtm1=      0*MASKn; % Fraction of leaf covered by snow [-]
e_snotm1=       0.97*MASKn; % Emissivity of the snow [-]
EGtm1=          0*MASKn; % Evaporation from Bare soil [mm/h]
EICEtm1=        0*MASKn; % Evaporation/sublimation from Ice [mm/h]
EIn_rocktm1=    0*MASKn; % Evaporation from rocks [mm/h]
EIn_urbtm1=     0*MASKn; % Evaporation from impervious surface (e.g. roads; Debris surface) [mm/h]
EKtm1=          0*MASKn; % Cumulate Kinetic Energy of Precipitation [J/mm2]
ELittertm1=     0*MASKn; % Evaporation from the Litter [mm/h]
ertm1=          0*MASKn; % Splash eriosion [kg/h m2]
ESN_Intm1=      0*MASKn; % Evaporation from intercepted snow [mm/h]
ESNtm1=         0*MASKn; % Evaporation from the snowpack at the ground [mm/h]
SSNtm1=          0*MASKn; % Evaporation from the snowpack at the ground [mm/h]
SSN_Intm1=       0*MASKn; % Snow surface sublimation [mm/h]
EWATtm1=        0*MASKn; % Evaporation from water and ponds [mm/h]
FROCKtm1=       0*MASKn; % Water Storage in fractured rocks [mm]
ftm1=           0*MASKn; 
Gfintm1=        0*MASKn; % Ground Heat Flux heat diffusion [W/m2]
Gtm1=           0*MASKn; % Ground Heat Flux force restore method [W/m2]
Htm1=           0*MASKn; % Sensible Heat Flux [W/m2]
HVtm1=          0*MASKn; % Sensible Heat Flux from vegetation in presence of snow, Csnow=1 [W/m2]
ICE_Dtm1=       reshape(GLH,num_cell,1); % Ice thickness [m]
ICEym1=1000*0.916*reshape(GLH,num_cell,1).*MASKn; % Ice water equivalent [mm]
ICEtm1=         1000*0.916*reshape(GLH,num_cell,1); % Ice water equivalent [mm]
Imelttm1=       0*MASKn; % Ice melt [mm]
In_Htm1=        0*ones(num_cell,cc_max); % Intercepted water (storage) High Vegetation [mm]
In_Littertm1=   0*MASKn; % Intercepted water in Litter [mm]
In_Ltm1=        0*ones(num_cell,cc_max); % Intercepted water (storage) Low Vegetation [mm]
In_rocktm1=     0*MASKn; % Intercepted water (storage) Rocks [mm]
In_SWEtm1=      0*MASKn; % Intercepted snow water equivalent (storage) [mm]
In_urbtm1=      0*MASKn; % Intercepted water (storage) in impervious surface (e.g. roads; Debris surface) [mm]
IP_wctm1=       0*MASKn; % Ice pack water content [mm]
Lk_rocktm1=     0*MASKn; % Leakage rock surface to bedrock (recharge) [mm/h]
Lk_wattm1=      0*MASKn; % Leakage water pond to bedrock (recharge) [mm/h]
Lktm1=          0*MASKn; % Bottom Leakage soil to bedrock (recharge) [mm/h]
Lphotm1=        0*MASKn; % Energy consumed in the photosynthesis process [W/m2]
NavlItm1=       zeros(num_cell,3); % Mineral Nutrient available in the soil mean of last 365 days for Nitrogen|Phosphorous|Potassium [gX/m2]
NDVItm1=        0*MASKn; % Normalized Difference Vegetation Index [-]
NIcetm1=        0*MASKn; % New formed Ice [mm]
NIn_SWEtm1=     0*MASKn; % New Intercepted Snow Water Equivalent [mm]
OFtm1=          0*MASKn; % Soil Moisture First Soil Layer [-]
Oicetm1=        zeros(num_cell,ms_max); % Frozen volumetric Water content [-]
OStm1=          0*MASKn; % Soil Moisture for Bare Evaporation Layers [-]
Otm1=           zeros(num_cell,ms_max); 
POTtm1=         zeros(num_cell,ms_max); % Soil Water Potential [mm]
Pr_liqtm1=      0*MASKn; % Liquid Precipitation [mm/h]
Pr_snotm1=      0*MASKn; % Solid (snow) Precipitation [mm/h]
Q_channel=      zeros(m_cell,n_cell); % Water in channels [mm]
Q_exit=         0;
q_runon=        0*MASKn; % Runon/water ponded [mm]
QEtm1=          0*MASKn; % Latent Heat [W/m2]
QEVtm1=         0*MASKn; % Latent Heat from vegetation in presence of snow, Csnow=1 [W/m2]
Qfmtm1=         0*MASKn; % Heat for freezing or melting [W/m2]
Qi_in=          zeros(num_cell,ms_max); % Incoming Lateral subsurface flow [mm/h]
Qi_in_Rout=     zeros(m_cell,n_cell,ms_max);
Qi_out=         zeros(num_cell,ms_max); % Outgoing Lateral subsurface flow [mm/h]
Qi_out_Rout=    zeros(m_cell,n_cell,ms_max);
Qi_outtm1=      zeros(num_cell,ms_max);
Qsub_exit=      0;
Qvtm1=          0*MASKn; % Heat advected by Precipitation [W/m2]
r_littertm1=    zeros(num_cell,cc_max); % Litter resistance [s/m]
r_soiltm1=      0*MASKn; % Soil resistance [s/m]
ratm1=          0*MASKn; % Aerodynamic resistance [s/m]
Rdtm1=          0*MASKn; % Saturation excess runoff [mm]
Rhtm1=          0*MASKn; % Infiltration excess runoff [mm]
Rntm1=          0*MASKn; % Net radiation [W/m2]
rostm1=         250*MASKn; % Snow density [kg/m3]
SE_rocktm1=     0*MASKn; % Runoff on rocks [mm]
SE_urbtm1=      0*MASKn; % Runoff on impervious surface (e.g. roads; Debris surface) [mm]
%Slo_head=       reshape(Slo_top_S,num_cell,1)*ones(1,ms_max);
Smelttm1=       0*MASKn; % Snow melt [mm]

SNDtm1=SNOWD;   SNDtm1(isnan(SNDtm1))=0;  
SNDtm1=         reshape(SNDtm1,num_cell,1);% Snow depth [m]
SWEym1=SNDtm1.*rostm1; % Snow water equivalent [mm]
snow_albedotm1=SNOWALB(:).*ones(num_cell,4);   snow_albedotm1(isnan(snow_albedotm1),:)=0;  % Initial Snow albedo from satellite data [-]
soil_albedotm1= 0.2*ones(num_cell,4);
surface_albedotm1=0.2*ones(num_cell,4);

SWEtm1=         SNDtm1.*rostm1; % Snow water equivalent [mm]
SP_wctm1=       SWEtm1*0.05; % Snow pack water content [mm]

SWE_avalanchedtm1=    0*MASKn; % Snow water equivalent of "avalanched" mass [mm]

t_slstm1=       0*MASKn; % Time since last snowfall [s]
tau_snotm1=     SNOWALB(:).*MASKn; % Here we provide the initial snow albedo, strange but for me (Achille) it was the only way to have it working
Tdamptm1=       0*MASKn; % Dampening depth soil temperature [�C]
Tdp_snowtm1=     zeros(num_cell,5);
Tdebtm1=        0*ones(num_cell,md_max-1); % Temperature of each Debris layer [�C]
Tdptm1=         0*ones(num_cell,ms_max); % % Soil Temperature of each soil layer [�C]
Ticetm1=        0*MASKn; % Ice Temperature in presence of an icepack (Upper temperature) [�C]
Tstm0=          0*MASKn; % Soil/snow Prognostic Temperature 0 for the energy balance [�C]
Tstm1=          0*MASKn; % Soil/snow Prognostic Temperature 1 for the energy balance [�C]
Ts_undertm1 =   NaN*Tstm1;
TsVEGtm1=       0*MASKn;
U_SWEtm1=       0*MASKn; % Unloaded snow water equivalent from intercepted snow [mm]
Vicetm1=        zeros(num_cell,ms_max); % Volume of frozen water stored in the soil layer [mm]
Vtm1=           zeros(num_cell,ms_max); % Volume of liquid water stored in the soil layer [mm]

for jk=1:ms_max
    Vtm1(:,jk)=vi(jk)*MASKn;
end

WATtm1=         0*MASKn; % Volume of water in the lakes/ponds [mm]
WIStm1=         0*MASKn; % Water flux incoming to the soil [mm]
WR_IPtm1=       0*MASKn; % Water released from the ice pack [mm]
WR_SPtm1=       0*MASKn; % Water released from the snow pack [mm]
Ws_undertm1=    1*MASKn; % Wind speed in the under-canopy [m/s]
ZWTtm1=         0*MASKn; % Water table depth [mm]



%% SPECIFICATIONS FOR HIGH & LOW VEGETATION 
%==========================================================================
%{

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AgeDL_Htm1=     zeros(num_cell,cc_max);
AgeDL_Ltm1=     zeros(num_cell,cc_max);
AgeL_Htm1=      zeros(num_cell,cc_max);
AgeL_Ltm1=      zeros(num_cell,cc_max);
AgePl_Htm1=     zeros(num_cell,cc_max);
AgePl_Ltm1=     zeros(num_cell,cc_max);
An_Htm1=        zeros(num_cell,cc_max);
An_Ltm1=        zeros(num_cell,cc_max);
ANPP_Htm1=      zeros(num_cell,cc_max);
ANPP_Ltm1=      zeros(num_cell,cc_max);
B_Htm1=         zeros(num_cell,cc_max,CP);
B_Ltm1=         zeros(num_cell,cc_max,CP);
BA_Htm1=        zeros(num_cell,cc_max);
BA_Ltm1=        zeros(num_cell,cc_max);
Bfac_dayHtm1=   zeros(num_cell,cc_max);
Bfac_dayLtm1=   zeros(num_cell,cc_max);
Bfac_weekHtm1=  ones(num_cell,cc_max);  %%>=1?
Bfac_weekLtm1=  ones(num_cell,cc_max);  %%>=1?
Citm1_shdH=     zeros(num_cell,cc_max);
Citm1_shdL=     zeros(num_cell,cc_max);
Citm1_sunH=     zeros(num_cell,cc_max);
Citm1_sunL=     zeros(num_cell,cc_max);
dflo_Htm1=      zeros(num_cell,cc_max);
dflo_Ltm1=      zeros(num_cell,cc_max);
Dr_Htm1=        zeros(num_cell,cc_max);
Dr_Ltm1=        zeros(num_cell,cc_max);
e_rel_Htm1=     ones(num_cell,cc_max);  %%>=1?
e_rel_Ltm1=     ones(num_cell,cc_max);  %%>=1?
e_relN_Htm1=	ones(num_cell,cc_max);  %%>=1?
e_relN_Ltm1=    ones(num_cell,cc_max);  %%>=1?
EIn_Htm1=       zeros(num_cell,cc_max);
EIn_Ltm1=       zeros(num_cell,cc_max);
fapar_Htm1=     zeros(num_cell,cc_max);
fapar_Ltm1=     zeros(num_cell,cc_max);
FNC_Htm1=       ones(num_cell,cc_max);  %%>=1?
FNC_Ltm1=       ones(num_cell,cc_max);  %%>=1?
gsr_Htm1=       zeros(num_cell,cc_max);
gsr_Ltm1=       zeros(num_cell,cc_max);
hc_Htm1=        zeros(num_cell,cc_max);
hc_Ltm1=        zeros(num_cell,cc_max);
In_Htm1=        zeros(num_cell,cc_max);
In_Ltm1=        zeros(num_cell,cc_max);
ISOIL_Htm1=     zeros(num_cell,cc_max,18);
ISOIL_Ltm1=     zeros(num_cell,cc_max,18);
Jsx_Htm1=       zeros(num_cell,cc_max);
Jsx_Ltm1=       zeros(num_cell,cc_max);
Jxl_Htm1=       zeros(num_cell,cc_max);
Jxl_Ltm1=       zeros(num_cell,cc_max);
Kleaf_Htm1=     zeros(num_cell,cc_max);
Kleaf_Ltm1=     zeros(num_cell,cc_max);
Kreserve_Htm1=	zeros(num_cell,cc_max);
Kreserve_Ltm1=	zeros(num_cell,cc_max);
Kuptake_Htm1=	zeros(num_cell,cc_max);
Kuptake_Ltm1=	zeros(num_cell,cc_max);
Kx_Htm1=        zeros(num_cell,cc_max);
Kx_Ltm1=        zeros(num_cell,cc_max);
LAI_Htm1=       zeros(num_cell,cc_max);
LAI_Ltm1=       zeros(num_cell,cc_max);
LAIdead_Htm1=   zeros(num_cell,cc_max);
LAIdead_Ltm1=   zeros(num_cell,cc_max);
ManIHtm1=       zeros(num_cell,cc_max);
ManILtm1=       zeros(num_cell,cc_max);
NBLeaf_Htm1=	zeros(num_cell,cc_max);
NBLeaf_Ltm1=	zeros(num_cell,cc_max);
NBLI_Htm1=      zeros(num_cell,cc_max);
NBLI_Ltm1=      zeros(num_cell,cc_max);
NPP_Htm1=       zeros(num_cell,cc_max);
NPP_Ltm1=       zeros(num_cell,cc_max);
NPPI_Htm1=      zeros(num_cell,cc_max);
NPPI_Ltm1=      zeros(num_cell,cc_max);
Nreserve_Htm1=  zeros(num_cell,cc_max);	
Nreserve_Ltm1=	zeros(num_cell,cc_max);
NuLit_Htm1=     zeros(num_cell,cc_max,3);
NuLit_Ltm1=     zeros(num_cell,cc_max,3);
NupI_Htm1=      zeros(num_cell,cc_max,3);
NupI_Ltm1=      zeros(num_cell,cc_max,3);
Nuptake_Htm1=	zeros(num_cell,cc_max);
Nuptake_Ltm1=	zeros(num_cell,cc_max);
OHtm1=          zeros(num_cell,cc_max);
OLtm1=          zeros(num_cell,cc_max);
PARI_Htm1=      zeros(num_cell,cc_max,3);
PARI_Ltm1=      zeros(num_cell,cc_max,3);
PHE_S_Htm1=     ones(num_cell,cc_max);    %%>=1, otherwise issues in PHENOLOGY_STATE.m
PHE_S_Ltm1=     ones(num_cell,cc_max);    %%>=1, otherwise issues in PHENOLOGY_STATE.m
Preserve_Htm1=	zeros(num_cell,cc_max);
Preserve_Ltm1=	zeros(num_cell,cc_max);
Psi_l_Htm1=     zeros(num_cell,cc_max);
Psi_l_Ltm1=     zeros(num_cell,cc_max);
Psi_s_Htm1=     zeros(num_cell,cc_max);
Psi_s_Ltm1=     zeros(num_cell,cc_max);
Psi_x_Htm1=     zeros(num_cell,cc_max);
Psi_x_Ltm1=     zeros(num_cell,cc_max);
Puptake_Htm1=	zeros(num_cell,cc_max);
Puptake_Ltm1=	zeros(num_cell,cc_max);
RA_Htm1=        zeros(num_cell,cc_max);
RA_Ltm1=        zeros(num_cell,cc_max);
rap_Htm1=       zeros(num_cell,cc_max);
rap_Ltm1=       zeros(num_cell,cc_max);
RB_Htm1=        zeros(num_cell,cc_max,7);
rb_Htm1=        zeros(num_cell,cc_max);
RB_Ltm1=        zeros(num_cell,cc_max,7);
rb_Ltm1=        zeros(num_cell,cc_max);
Rdark_Htm1=     zeros(num_cell,cc_max);
Rdark_Ltm1=     zeros(num_cell,cc_max);
Rexmy_Htm1=     zeros(num_cell,cc_max,3);
Rexmy_Ltm1=     zeros(num_cell,cc_max,3);
Rg_Htm1=        zeros(num_cell,cc_max);
Rg_Ltm1=        zeros(num_cell,cc_max);
rKc_Htm1=       zeros(num_cell,cc_max);
rKc_Ltm1=       zeros(num_cell,cc_max);
Rmc_Htm1=       zeros(num_cell,cc_max);
Rmc_Ltm1=       zeros(num_cell,cc_max);
Rmr_Htm1=       zeros(num_cell,cc_max);
Rmr_Ltm1=       zeros(num_cell,cc_max);
Rms_Htm1=       zeros(num_cell,cc_max);
Rms_Ltm1=       zeros(num_cell,cc_max);
rNc_Htm1=       zeros(num_cell,cc_max);
rNc_Ltm1=       zeros(num_cell,cc_max);
rPc_Htm1=       zeros(num_cell,cc_max);
rPc_Ltm1=       zeros(num_cell,cc_max);
Rrootl_Htm1=    zeros(num_cell,cc_max);
Rrootl_Ltm1=    zeros(num_cell,cc_max);
rs_shdHtm1=     zeros(num_cell,cc_max);
rs_shdLtm1=     zeros(num_cell,cc_max);
rs_sunHtm1=     zeros(num_cell,cc_max);
rs_sunLtm1=     zeros(num_cell,cc_max);
SAI_Htm1=       zeros(num_cell,cc_max);
SAI_Ltm1=       zeros(num_cell,cc_max);
Sfr_Htm1=       zeros(num_cell,cc_max);
Sfr_Ltm1=       zeros(num_cell,cc_max);
SIF_Htm1=       zeros(num_cell,cc_max);
SIF_Ltm1=       zeros(num_cell,cc_max);
Slf_Htm1=       zeros(num_cell,cc_max);
Slf_Ltm1=       zeros(num_cell,cc_max);
Sll_Htm1=       zeros(num_cell,cc_max);
Sll_Ltm1=       zeros(num_cell,cc_max);
Sr_Htm1=        zeros(num_cell,cc_max);
Sr_Ltm1=        zeros(num_cell,cc_max);
SupK_Htm1=      zeros(num_cell,cc_max);
SupK_Ltm1=      zeros(num_cell,cc_max);
SupN_Htm1=      zeros(num_cell,cc_max);
SupN_Ltm1=      zeros(num_cell,cc_max);
SupP_Htm1=      zeros(num_cell,cc_max);
SupP_Ltm1=      zeros(num_cell,cc_max);
Swm_Htm1=       zeros(num_cell,cc_max);
Swm_Ltm1=       zeros(num_cell,cc_max);
T_Htm1=         zeros(num_cell,cc_max);
T_Ltm1=         zeros(num_cell,cc_max);
TBio_Htm1=      zeros(num_cell,cc_max);
TBio_Ltm1=      zeros(num_cell,cc_max);
Tden_Htm1=      zeros(num_cell,cc_max);
Tden_Ltm1=      zeros(num_cell,cc_max);
Tdp_Htm1=       zeros(num_cell,cc_max);
Tdp_Ltm1=       zeros(num_cell,cc_max);
TdpI_Htm1=      zeros(num_cell,cc_max);
TdpI_Ltm1=      zeros(num_cell,cc_max);
TexC_Htm1=      zeros(num_cell,cc_max);
TexC_Ltm1=      zeros(num_cell,cc_max);
TexK_Htm1=      zeros(num_cell,cc_max);
TexK_Ltm1=      zeros(num_cell,cc_max);
TexN_Htm1=      zeros(num_cell,cc_max);
TexN_Ltm1=      zeros(num_cell,cc_max);
TexP_Htm1=      zeros(num_cell,cc_max);
TexP_Ltm1=      zeros(num_cell,cc_max);
TNIT_Htm1=      zeros(num_cell,cc_max);
TNIT_Ltm1=      zeros(num_cell,cc_max);
TPHO_Htm1=      zeros(num_cell,cc_max);
TPHO_Ltm1=      zeros(num_cell,cc_max);
TPOT_Htm1=      zeros(num_cell,cc_max);
TPOT_Ltm1=      zeros(num_cell,cc_max);
Vl_Htm1=        zeros(num_cell,cc_max);
Vl_Ltm1=        zeros(num_cell,cc_max);
Vx_Htm1=        zeros(num_cell,cc_max);
Vx_Ltm1=        zeros(num_cell,cc_max);

%% CLASSIFICATION OF STUDY AREA
%==========================================================================
% 1 Fir (everg.) H  |  2 Larch (decid.) H 
% 3 Grass C3     L  |  4 Shrub (decid.) L
% 5 Grass C3     H  |  6 Shrub (decid.) H
%
% ksv vector maps CORINE classification to the aggregated classes from the 
% file Legendkey. 
% The aggregated classes are 8
% 1 to 6 are vegetation.
% 7: Rock
% 8: Water
% 9: Urban
%==========================================================================
%disp('here')

%k = 90
for k = 1:length(ij_vect)
%disp(k)
valid = find(II_vect(k,:)); %Vector with valid parameters

if ~isempty(find(II_vect(k,:))) 
% if isempty(find(II_vect(k,:) == 1)) is empty, it means I am in a pixel
% that is a rock, water or bare. It is not vegetation, and therefore, I do
% not assign any parameter of vegetation to that pixel and the next lines
% are skipped. 

%% AgeL: Leaf Age [days] (FROM REFERENCE)
% OLD categories   [fir    larch   grass  shrub  BLever    BLdec  NoVeg]
% categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low   BLdec_high NoVeg]  
AgeL_class_H = [979         0         0          0        0           0         0          0        305          80          0          NaN];
AgeL_class_L = [0           0         0          0        0           0         160        644      0            89          0          NaN];

AgeL_Htm1(ij_vect(k),1:length(valid)) = AgeL_class_H(valid);
AgeL_Ltm1(ij_vect(k),1:length(valid)) = AgeL_class_L(valid);

%AgeL_Htm1(ksv==1,1)=979; 
%AgeL_Htm1(ksv==2,1)=220; 
%AgeL_Ltm1(ksv==3,1)=0; 
%AgeL_Ltm1(ksv==4,1)=0; 
%AgeL_Htm1(ksv==5,1)=305; 
%AgeL_Htm1(ksv==6,1)=0; 

%% AgeDL: Dead leaf Age [days] (NO REFERENCE)
% categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low   BLdec_high NoVeg]  
AgeDL_class_H = [0           0         0          0        0            0         0          0        0           0           0          NaN  ]; 
AgeDL_class_L = [0           0         0          0        0            0         0          0        0           0           0          NaN  ]; 

AgeDL_Htm1(ij_vect(k),1:length(valid)) = AgeDL_class_H(valid);
AgeDL_Ltm1(ij_vect(k),1:length(valid)) = AgeDL_class_L(valid);

%disp(AgeDL_Htm1(ij,:))

%AgeDL_Htm1(ksv==1,1)=0; 
%AgeDL_Htm1(ksv==2,1)=0; 
%AgeDL_Ltm1(ksv==3,1)=0; 
%AgeDL_Ltm1(ksv==4,1)=0; 
%AgeDL_Htm1(ksv==5,1)=0; 
%AgeDL_Htm1(ksv==6,1)=0; 

%% AgePl: Age of the forest stand or plantation [days] (NO REFERENCE)
% categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
AgePl_class_H = [0           0         0          0        0           0         0          0        0            0          0          NaN  ]; 
AgePl_class_L = [0           0         0          0        0           0         0          0        0            0          0          NaN  ]; 

AgePl_Htm1(ij_vect(k),1:length(valid)) = AgePl_class_H(valid);
AgePl_Ltm1(ij_vect(k),1:length(valid)) = AgePl_class_L(valid);

%AgePl_Htm1(ksv==1,1)=0; 
%AgePl_Htm1(ksv==2,1)=0; 
%AgePl_Ltm1(ksv==3,1)=0; 
%AgePl_Ltm1(ksv==4,1)=0;
%AgePl_Htm1(ksv==5,1)=0;
%AgePl_Htm1(ksv==6,1)=0;

%% Bfac_week: Plant stress factor integrated at the weekly scale [0-1] (NO REFERENCE)
% categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Bfac_class_H =  [0           0         0          0        0           0         0          0        0            0          0          NaN  ]; 
Bfac_class_L =  [0           0         0          0        0           0         0          0        0            0          0          NaN  ]; 

Bfac_weekHtm1(ij_vect(k),1:length(valid)) = Bfac_class_H(valid);
Bfac_weekLtm1(ij_vect(k),1:length(valid)) = Bfac_class_L(valid);

%Bfac_weekHtm1(ksv==1,1)=0; 
%Bfac_weekHtm1(ksv==2,1)=0; 
%Bfac_weekLtm1(ksv==3,1)=0; 
%Bfac_weekLtm1(ksv==4,1)=0;
%Bfac_weekHtm1(ksv==5,1)=0;
%Bfac_weekHtm1(ksv==6,1)=0;

%% Ci_sun: CO2 sunlit leaf internal concentration [umolCO2/mol] (NO REFERENCE)
% categories     [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Citm1_class_H =  [Ca(1)       0         0          0        0           0         0          0        Ca(1)        Ca(1)      Ca(1)      NaN  ]; 
Citm1_class_L =  [0           Ca(1)     Ca(1)      Ca(1)    Ca(1)       Ca(1)     Ca(1)      Ca(1)    0            0          0          NaN  ]; 

Citm1_sunH(ij_vect(k),1:length(valid)) = Citm1_class_H(valid);
Citm1_sunL(ij_vect(k),1:length(valid)) = Citm1_class_L(valid);

%Citm1_sunH(ksv==1,1)=Ca(1); 
%Citm1_sunH(ksv==2,1)=Ca(1); 
%Citm1_sunL(ksv==3,1)=Ca(1); 
%Citm1_sunL(ksv==4,1)=Ca(1);  
%Citm1_sunH(ksv==5,1)=Ca(1);  
%Citm1_sunH(ksv==6,1)=Ca(1);  

%% Ci_shd:  CO2 sunlit leaf internal concentration [umolCO2/mol] (NO REFERENCE)
Citm1_shdH(ij_vect(k),1:length(valid)) = Citm1_class_H(valid);
Citm1_shdL(ij_vect(k),1:length(valid)) = Citm1_class_L(valid);

%Citm1_shdH(ksv==1,1)=Ca(1); 
%Citm1_shdH(ksv==2,1)=Ca(1); 
%Citm1_shdL(ksv==3,1)=Ca(1); 
%Citm1_shdL(ksv==4,1)=Ca(1); 
%Citm1_shdH(ksv==5,1)=Ca(1); 
%Citm1_shdH(ksv==6,1)=Ca(1); 

%% dflo:  Days from leaf onset [days] (FROM REFERENCE)
% categories     [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
dflo_class_H =   [1           0         0          0        0           0         0          0        1            102        0          NaN  ]; 
dflo_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ]; 

dflo_Htm1(ij_vect(k),1:length(valid)) = dflo_class_H(valid);
dflo_Ltm1(ij_vect(k),1:length(valid)) = dflo_class_L(valid);

%dflo_Htm1(ksv==1,1)=1; 
%dflo_Htm1(ksv==2,1)=1; 
%dflo_Ltm1(ksv==3,1)=0; 
%dflo_Ltm1(ksv==4,1)=0;   
%dflo_Htm1(ksv==5,1)=1;   
%dflo_Htm1(ksv==6,1)=1;   


%% e_rel: Relative Efficiency of the photosynthesis apparatus due to Age/Day-length (FROM REFERENCE)
% categories     [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
e_rel_class_H =  [1           0         0           0       0           0         0          0        1            1          1          NaN  ]; 
e_rel_class_L =  [0           1         1           1       1           1         1          1        0            1          1          NaN  ]; 

e_rel_Htm1(ij_vect(k),1:length(valid)) = e_rel_class_H(valid);
e_rel_Ltm1(ij_vect(k),1:length(valid)) = e_rel_class_L(valid);

%e_rel_Htm1(ksv==1,1)=1; 
%e_rel_Htm1(ksv==2,1)=1; 
%e_rel_Ltm1(ksv==3,1)=1; 
%e_rel_Ltm1(ksv==4,1)=1;   
%e_rel_Htm1(ksv==5,1)=1;   
%e_rel_Htm1(ksv==6,1)=1;   


%% e_relN:       Relative Efficiency of the photosynthesis apparatus due to N Limitations (NO REFERENCE)
% categories      [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
e_relN_class_H =  [1           0         0          0        0           0         0          0        1            0          1          NaN  ]; 
e_relN_class_L =  [0           1         1          1        1           1         1          1        0            0          0          NaN  ]; 

e_relN_Htm1(ij_vect(k),1:length(valid)) = e_relN_class_H(valid);
e_relN_Ltm1(ij_vect(k),1:length(valid)) = e_relN_class_L(valid);

%e_relN_Htm1(ksv==1,1)=1; 
%e_relN_Htm1(ksv==2,1)=1; 
%e_relN_Ltm1(ksv==3,1)=1; 
%e_relN_Ltm1(ksv==4,1)=1;  
%e_relN_Htm1(ksv==5,1)=1; 
%e_relN_Htm1(ksv==6,1)=1;   


%% FNC: Nitrogen Stress Factor for vegetation [0-1] (FROM REFERENCE)
% categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
FNC_class_H =   [1           1         1          1        1           0         1          0        1            1          1          NaN  ];
FNC_class_L =   [0           1         1          1        1           1         1          1        0            1          1          NaN  ];

FNC_Htm1(ij_vect(k),1:length(valid)) = FNC_class_H(valid);
FNC_Ltm1(ij_vect(k),1:length(valid)) = FNC_class_L(valid);

%FNC_Htm1(ksv==1,1)=1; 
%FNC_Htm1(ksv==2,1)=1; 
%FNC_Ltm1(ksv==3,1)=1; 
%FNC_Ltm1(ksv==4,1)=1;  
%FNC_Htm1(ksv==5,1)=1; 
%FNC_Htm1(ksv==6,1)=1; 

%% hc: Vegetation Height [m] (FROM REFERENCE)
% categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R      grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
hc_class_H =    [10          0         0          0        0            0         0          0.0      15           13.5       22         NaN  ]; 
hc_class_L =    [0           0.15      0.15       0.15     0.15         0.1       0.15       0.5      0            0.1        0          NaN  ]; 

hc_Htm1(ij_vect(k), 1:length(valid)) = hc_class_H(valid);
hc_Ltm1(ij_vect(k), 1:length(valid)) = hc_class_L(valid);

%hc_Htm1(ksv==1,1)=10; 
%hc_Htm1(ksv==2,1)=10; 
%hc_Ltm1(ksv==3,1)=0.1; 
%hc_Ltm1(ksv==4,1)=0.5;         
%hc_Htm1(ksv==5,1)=15; 
%hc_Htm1(ksv==6,1)=15;         

%% Kreserve: Mobile Reserve of Potassium in vegetation [gK/m2 PFT] (FROM REFERENCE)
% categories         [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Kreserve_class_H =   [1000        0         0          0        0           1000      1000       1000     1000         1000       1000       NaN  ]; 
Kreserve_class_L =   [1000        1000      1000       1000     1000        1000      1000       1000     1000         1000       1000       NaN  ]; 

Kreserve_Htm1(ij_vect(k),1:length(valid)) = Kreserve_class_H(valid);
Kreserve_Ltm1(ij_vect(k),1:length(valid)) = Kreserve_class_L(valid);

%Kreserve_Htm1(ksv==1,1)=1000; 
%Kreserve_Htm1(ksv==2,1)=1000; 
%Kreserve_Ltm1(ksv==3,1)=1000; 
%Kreserve_Ltm1(ksv==4,1)=1000;
%Kreserve_Htm1(ksv==5,1)=1000; 
%Kreserve_Htm1(ksv==6,1)=1000; 

%% LAI: Leaf area index [-] (NO REFERENCE)
% categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
LAI_class_H =   [3           0         0          0        0           0         0          0        2            0          0          NaN  ]; 
LAI_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

LAI_Htm1(ij_vect(k),1:length(valid)) = LAI_class_H(valid);
LAI_Ltm1(ij_vect(k),1:length(valid)) = LAI_class_L(valid);

%LAI_Htm1(ksv==1,1)=3;
%LAI_Htm1(ksv==2,1)=2;
%LAI_Ltm1(ksv==3,1)=0;
%LAI_Ltm1(ksv==4,1)=0;
%LAI_Htm1(ksv==5,1)=2;
%LAI_Htm1(ksv==6,1)=0;

%% NBLeaf: New Leaf Biomass [gC/m2 day] (NO REFERENCE)
% categories       [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
NBLeaf_class_H =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ]; 
NBLeaf_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

NBLeaf_Htm1(ij_vect(k),1:length(valid)) = NBLeaf_class_H(valid);
NBLeaf_Ltm1(ij_vect(k),1:length(valid)) = NBLeaf_class_L(valid);

%NBLeaf_Htm1(ksv==1,1)=0; 
%NBLeaf_Htm1(ksv==2,1)=0; 
%NBLeaf_Ltm1(ksv==3,1)=0; 
%NBLeaf_Ltm1(ksv==4,1)=0; 
%NBLeaf_Htm1(ksv==5,1)=0; 
%NBLeaf_Htm1(ksv==6,1)=0; 

%% NBLI: Integral of New Leaf Biomass over 30 days [gC/m2 day] (NO REFERENCE)
% categories       [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
NBLI_class_H =     [0           0         0          0        0           0         0          0        0            0          0          NaN  ]; 
NBLI_class_L =     [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

NBLI_Htm1(ij_vect(k),1:length(valid)) = NBLI_class_H(valid);
NBLI_Ltm1(ij_vect(k),1:length(valid)) = NBLI_class_L(valid);

%NBLI_Htm1(ksv==1,1)=0; 
%NBLI_Htm1(ksv==2,1)=0; 
%NBLI_Ltm1(ksv==3,1)=0; 
%NBLI_Ltm1(ksv==4,1)=0; 
%NBLI_Htm1(ksv==5,1)=0;     
%NBLI_Htm1(ksv==6,1)=0;     

%% NPP:  Net Primary Production [gC/m2 PFT day] (NO REFERENCE)
% categories    [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
NPP_class_H =   [0           0         0          0         0          0         0          0        0            0          0          NaN  ]; 
NPP_class_L =   [0           0         0          0         0          0         0          0        0            0          0          NaN  ];

NPP_Htm1(ij_vect(k),1:length(valid)) = NPP_class_H(valid);
NPP_Ltm1(ij_vect(k),1:length(valid)) = NPP_class_L(valid);

%NPP_Htm1(ksv==1,1)=0; 
%NPP_Htm1(ksv==2,1)=0; 
%NPP_Ltm1(ksv==3,1)=0; 
%NPP_Ltm1(ksv==4,1)=0;  
%NPP_Htm1(ksv==5,1)=0;    
%NPP_Htm1(ksv==6,1)=0;    

%% NPII: Integral of Net Primary Production over 7 days [gC/m2 PFT day] (NO REFERENCE)
% categories     [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
NPPI_class_H =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ]; 
NPPI_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

NPPI_Htm1(ij_vect(k),1:length(valid)) = NPPI_class_H(valid);
NPPI_Ltm1(ij_vect(k),1:length(valid)) = NPPI_class_L(valid);

%NPPI_Htm1(ksv==1,1)=0; 
%NPPI_Htm1(ksv==2,1)=0; 
%NPPI_Ltm1(ksv==3,1)=0; 
%NPPI_Ltm1(ksv==4,1)=0;
%NPPI_Htm1(ksv==5,1)=0;
%NPPI_Htm1(ksv==6,1)=0;

%% Nreserve: Mobile Reserve of Nitrogen in vegetation [gN/m2 PFT] (FROM REFERENCE)
% categories        [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R      grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Nreserve_class_H =  [1000        0         0          0        0            0         0          0        1000         1000       1000       NaN  ];
Nreserve_class_L =  [0           1000      1000       1000     1000         1000       1000       1000     0           1000       0          NaN  ];

Nreserve_Htm1(ij_vect(k),1:length(valid)) = Nreserve_class_H(valid);
Nreserve_Ltm1(ij_vect(k),1:length(valid)) = Nreserve_class_L(valid);

%Nreserve_Htm1(ksv==1,1)=1000; 
%Nreserve_Htm1(ksv==2,1)=1000; 
%Nreserve_Ltm1(ksv==3,1)=1000; 
%Nreserve_Ltm1(ksv==4,1)=1000; 
%Nreserve_Htm1(ksv==5,1)=1000; 
%Nreserve_Htm1(ksv==6,1)=1000; 

%% PHE: Phenology State [#] (FROM REFERENCE)
% categories      [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
PHE_S_class_H =   [1           0         0          0        0           0         0          0        1            3          1          NaN  ];
PHE_S_class_L =   [0           1         1          1        1           1         1          1        0            3          0          NaN  ];

PHE_S_Htm1(ij_vect(k),1:length(valid)) = PHE_S_class_H(valid);
PHE_S_Ltm1(ij_vect(k),1:length(valid)) = PHE_S_class_L(valid);

%PHE_S_Htm1(ksv==1,1)=1; 
%PHE_S_Htm1(ksv==2,1)=1; 
%PHE_S_Ltm1(ksv==3,1)=1; 
%PHE_S_Ltm1(ksv==4,1)=1;
%PHE_S_Htm1(ksv==5,1)=1;
%PHE_S_Htm1(ksv==6,1)=1;

%% Preserve: Mobile Reserve of Phosporus in vegetation [gP/m2 PFT] (FROM REFERENCE)
% categories        [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Preserve_class_H =  [1000        0         0          0        0           0         0          0        1000         1000       1000       NaN  ];
Preserve_class_L =  [0           1000      1000       1000     1000        1000      1000       1000     0            1000       0          NaN  ];

Preserve_Htm1(ij_vect(k),1:length(valid)) = Preserve_class_H(valid);
Preserve_Ltm1(ij_vect(k),1:length(valid)) = Preserve_class_L(valid);

%Preserve_Htm1(ksv==1,1)=1000; 
%Preserve_Htm1(ksv==2,1)=1000; 
%Preserve_Ltm1(ksv==3,1)=1000; 
%Preserve_Ltm1(ksv==4,1)=1000;
%Preserve_Htm1(ksv==5,1)=1000;
%Preserve_Htm1(ksv==6,1)=1000;

%% Psi_x:        Soil water potential in the stem xylem [MPa] (NO REFERENCE)
% categories      [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Psi_x_class_H =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];
Psi_x_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

Psi_x_Htm1(ij_vect(k),1:length(valid)) = Psi_x_class_H(valid);
Psi_x_Ltm1(ij_vect(k),1:length(valid)) = Psi_x_class_L(valid);

%Psi_x_Htm1(ksv==1,1)=0; 
%Psi_x_Htm1(ksv==2,1)=0; 
%Psi_x_Ltm1(ksv==3,1)=0; 
%Psi_x_Ltm1(ksv==4,1)=0; 
%Psi_x_Htm1(ksv==5,1)=0; 
%Psi_x_Htm1(ksv==6,1)=0; 

%% Psi_l:        Leaf water potential in the leaves [MPa] (NO REFERENCE)
% categories      [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Psi_l_class_H =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];
Psi_l_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

Psi_l_Htm1(ij_vect(k),1:length(valid)) = Psi_l_class_H(valid);
Psi_l_Ltm1(ij_vect(k),1:length(valid)) = Psi_l_class_L(valid);

%Psi_l_Htm1(ksv==1,1)=0; 
%Psi_l_Htm1(ksv==2,1)=0; 
%Psi_l_Ltm1(ksv==3,1)=0; 
%Psi_l_Ltm1(ksv==4,1)=0;    
%Psi_l_Htm1(ksv==5,1)=0;    
%Psi_l_Htm1(ksv==6,1)=0;    

%% SAI: Stem area index [-] (FROM REFERENCE)
% categories      [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
SAI_class_H =     [0.1         0         0          0        0            0         0          0        0.1         0.2        0.2        NaN  ];
SAI_class_L =     [0           0.001     0.001      0.001    0.001        0.01      0.01       0.01     0           0.001      0          NaN  ];

SAI_Htm1(ij_vect(k),1:length(valid)) = SAI_class_H(valid);
SAI_Ltm1(ij_vect(k),1:length(valid)) = SAI_class_L(valid);

%SAI_Htm1(ksv==1,1)=0.1; 
%SAI_Htm1(ksv==2,1)=0.1; 
%SAI_Ltm1(ksv==3,1)=0.01; 
%SAI_Ltm1(ksv==4,1)=0.01;
%SAI_Htm1(ksv==5,1)=0.1;
%SAI_Htm1(ksv==6,1)=0.2;

%% TBio: Total standing biomass temporally variable [ton DM / ha ] (FROM REFERENCE)
% categories      [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low   BLdec_high NoVeg]  
TBio_class_H =    [1           0         0          0        0           0         0          0        1            200         1          NaN  ];
TBio_class_L =    [0           1         1          1        1           1         1          1        0            1           0          NaN  ];

TBio_Htm1(ij_vect(k),1:length(valid)) = TBio_class_H(valid);
TBio_Ltm1(ij_vect(k),1:length(valid)) = TBio_class_L(valid);

%TBio_Htm1(ksv==1,1)=1; 
%TBio_Htm1(ksv==2,1)=1; 
%TBio_Ltm1(ksv==3,1)=1; 
%TBio_Ltm1(ksv==4,1)=1; 
%TBio_Htm1(ksv==5,1)=1; 
%TBio_Htm1(ksv==6,1)=1; 

%% Tden: Tree density, High vegetation [n*ind/ ha] (NO REFERENCE)
% categories     [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Tden_class_H =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];
Tden_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

Tden_Htm1(ij_vect(k),1:length(valid)) = Tden_class_H(valid);
Tden_Ltm1(ij_vect(k),1:length(valid)) = Tden_class_L(valid);

%Tden_Htm1(ksv==1,1)=0; 
%Tden_Htm1(ksv==2,1)=0; 
%Tden_Ltm1(ksv==3,1)=0; 
%Tden_Ltm1(ksv==4,1)=0; 
%Tden_Htm1(ksv==5,1)=0; 
%Tden_Htm1(ksv==6,1)=0; 

%% Vx:           Water Volume in the xylem [mm m2 ground/m2 PFT] (FROM REFERENCE)
% categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Vx_class_H =   [0           0         0          0        0           0         0          0        0            10         1          NaN  ];
Vx_class_L =   [0           0         0          0        0           0         0          0        0            10         1          NaN  ];

Vx_Htm1(ij_vect(k),1:length(valid)) = Vx_class_H(valid);
Vx_Ltm1(ij_vect(k),1:length(valid)) = Vx_class_L(valid);

%Vx_Htm1(ksv==1,1)=0; 
%Vx_Htm1(ksv==2,1)=0; 
%Vx_Ltm1(ksv==3,1)=0; 
%Vx_Ltm1(ksv==4,1)=0;   
%Vx_Htm1(ksv==5,1)=0;          
%Vx_Htm1(ksv==6,1)=0;          

%% Vl:           Water Volume in the leaves [mm m2 ground/m2 PFT] (FROM REFERENCE)
% categories   [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R    grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
Vl_class_H =   [0           0         0          0        0          0         0          0        0            10         1          NaN  ];
Vl_class_L =   [0           1         1          1        1          0         0          0        0            1          1          NaN  ];

Vl_Htm1(ij_vect(k),1:length(valid)) = Vl_class_H(valid);
Vl_Ltm1(ij_vect(k),1:length(valid)) = Vl_class_L(valid);

%Vl_Htm1(ksv==1,1)=0; 
%Vl_Htm1(ksv==2,1)=0; 
%Vl_Ltm1(ksv==3,1)=0; 
%Vl_Ltm1(ksv==4,1)=0;
%Vl_Htm1(ksv==5,1)=0;
%Vl_Htm1(ksv==6,1)=0;

%% NupI:         Integrated Nitrogen|Phosphorous|Potassium uptake over 365 days [gX/m2 PFT day] (FROM REFERENCE)
% categories     [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
NupI_class_H =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];
NupI_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

NupI_Htm1(ij_vect(k),1:length(valid),1) = NupI_class_H(valid);
NupI_Ltm1(ij_vect(k),1:length(valid),1) = NupI_class_L(valid);

NupI_Htm1(ij_vect(k),1:length(valid),2) = NupI_class_H(valid);
NupI_Ltm1(ij_vect(k),1:length(valid),2) = NupI_class_L(valid);

NupI_Htm1(ij_vect(k),1:length(valid),3) = NupI_class_H(valid);
NupI_Ltm1(ij_vect(k),1:length(valid),3) = NupI_class_L(valid);

%disp(NupI_Ltm1(ij,:,3))

%NupI_Htm1(ksv==1,1,1)= 0; NupI_Htm1(ksv==1,1,2)= 0; NupI_Htm1(ksv==1,1,3)= 0; 
%NupI_Htm1(ksv==2,1,1)= 0; NupI_Htm1(ksv==2,1,2)= 0; NupI_Htm1(ksv==2,1,3)= 0;
%NupI_Ltm1(ksv==3,1,1)= 0; NupI_Ltm1(ksv==3,1,2)= 0; NupI_Ltm1(ksv==3,1,3)= 0;
%NupI_Ltm1(ksv==4,1,1)= 0; NupI_Ltm1(ksv==4,1,2)= 0; NupI_Ltm1(ksv==4,1,3)= 0;
%NupI_Htm1(ksv==5,1,1)= 0; NupI_Htm1(ksv==5,1,2)= 0; NupI_Htm1(ksv==5,1,3)= 0;
%NupI_Htm1(ksv==6,1,1)= 0; NupI_Htm1(ksv==6,1,2)= 0; NupI_Htm1(ksv==6,1,3)= 0;

%% NupLit:       Nitrogen|Phosphorous|Potassium content in the litter [gN/m2 PFT] (NO REFERENCE)
% categories       [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R      grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
NupLit_class_H =   [0           0         0          0        0            0         0          0        0            0          0          NaN  ];
NupLit_class_L =   [0           0         0          0        0            0         0          0        0            0          0          NaN  ];

NupLit_Htm1(ij_vect(k),1:length(valid),1) = NupLit_class_H(valid);
NupLit_Ltm1(ij_vect(k),1:length(valid),1) = NupLit_class_L(valid);

NupLit_Htm1(ij_vect(k),1:length(valid),2) = NupLit_class_H(valid);
NupLit_Ltm1(ij_vect(k),1:length(valid),2) = NupLit_class_L(valid);

NupLit_Htm1(ij_vect(k),1:length(valid),3) = NupLit_class_H(valid);
NupLit_Ltm1(ij_vect(k),1:length(valid),3) = NupLit_class_L(valid);

%NupLit_Htm1(ksv==1,1,1)= 0; NupLit_Htm1(ksv==1,1,2)= 0; NupLit_Htm1(ksv==1,1,3)= 0; 
%NupLit_Htm1(ksv==2,1,1)= 0; NupLit_Htm1(ksv==2,1,2)= 0; NupLit_Htm1(ksv==2,1,3)= 0;
%NupLit_Ltm1(ksv==3,1,1)= 0; NupLit_Ltm1(ksv==3,1,2)= 0; NupLit_Ltm1(ksv==3,1,3)= 0; 
%NupLit_Ltm1(ksv==4,1,1)= 0; NupLit_Ltm1(ksv==4,1,2)= 0; NupLit_Ltm1(ksv==4,1,3)= 0; 
%NupLit_Htm1(ksv==5,1,1)= 0; NupLit_Htm1(ksv==5,1,2)= 0; NupLit_Htm1(ksv==5,1,3)= 0; 
%NupLit_Htm1(ksv==6,1,1)= 0; NupLit_Htm1(ksv==6,1,2)= 0; NupLit_Htm1(ksv==6,1,3)= 0; 

%% PARI:         Smoothed average of PAR radiation over 30 days (1); 45 days (2); (NO REFERENCE)
%               (3) 10 days average of the difference of (2) and (1) [W/m2] 
% categories     [fir_high    Crops_WW  Crops_WB   Crops_S  Crops_R     grass_A   grass_B    shrub    BLever_high  BLdec_low  BLdec_high NoVeg]  
PARI_class_H =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];
PARI_class_L =   [0           0         0          0        0           0         0          0        0            0          0          NaN  ];

PARI_Htm1(ij_vect(k),1:length(valid),1) = PARI_class_H(valid);
PARI_Ltm1(ij_vect(k),1:length(valid),1) = PARI_class_L(valid);

PARI_Htm1(ij_vect(k),1:length(valid),1) = PARI_class_H(valid);
PARI_Ltm1(ij_vect(k),1:length(valid),1) = PARI_class_L(valid);

PARI_Htm1(ij_vect(k),1:length(valid),1) = PARI_class_H(valid);
PARI_Ltm1(ij_vect(k),1:length(valid),1) = PARI_class_L(valid);

%PARI_Htm1(ksv==1,1,1)=0; PARI_Htm1(ksv==1,1,2)=0; PARI_Htm1(ksv==1,1,3)=0;
%PARI_Htm1(ksv==2,1,1)=0; PARI_Htm1(ksv==2,1,2)=0; PARI_Htm1(ksv==2,1,3)=0;  
%PARI_Ltm1(ksv==3,1,1)=0; PARI_Ltm1(ksv==3,1,2)=0; PARI_Ltm1(ksv==3,1,3)=0; 
%PARI_Ltm1(ksv==4,1,1)=0; PARI_Ltm1(ksv==4,1,2)=0; PARI_Ltm1(ksv==4,1,3)=0; 
%PARI_Htm1(ksv==5,1,1)=0; PARI_Htm1(ksv==5,1,2)=0; PARI_Htm1(ksv==5,1,3)=0; 
%PARI_Htm1(ksv==6,1,1)=0; PARI_Htm1(ksv==6,1,2)=0; PARI_Htm1(ksv==6,1,3)=0; 


%% B:            Carbon Pools [gC/m^2] (n=CP)
%(ksv==1,1,1) 
% first: Type of vegetation (6 classifications)
% Second: iteration (thousands)
% Third: Type of biomass - soil (8 classifications)

B_Htm1(ij_vect(k),:,1)= 310;  % Foliage 
B_Htm1(ij_vect(k),:,2)= 556;  % Liv. Sapwwod
B_Htm1(ij_vect(k),:,3)= 406;  % Carbohydrate reserve
B_Htm1(ij_vect(k),:,4)= 411;  % Fruit and Flowers
B_Htm1(ij_vect(k),:,5)= 7;    % Heartwood/Dead sapwood
B_Htm1(ij_vect(k),:,6)= 0;    % Standing dead foliage
B_Htm1(ij_vect(k),:,7)= 16;   % Aux
B_Htm1(ij_vect(k),:,8)= 0;    % High vegetation

%%%
B_Ltm1(ij_vect(k),:,1)= 0;    
B_Ltm1(ij_vect(k),:,2)= 0;
B_Ltm1(ij_vect(k),:,3)= 303;
B_Ltm1(ij_vect(k),:,4)= 235;
B_Ltm1(ij_vect(k),:,5)= 0.3;  
B_Ltm1(ij_vect(k),:,6)= 0;  
B_Ltm1(ij_vect(k),:,7)= 13;
B_Ltm1(ij_vect(k),:,8)= 0;


%%%  

%{
B_Htm1(ij_vect,1,1)= 0;      
B_Htm1(ij_vect,1,2)= 330;
B_Htm1(ij_vect,1,3)= 178;
B_Htm1(ij_vect,1,4)= 250;
B_Htm1(ij_vect,1,5)= 1;  
B_Htm1(ij_vect,1,6)= 0;
B_Htm1(ij_vect,1,7)= 0.3;
B_Htm1(ij_vect,1,8)= 0;

%%%
B_Ltm1(ij_vect,1,1)= 0;   
B_Ltm1(ij_vect,1,2)= 141;
B_Ltm1(ij_vect,1,3)= 61;
B_Ltm1(ij_vect,1,4)= 101;
B_Ltm1(ij_vect,1,5)= 0;  
B_Ltm1(ij_vect,1,6)= 0;
B_Ltm1(ij_vect,1,7)= 1;
B_Ltm1(ij_vect,1,8)= 0;  
%%%
B_Htm1(ij_vect,3,1)= 160;   
B_Htm1(ij_vect,3,2)= 77;
B_Htm1(ij_vect,3,3)= 456;
B_Htm1(ij_vect,3,4)= 55;
B_Htm1(ij_vect,3,5)= 2;  
B_Htm1(ij_vect,3,6)= 0;
B_Htm1(ij_vect,3,7)= 16;
B_Htm1(ij_vect,3,8)= 0;
%%%  
B_Htm1(ij_vect,3,1)= 0;      
B_Htm1(ij_vect,1,2)= 195;
B_Htm1(ij_vect,1,3)= 352;
B_Htm1(ij_vect,1,4)= 154;
B_Htm1(ij_vect,1,5)= 2;  
B_Htm1(ij_vect,1,6)= 0;
B_Htm1(ij_vect,1,7)= 1;
B_Htm1(ij_vect,1,8)= 0;

%}

end
end


%Save data to mat file
save(out);

