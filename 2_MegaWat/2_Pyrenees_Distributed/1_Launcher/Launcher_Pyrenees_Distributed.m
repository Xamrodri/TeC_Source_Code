%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
% Created on Nov 25, 2024
% Author: MAXIMILIANO RODRIGUEZ
% Code originally from: ACHILLE JOUBERTON
% Area of Study: Pyrenees
% Region: Cinca Subcatchment
% Code explanation: This code launches TC model.
% Output variables are managed by the file OUTPUT_MANAGER_DIST_LABEL.
% Some code and names still referes to previous versions of the code.
%==========================================================================

%% INITIALIZING
close all
clc
% clear all   % Don't clear all such that it can be run on the HPC cluster without any issues 
% delete(gcp('nocreate'))

%% DIRECTORIES
study_name = '2_Pyrenees_distributed';
%local_machine = 'WSL28243'; % put your computer name here (to differentiate it with the HPC cluster)

%% Initial inputs by the user
% Select cluster or personal computer
ISTA = 1; % Options: 1: yes for local computer, 2: Yes for cluster
% Select catchment
s = 2;  % This launcher can be used for different catchment. Here in this case study, only sitenumber = 2 works.

%% 
if ISTA == 1 % << Local computer of MR >>
path_output = 'M:/19_ISTA/1_TC/3_Model_Source/2_MegaWat/2_Pyrenees_Distributed/5_Output_files/1_Model_outputs';
folder_path = 'M:/19_ISTA/1_TC/3_Model_Source/2_MegaWat'; % Put here the path of where you downloaded the repository
else % << Or Cluster - MR login >>
path_output = '/home/mrodrigu/1_TC/2_Outputs/2_Model_outputs';  %'/home/jouberto/TC_outputs/Kyzylsu/Distributed';
folder_path = '/home/mrodrigu/1_TC/1_Model'; % Put here the path of where you downloaded the repository. '/home/jouberto/TeC_Source_Code'
end

%% Site specifications

restart = 0; % Set to 1 to continue a un-completed T&C run (only works if you used OUTPUT_MANAGER_DIST_LABEL)
fn_restart = '250924_1999_2023_PG00_Tmod0_2000mmSWEcap_a14c145';

%{
% on the HPC, 's' comes from SLURM script
machine='WSL28243'; %getenv('WSL28243');
if strcmp(machine,local_machine) %%% << Achille's laptop >>
    s=sitenumber; 
end
%}

SITEs = ["Langtang", "Cinca_Mid", "Parlung24K","Rolwaling","Parlung4","Mugagangqiong"]; % All of my study catchments
FORCING = 'Forcing_(Chelsa)';
Lats = [28.2108, 42.4269, 29.773,27.836191,29.233504,32.235];
Lons = [85.5695,0.1776,95.699, 86.531051, 96.923511,87.486];
DeltaGMTs=[5.75,5,8,5.45,8,8];
Zbass=[3862,3579,3800,5449,4600,5850]; %in this setup: elevation of the reference pressure logger

% dtm file name
dtm_files = ["dtm_dtSMB_Langtang_100m.mat", "dtm_Cinca_Mid_250m.mat", "dtm_24K_100m.mat","dtm_Rolwaling_200m.mat","dtm_Parlung4_100m.mat","dtm_Mugagangqiong_100m.mat"];

% Catchment outlet POI name
outlet_names = ["qqq","Cinca_Mid","ccc","Streamgauge_lake","Further_outlet_AWS4600","Outlet_gaging_station"];

% Parameters
Tmods = [0,0,0,0,0,0]; % temperature modification above clean ice [°C];
Pmods = [0,0,0,0.0,0,0,0]; % factor Pmod preciptation modification, e.g. 0.3 means 30% more precipitation at highest elevation
Z_mins = [3000,500,3000,4000,4400,4400]; % lowest elevation for linear precipitation modification (min factor -> 0)
Z_maxs = [6000,5000,6000,6000,6000,6000]; % highest elevation for linear precipitation modification (max factor -> Pmod)
Tmaxs = [2, 2, 1.8, 2, 1.8, 1.8]; % Maximum air temperature for precipitation phase scheme (2-dual thresholds)
Tmins = [0, 0, -0.5, 0, -0.5, -0.5]; % Minimum air temperature for precipitation phase scheme (2-dual thresholds)
 

%% Simulation period
% Starting point of the simulation
x1s = ["01-Oct-2018 00:00:00", "01-Sep-2021 00:00:00", "01-Oct-2018 00:00:00","01-Oct-2018 00:00:00","01-Jan-2015 00:00:00","01-Jan-2015 00:00:00"];
% Ending point of the simulation
x2s = ["30-Sep-2020 23:00:00", "02-Sep-2021 23:00:00", "28-Apr-2022 23:00:00","01-Dec-2019 00:00:00","31-Jan-2015 23:00:00","31-Jan-2015 23:00:00"];

if exist('date_start','var')
    x1 = datetime(date_start);
    x2 = datetime(date_end);
  else 
    x1=datetime(x1s(s));
    x2=datetime(x2s(s));
end 

%% OUTPUT MANAGER : decide which aggregated results to output
% Veg: Vegetation
% LC: Land Cover class

           % Catch_av  Catch_std   Veg_avg  Veg_std  LC_avg   LC_std  Maps  POI
output_manag = [1,         0,         1,        0,      0,       0,     1    1];

%% GLACIER AND AVALANCHE PARAMETERS
% Switch for glacier dynamics (keep to 0 for this case study)
Idyns = [0,0,0,0,0,0];

% switch on/off avalanching
Avals = [1,1,1,1,1,1];  % 1 to turn on avalanching, 0 to turn off avalanching
a_avals = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]; % avalanche parameters a (Bernhart & Schulz 2010)
C_avals = [145, 145, 145, 145, 145, 145]; % avalanche parameters C 

%% SNOW INITIAL CONDITION
% Initial snow depth
% Set at 2000 masl initially.
diff_IniSND = 1;
fn_IniSnowDepth = 'Cinca_Init_Snow_Depth_virtual.mat';
fn_IniSnowAlbedo = 'Cinca_Init_Snow_Albedo_virtual.mat';

%% Precipitation phase partitionning

%1 = 2-threshold, 2 = Ding 2017, 3 = single-threshold, 4 = Pomeroy 2013, 5
%= Wang 2019, 6 = Jennings 2018

parameterize_phase.OPT_Pr_Part = 2; % Choice of the precipitation phase scheme
parameterize_phase.Tmax = Tmaxs(s); % Upper air temperature for dual temperature threshold
parameterize_phase.Tmin = Tmins(s); % Lower air temperature for dual temperature threshold
parameterize_phase.Tconst = 2; % Air temperature for constant thresholds

parameterize_phase_labels = {'2-Ta','Ding','1-Ta','Pomeroy','Wang','Jennings'};
parameterize_phase_label = parameterize_phase_labels(parameterize_phase.OPT_Pr_Part);

% Skin layer thickness for the 2-layer snowpack module
hSTL = 0.003; %m

% Choice of the snow albedo scheme 
Albsno_method = 5; % 3 doesn't work, 4 is Brock 2000, 5 is Ding 2017

%% Simulation comment
sim_comment = "Distributed_Ebro";
if ~exist('sim_comment','var')
 sim_comment = "Distributed_Pyrenees";
end 

%% Simulation name
simnm = strcat('Region_(',sim_comment,')_SubBasin_(',outlet_names(s),')_', ...
    FORCING, "_RunningDate_(", datestr(datetime("now"),'ddmmyy_HHMM'),")_RunningPeriod_(",num2str(year(x1)),"_", ...
    num2str(year(x2)) ,")");

%% Select which forcing variable is bias-corrected or not

% vars: define which reanalysis variables to load as monthly spatial input 
% vars_DS: 0 downscaled and bias-corrected, 1 only downscaled     
vars = ["FF", "LWIN", "PP", "RH", "TA", "SWPART", "PRESS"];
vars_DS = [
     0,0,0,0,0,0,0;
     0,0,0,0,0,0,0;
     0,0,0,0,0,0,0;
     0,1,1,0,0,1,0
     0,0,1,0,0,0,0
     0,0,0,0,0,1,0];

% Solar variables 
t_befs = [1,0.5,0.75,0.75,1,1]; % Integration interval for solar variables (hour)
t_afts = [0,0.5,0.25,0.25,0,0]; % Integration interval for solar variables (hour)

%% Prepare launcher variable based on chosen study site:
SITE = char(SITEs(s)); Lat = Lats(s); Lon = Lons(s);
DeltaGMT=DeltaGMTs(s); t_bef=t_befs(s); t_aft=t_afts(s); 
dtm_file = char(dtm_files(s)); Tmod = Tmods(s); 
outlet_name = char(outlet_names(s));
Pmod = Pmods(s); Z_min = Z_mins(s); Z_max = Z_maxs(s);
vars_load = vars_DS(s,:); 
Zbas=Zbass(s); 
Idyn = Idyns(s);
Aval = Avals(s); a_aval = a_avals(s); C_aval = C_avals(s);
TITLE_SAVE = SITE;
resol = str2num(string(extractBetween(dtm_file,[SITE '_'],'m.mat'))); % Spatial resolution

%% Diplay setting of the incoming T&C model runs:
disp(['Site selected: ' TITLE_SAVE])
disp(['Simulation period: ' datestr(x1,'dd-mmm-yyyy HH:MM') ' to ' datestr(x2,'dd-mmm-yyyy HH:MM')])
disp(['Precipitation phase scheme: ' parameterize_phase_label{:}])

if Idyn == 0; disp('Ice dynamics: off'); else; disp('Ice dynamics: on'); end
if Aval == 0; disp('Avalanching: off'); else; disp('Avalanching: on'); end

%% PATHS
%==========================================================================
% Attach folders, launch parallel pool, generate paths
%==========================================================================

addpath(genpath([folder_path, '/' study_name '/2_Inputs'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path, '/' study_name '/3_Forcing'])); % Where is located the meteorological forcing and Shading matrix 
addpath(genpath([folder_path, '/5_Common_inputs'])); % Where is located common inputs of the distributed and point scale models 
addpath(genpath([folder_path, '/' study_name '/2_Inputs/1_Input_Carbon'])); % Add path to Ca_Data    
addpath(genpath([folder_path, '/1_Functions'])); % Where are distributed model set-up files (needed ? yes to load dtm)
outlocation = [path_output '/' char(simnm)];


% Create output directory if it doesn't exist already
if ~exist(outlocation, 'dir'); mkdir(outlocation); addpath(genpath(outlocation)); end

%Deactivate hyperthreading
N=1;
LASTN=maxNumCompThreads(N);
%%%

%% If we are NOT restarting. Normal case.
if restart ~=1

checks = [outlocation, '/checks'];
mkdir(outlocation); %Create main folder
mkdir(strcat(outlocation,'/Spatial_data')); %Create main subfolder where to store spatial results
mkdir(strcat(outlocation,'/Cell_data_final')); %Create main subfolder where to store results at a cell.
mkdir(strcat(outlocation,'/Spatial_data_intermediate')); %Create main subfolder where to store spatial results at specific intervals.

%mkdir(checks);
addpath(genpath(outlocation));
end 


%% INPUTS
%==========================================================================
% Load geodata, time handling, carbon data
%==========================================================================
load(dtm_files(s)) %This file is in 2_nputs

if Idyn == 0
    GLH(GLH>0) = GLH(GLH>0) +400;  % Only use when ice dynamics are off, to avoid glacier disappearance
end

% Load different initial snow depth and snow albedo from the one in the dtm_file
path_fn_IniSND = [folder_path '/2_Pyrenees_distributed/2_Inputs/2_Ini_SnowDepth/' fn_IniSnowDepth];
path_fn_IniSNOWALB = [folder_path '/2_Pyrenees_distributed/2_Inputs/2_Ini_SnowDepth/' fn_IniSnowAlbedo];

if (diff_IniSND == 1) && exist(path_fn_IniSND,'file')>0
    load(path_fn_IniSND)
end 

if (diff_IniSND == 1) && exist(path_fn_IniSNOWALB,'file')>0
    load(path_fn_IniSNOWALB)
end 

% PATCH, check pre-processing!
GLA_ID(MASK==0) = NaN;
GLA_ID(GLH==0)=NaN;

% Compute bedrock DEM for ice flow
DTM_Bedrock = DTM_orig-GLH; 

load([folder_path, '/2_Pyrenees_distributed/2_Inputs/1_Input_Carbon/Ca_Data.mat']);

clear Xvs Yvs
%%%

Date = x1:hours(1):x2;
N_time_step=length(Date);
%%%
d1 = find(abs(Date_CO2-datenum(Date(1)))<1/36);d2 = find(abs(Date_CO2-datenum(Date(end)))<1/36);
Ca=Ca(d1:d2);
clear d1 d2 Date_CO2
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%
Nd_time_step = ceil(N_time_step/24)+1;

%% find last timestep per month for intermediate output
tstore = find(diff(month(Date)) ~= 0); % general spatial output
tstore_alb = find((day(Date) == 1 & hour(Date)==12) | (day(Date) == 15 & hour(Date)==12)) ; %half-monthly albedo output
tstore_snow = find(hour(Date)==12); %daily snow map output

t1_reinit=0; %%only placeholder for re-initialisation runs
plot_soil = 0; %plot soil property maps?


%% GENERAL PARAMETER
%==========================================================================
%
%==========================================================================

%% Time 
dt=3600; %% [s]
dth=1; %% [h]
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE
%%%
L_day=zeros(length(Datam),1);
for j=2:24:length(Datam)
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')
%%%
a_dis=NaN; pow_dis=NaN;
clear a0 gam1 pow0 k2 DTii

%% Other parameters
rho_g = 0.35; %%% Spatial Albedo
cc_max = 1; %% Number of vegetation ??
ms_max = 10; %% Number of soil layers
md_max = 10; %% % Number of debris layers

%Set zatm
zatm_surface = 2.0; %% Reference Height single value
% zatm_hourly_on=0; %Switch for using variable instrument heights, 1=on
% if zatm_hourly_on==1
%     zatm_hourly=((rand(N_time_step,1)*10)+250)/100; % Instrument height timeseries (m), this is turned into zatm in Main_Frame [Dummy variable made here for testing]
% else
%     zatm_surface = 2.0; %% Reference Height single value
% end

%% Albedo vs Elevation
% Load mean glacier albedo eleavtion profile, development stage

fn_alb_elev = [SITE '_Albedo_vs_elev.mat'];

if exist(fn_alb_elev,'file')>0
    load([SITE '_Albedo_vs_elev'])

    Afirn = DTM.*0 + 0.28;
    dem_inc = DTM <= Alb_el_tt{1,2};
    Afirn(dem_inc) = Alb_el_tt{1,1};

        for ii = 1:size(Alb_el_tt,1)
              if ii == size(Alb_el_tt,1)
                 dem_inc = DTM >= Alb_el_tt{ii,2};
                 Afirn(dem_inc) = Alb_el_tt{ii,1};
              else
                 dem_inc = DTM > Alb_el_tt{ii,2} & DTM < Alb_el_tt{ii+1,2};
                 Afirn(dem_inc) = 0.5*(Alb_el_tt{ii,1} + Alb_el_tt{ii+1,1});
              end 
        end 
    Afirn(Afirn > 0.4) = 0.4; %Limit bare ice albedo to 0.4, as above it's firn. 
    disp('Measured bare-ice abledo is used')
else 
    Afirn = DTM.*0 + 0.28;
    disp('Constant bare-ice albedo of 0.28 is used')
end

%%
Pmod_S = MASK;
rate = Pmod/(Z_max-Z_min);
Pmod_S(DTM>Z_min) = 1+rate.*(DTM(DTM>Z_min)-Z_min);


%% TOPOGRAPHIC PARAMETER 
%==========================================================================
%
%==========================================================================
%%%   Matrix [m_cell x n_cell]
%DTM ;T_flow ; cellsize; xllcorner ; yllcorner ; SN ; outlet ; Aacc

GLH=GLH.*(GLA_MAP2>0);
clear m n %MASK
[m_cell,n_cell]=size(DTM);
num_cell=numel(DTM);
x_cell=xllcorner:cellsize:(xllcorner+cellsize*(n_cell-1));
y_cell=yllcorner:cellsize:(yllcorner+cellsize*(m_cell-1));
% MASK=ones(m_cell,n_cell); MASK(isnan(DTM))=0;
MASKn=reshape(MASK,num_cell,1);
Kinde = find(MASK==1);
%"Xoutlet" & "Youtlet": outlet point, predefined in dtm_XXX.mat-file

%%% Slo_top [Fraction] %%% Aspect [rad] from N
[Slo_top,Aspect]=Slope_Aspect_indexes(DTM_orig,cellsize,'mste');
Aspect(isnan(Aspect))=0;
Slo_top(Slo_top<0.001)=0.001;  %to avoid flow routing issues
%%%
Asur=(1./cos(atan(Slo_top))); %% Effective Area / Projected Area
Asur=reshape(Asur,num_cell,1);
aTop= 1000*ones(m_cell,n_cell)*(cellsize^2)/cellsize; %%[mm] Area/Contour length ratio
Ared=ones(num_cell,1);
%%%check: >> imagesc(DTM),axis xy,axis equal, colorbar()
%"Xout" & "Yout": tracking points, predefined in dtm_XXX.mat-file

%%% Flow Boundary Condition
Slo_top(Youtlet,Xoutlet)=0.05;
npoint = length(Xout);
Area= (cellsize^2)*sum(sum(MASK)); %% Projected area [m^2]

%%% Flow potential
T_pot=cell(1,ms_max);
for jk=1:ms_max
    T_pot{jk}= T_flow;
end
%%%
WC = cellsize*ones(m_cell,n_cell); %% [m]  Width channel
WC(SN==1)=0.0018*sqrt((cellsize^2)*Aacc(SN==1)); %% [m]
WC=WC.*MASK;
SN(isnan(SN))=0; %% [Stream Identifier]
SNn=reshape(SN,num_cell,1);
NMAN_C=SN*0.040; NMAN_H=0.1; %%[s/(m^1/3)] manning coefficient
MRough = 0.01*(1-SN); %%[m] Microroughness
NMAN_C=NMAN_C.*MASK;
NMAN_H=NMAN_H.*MASK;
MRough=MRough.*MASK;
%%%
Kres_Rock =8760; %%[h] Bedrock aquifer constant
SPRINGn =SNn; %% Spring Location


%% VEGETATION PARAMETERS LOOK-UP TABLE 
%==========================================================================
%
%==========================================================================
%ksv=1*ones(num_cell,1);
ksv=reshape(VEG_CODE,num_cell,1);
%%% 1 Fir (evergr.)
%%% 2 Larch (decid.)
%%% 3 Grass C3
%%% 4 Shrub (decid.)
%%% 5 Broadleaf evergreen
%%% 6 Broadleaf deciduous
%%% 7 Rock
ksv(ksv==8)=7; %%% 8 Ice = 7 Rock

%Ccrown_OUT and EVcode are used in the OUTPUT_MANAGER_DIST_LABEL.m
Ccrown_OUT =[ 1 ; 1; 1 ; 0.9 ; 1; 1; 0];  %% Ccrown fraction for Plant Functional Types (PFT)
EVcode = [ 1 2 3 4 5 6 7 ];             %% code of each PFTs

%{
%Mike
Ccrown_OUT = [ 1   0 ; 0.5  0.5 ; 1   0 ; 1   0];  %% Ccrown fraction for PFT
EVcode     = [ 1      2           3       4    ];  %% Code of each PFT
cc_max = size(Ccrown_OUT,2);                       %% Number of vegetation types per cell (max.)
%}

%% SPATIAL INDICES PER LAND COVER CLASS

Kinde = find(MASK==1);  %%% basin index
idx_Veg1 = find(VEG_CODE == 1 & MASK == 1); %%% veg 1 index (Fir)
idx_Veg2 = find(VEG_CODE == 2 & MASK == 1); %%% veg 2 index (Larch)
idx_Veg3 = find(VEG_CODE == 3 & MASK == 1); %%% veg 3 index (Grass)
idx_Veg4 = find(VEG_CODE == 4 & MASK == 1); %%% veg 4 index (Shrub)
idx_Veg5 = find(VEG_CODE == 5 & MASK == 1); %%% veg 5 index (Broadleaf evergreen)
idx_Veg6 = find(VEG_CODE == 6 & MASK == 1); %%% veg 6 index (Broadleaf deciduous)
idx_Rock = find(VEG_CODE == 7 & MASK == 1); %%% rock index
idx_Ice = find(VEG_CODE == 8 & MASK == 1);  %%% ice index
idx_Cleanice = find(VEG_CODE == 8 & MASK == 1 & DEB_MAP == 0);%%% clean-ice index
idx_Debice = find(VEG_CODE == 8 & MASK == 1 & DEB_MAP > 0);%%% debris-covered ice index
LCinde = {idx_Veg1,idx_Veg2,idx_Veg3,idx_Veg4,idx_Veg5,idx_Veg6,...SOIL_PARAMETERS
    idx_Rock,idx_Ice,idx_Cleanice,idx_Debice};


%% SOIL PARAMETER 
%==========================================================================
%
%==========================================================================

Pss = [800]; Pwp = [3500]; %% [kPa]
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
[Osat,L,Pe,Ks,O33]=Soil_parameters_spatial(PSAN/100,PCLA/100,PORG/1000);
[Ofc,Oss,Owp,Ohy]=Soil_parametersII_spatial(Osat,L,Pe,Ks,O33,Kfc,Pss,Pwp,Phy);
%plot soil properties

%% Storing in the normal case.
if restart ~=1
    save(strcat(outlocation,"/Spatial_data/", "Ohy_1.mat"), 'Ohy')
    save(strcat(outlocation,"/Spatial_data/", "Osat_1.mat"), 'Osat')
end 

if plot_soil==1
    figure('Name', 'PSAN');imagesc(PSAN);title('PSAN');colorbar;saveas(gcf,[checks, '/PSAN.png'])
    figure('Name', 'PORG');imagesc(PORG);title('PORG');colorbar;saveas(gcf,[checks, '/PORG.png'])
    figure('Name', 'PCLA');imagesc(PCLA);title('PCLA');colorbar;saveas(gcf,[checks, '/PCLA.png'])
    figure('Name', 'PSOIL_SUM');imagesc(PCLA+PORG+PSAN);title('PSOIL_SUM');colorbar;saveas(gcf,[checks, '/PSOIL_SUM.png'])
    figure('Name', 'Osat');imagesc(Osat);title('Osat');colorbar;saveas(gcf,[checks, '/Osat.png'])
    figure('Name', 'Ohy');imagesc(Ohy);title('Ohy');colorbar;saveas(gcf,[checks, '/Ohy.png'])
    figure('Name', 'Ks');imagesc(Ks);title('Ks');colorbar;saveas(gcf,[checks, '/Ks.png'])
end
%%
clear Pss Pwp Kfc Phy L Pe Ks O33 Ofc Oss Owp
Osat_OUT =   Osat.*MASK; clear Osat
Osat_OUT=reshape(Osat_OUT,num_cell,1);
Ohy_OUT =   Ohy.*MASK; clear Ohy
Ohy_OUT=reshape(Ohy_OUT,num_cell,1);


%Mdep1 = 50; Mdep2= 100; Mdep3=250;
%%%
PSAN=reshape(PSAN/100,num_cell,1);
PCLA=reshape(PCLA/100,num_cell,1);
PORG=reshape(PORG/1000,num_cell,1);
%Zs_OUT=reshape(SOIL_TH*(2/3)*10,num_cell,1); %%[mm]  %%orig: 
Zs_OUT=800*ones(num_cell,1);


%% SOLAR PARAMETER 
%==========================================================================
%
%==========================================================================
% Computation Horizon Angle
%[HZ,Zasp] = Horizon_angle_polar(DTM,cellsize);
[HZ,Zasp] = Horizon_Angle(DTM_orig,cellsize);
%%% HZ Horizon angle array [angular degree]
%%% Z Azimuth directions  [angular degree] from N
%%% Sky View Factor and Terrain Configuration Factor
[SvF,Ct] = Sky_View_Factor(DTM_orig,atan(Slo_top)*180/pi,Aspect,HZ,Zasp);


%% TOTAL INITIAL CONDITION 

%SNOW ALBEDO
% Check if initial snow albedo is given
if ~exist('SNOWALB','var')
    SNOWALB = SNOWD;
    SNOWALB(SNOWD>0) = 0.6;
end 

%% If I am NOT restarting. So, in the normal case.
if restart ~=1
out = [outlocation,'/Spatial_data/', 'INITIAL_CONDITIONS_', SITE, '.mat'];
INIT_COND_v2(num_cell,m_cell,n_cell,...
    cc_max,ms_max,md_max,...
    MASKn,GLH,Slo_top,ksv,Ca,SNOWD,SNOWALB,out);
load(out);
end

Fstep=strcat('/Spatial_data/Final_reached_step_',SITE);
ij=0;
tday=0;

%% NUMERICAL METHODS OPTIONS
%==========================================================================
%
%==========================================================================
OPT_ALLOME=                 0;
Opt_CR=                     optimset('TolFun',1);%,'UseParallel','always');
OPT_EnvLimitGrowth=         0;
OPT_FR_SOIL=                1;  % Option for freezing soil
OPT_HEAD=                   0;
OPT_PH=                     odeset('AbsTol',0.01);
OPT_PlantHydr=              0;
OPT_Pr=                     2; % Option for precipitation phase calculation, 1: Wigmosta et.al 1994, 2: Ding et.al 2014
OPT_SM=                     odeset('AbsTol',0.05,'MaxStep',dth);
OPT_SoilBiogeochemistry=    0;
OPT_SoilTemp=               1;
Opt_ST=                     optimset('TolFun',0.1);%,'UseParallel','always');
Opt_ST2 =                   optimset('TolFun',0.1,'Display','off');
OPT_STh=                    odeset('AbsTol',5e+3);
OPT_VCA=                    0;
OPT_VD=                     odeset('AbsTol',0.05);
OPT_VegSnow=                1;
OPT_min_SPD =               0.006;            %% [m] minimum snow pack depth to have a multilayer snow 

tic ;
%profile on
%bau = waitbar(0,'Waiting...');

%% Restart condition
fts = 2; 
tinp = 0;
root = 'M:';

if restart == 1  % Restart simulation from a previously saved step
    outlocation = [root,'/TC_outputs/',SITE,'/Distributed/',fn_restart,'/'];
    load([outlocation '/Spatial_data/Final_reached_step_' SITE '.mat'])
    fts = t; 
    t1_reinit = t;
    %Find the right value of tinp
    Datam_reinit = Datam(t1_reinit-1,:);
    tinp = (Datam_reinit(3)-1)*24 + Datam_reinit(4);
end

%% Iterating on time
t=3;
for t=fts:N_time_step
    tinp = tinp+1;
    
    %waitbar(t/N_time_step,bau)
   disp('Iter:'); disp(t);

    %pdind = [max(1,tinp-24):tinp-1]; % previous day indexes MAYBE NEEDS TO BE REPLACED WITH SOMETHING ELSE BECAUSE OF MONTHLY INPUTS
    Datam_S=Datam(t-1,:);
    if (Datam_S(4)==1)
        tday = tday+1;
    end
    
    if (Datam_S(3)==1 && Datam_S(4)==0)
        tinp=1;
    end
   
    [jDay]= JULIAN_DAY(Datam_S);
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day]= SetSunVariables(Datam_S,DeltaGMT,Lon,Lat,t_bef,t_aft);
    [ShF] = Shadow_Effect(DTM,h_S,zeta_S,HZ,Zasp);

    
    %% DISTRIBUTED FORCINGS
    %======================================================================
    % Forcing from stations and distributed across the catchment.
    % load spatial input data on first day and hour of month 
    %======================================================================      

    %% Stations
    Elev = [3371, 3900]; % Elevation of the two meteorological stations

    Station_1 = load('Pluvio.mat'); Station_1 = Station_1.forcing_all;
    Station_2 = load('AWS_3900m.mat'); Station_2 = Station_2.forcing_all;

    t_bef=t_befs(s); t_aft=t_afts(s); % otherwise problems when loading SWPART
        
    %% Selection of time period
    % Find the indices in the meteorological time-series corresponding 
    % to the timestamp of the simulation timestep.
    ind_meteo = find(datetime(Datam_S(1), Datam_S(2), Datam_S(3), Datam_S(4),0,0) == Station_1.Time);
    
    %% Temperature
    % Air temperature from stations and selected period
    Ta_1 = Station_1.TA(ind_meteo);
    Ta_2 = Station_2.TA(ind_meteo);
    
    % Regression using two stations
    AUX = polyfit(Elev,[Ta_1 Ta_2],1);
    LR = AUX(1); clear AUX; % air tempetature lapse-rate
    
    % Using stations to obtain distributed temperature along the catchment
    Ta_S = DTM; Ta_S = Ta_1 + LR.*(DTM-Elev(1)); % Distributed air temperature, extrapolated from pluviometer station
    Ta_S(GLH>0 & DEB_MAP == 0) = Ta_S(GLH > 0 & DEB_MAP == 0)+Tmod; % Temperature modification above clean-ice via Tmod
    Ta_S = reshape(Ta_S,num_cell,1);
    Ta_S(MASKn==0) = NaN;

    %%% define Ta_day array (Ta-maps of last 24 hours), needed for albedo scheme
    %%% parameterization

    if t == fts
        Ta_day = Ta_S;
    elseif t > fts && t <= fts+24
        Ta_day = [Ta_day,Ta_S];
    elseif t > fts+24
        Ta_day = [Ta_day(:,2:24),Ta_S];
    end
    
    %% Longwave radiation (for this case study, it's constant in space)
    N_1 = Station_1.LWIN(ind_meteo);
    N_S = DTM.*0 + N_1;
    N_S = reshape(N_S,num_cell,1);
    N_S(MASKn==0) = NaN;

    %% Relative humidity (ERA5)
    U_1 = Station_1.RH(ind_meteo)./100;
    U_S = DTM.*0 + U_1;
    U_S = reshape(U_S,num_cell,1);
    U_S(U_S<0)=0; U_S(U_S>1)=1;
    U_S(MASKn==0) = NaN;

    %% Wind speed
    Ws_1 = Station_1.FF(ind_meteo);
    Ws_S = DTM.*0 + Ws_1;
    Ws_S = reshape(Ws_S,num_cell,1);
    Ws_S(isnan(Ws_S)) = 0.01;
    Ws_S(Ws_S < 0.01) = 0.01;
    Ws_S = Ws_S.*MASKn;

    %% Precipitation
    Pr_1 = Station_1.PP(ind_meteo);
    Pr_S = DTM.*0 + Pr_1;
    Pr_S = reshape(Pr_S,num_cell,1);
    Pr_S(isnan(Pr_S)) = 0;
    Pr_S(Pr_S<0.001)=0;

    %% Air pressure
    Pre_1 = Station_1.PRESS(ind_meteo);
    Pre_S = DTM.*0 + Pre_1;
    Pre_S = reshape(Pre_S,num_cell,1);
    
    %needed, if terrain effects have not been considered during pre-processing
    cos_fst = cos(atan(Slo_top))*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;

    %% Radiation (SW)
    % ERA5 downscaled/bias-corrected shortwave radiation (6 components)
    if sin(h_S) <= 0.10 %%[5.73�] Numerical problems
        SAB1_S  =  0*MASKn; %First band direct radiation
        SAB2_S  =  0*MASKn; %Second band direct radiation
        SAD1_S  =  0*MASKn; %First band diffuse radiation
        SAD2_S  =  0*MASKn; %Second band diffuse radiation
        PARB_S  =  0*MASKn; %PAR radiation direct. PAR: Photosynthetically Active Radiation.
        PARD_S  =  0*MASKn; %PAR radiation diffuse
    else

        SAD1_1 = Station_1.SAD1(ind_meteo); SAD1_S = DTM.*0 + SAD1_1;
        SAD2_1 = Station_1.SAD2(ind_meteo); SAD2_S = DTM.*0 + SAD2_1;
        SAB1_1 = Station_1.SAB1(ind_meteo); SAB1_S = DTM.*0 + SAB1_1;
        SAB2_1 = Station_1.SAB2(ind_meteo); SAB2_S = DTM.*0 + SAB2_1;
        PARD_1 = Station_1.PARD(ind_meteo); PARD_S = DTM.*0 + PARD_1;
        PARB_1 = Station_1.PARB(ind_meteo); PARB_S = DTM.*0 + PARB_1;


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SAD1_S = SAD1_S.*SvF + Ct.*rho_g.*(SAD1_S./sin(h_S).*cos_fst + (1-SvF).*SAD1_S);
        SAD2_S = SAD2_S.*SvF + Ct.*rho_g.*(SAD2_S./sin(h_S).*cos_fst + (1-SvF).*SAD2_S);
        PARD_S = PARD_S.*SvF + Ct.*rho_g.*(PARD_S./sin(h_S).*cos_fst + (1-SvF).*PARD_S);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SAB1_S =(SAB1_S./sin(h_S)).*cos_fst.*ShF;
        SAB2_S =(SAB2_S./sin(h_S)).*cos_fst.*ShF;
        PARB_S = (PARB_S./sin(h_S)).*cos_fst.*ShF;        
	    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %grad1= 0.0049; % gradient [[W/m^2]/m]
        %SAB1_S= (SAB1_S + grad1*(DTM-qta_Sta));
        %grad1= 0.0080; % gradient  [[W/m^2]/m]
        %SAB2_S= (SAB2_S + grad1*(DTM-qta_Sta));
        %grad1=  0.0047; % gradient [[W/m^2]/m]
        %PARB_S= (PARB_S + grad1*(DTM-qta_Sta));
        %%%%%%%%%%%%%%%%%%temporary, remove after next version of input SW
        SAB1_S(SAB1_S<0)=0;
        SAB2_S(SAB2_S<0)=0;
        PARB_S(PARB_S<0)=0;
        PARD_S(PARD_S<0)=0;
        SAD1_S(SAD1_S<0)=0;
        SAD2_S(SAD2_S<0)=0;
        SAB1_S(isnan(SAB1_S)) = 0;
        SAB2_S(isnan(SAB2_S)) = 0;
        SAD1_S(isnan(SAD1_S)) = 0;
        SAD2_S(isnan(SAD2_S)) = 0;
        PARB_S(isnan(PARB_S)) = 0;
        PARD_S(isnan(PARD_S)) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SAB1_S  =  reshape(SAB1_S,num_cell,1);
        SAB2_S  =  reshape(SAB2_S,num_cell,1);
        SAD1_S  =  reshape(SAD1_S,num_cell,1);
        SAD2_S  =  reshape(SAD2_S,num_cell,1);
        PARB_S  =  reshape(PARB_S,num_cell,1);
        PARD_S  =  reshape(PARD_S,num_cell,1);

    end
       
    %% Carbon
    Ca_S = Ca(t)*MASKn;
    IrD_S =  MASK*0; IrD_S=reshape(IrD_S,num_cell,1);
    Salt_S =  MASK*0; Salt_S=reshape(Salt_S,num_cell,1);
    %N_S = N(t)*MASKn;
    
    %% Vapor pressure
    % esat/ea/Ds/Tdew maps
    a=17.27; b=237.3;
    esat_S=611*exp(a*Ta_S./(b+Ta_S)); %Vapour pressure at saturation (Pa)
    ea_S=U_S.*esat_S;                 %Vapour pressure (Pa)
    Ds_S= esat_S - ea_S;              %Vapor Pressure Deficit (Pa)
    Ds_S(Ds_S<0)=0; 
    xr=a*Ta_S./(b+Ta_S)+log10(U_S);
    Tdew_S=b*xr./(a-xr);               %Presumed dewpoint temperature (�C)
    clear a b xr;
    esat_S=reshape(esat_S,num_cell,1);
    ea_S=reshape(ea_S,num_cell,1);
    Ds_S=reshape(Ds_S,num_cell,1);
    Tdew_S=reshape(Tdew_S,num_cell,1);
      
    
    %% SPATIAL INITIALIZATION VECTOR PREDEFINING
    %======================================================================
    %
    %======================================================================
    
    if t == 2
        %%%%%%%%%%%%%%%%% GENERAL VEGETATION / HYDROLOGY
        alp_soil=	alp_soiltm1;
        b_soil=     b_soiltm1;
        Bam=        Bamtm1;
        Bem=        Bemtm1;
        BLit=       BLittm1;
        Ccrown_t=	Ccrown_t_tm1;
        Cice=       Cicetm1;
        Cicew=      Cicewtm1;
        CK1=        CK1tm1;
        Csno=       Csnotm1;
        Csnow=      Csnowtm1;
        dQ_S=       dQ_Stm1;
        DQ_S=       DQ_Stm1;
        dQVEG=      dQVEGtm1;
        DT_S=       DT_Stm1;
        dw_SNO=     dw_SNOtm1;
        e_sno=      e_snotm1;
        EG=         EGtm1;
        EICE=       EICEtm1;
        EIn_rock=	EIn_rocktm1;
        EIn_urb=	EIn_urbtm1;
        EK=         EKtm1;
        ELitter=	ELittertm1;
        er=         ertm1;
        ESN_In=     ESN_Intm1;
        SSN_In =    SSN_Intm1;
        ESN=        ESNtm1; 
        SSN=        SSNtm1;  
        EWAT=       EWATtm1;
        FROCK=      FROCKtm1;
        f=          ftm1;
        Gfin=       Gfintm1;
        G=          Gtm1;
        H=          Htm1;
        HV=         HVtm1;
        ICE_D=      ICE_Dtm1;
        ICE=        ICEtm1;
        Imelt=      Imelttm1;
        In_H=       In_Htm1;
        In_Litter=	In_Littertm1;
        In_L=       In_Ltm1;
        In_rock=	In_rocktm1;
        In_SWE=     In_SWEtm1;
        In_urb=     In_urbtm1;
        IP_wc=      IP_wctm1;
        Lk_rock=	Lk_rocktm1;
        Lk_wat=     Lk_wattm1;
        Lk=         Lktm1;
        Lpho=       Lphotm1;
        NavlI=      NavlItm1;
        NDVI=       NDVItm1;
        NIce=       NIcetm1;
        NIn_SWE=	NIn_SWEtm1;
        OF=         OFtm1;
        Oice=       Oicetm1;
        OS=         OStm1;
        O=          Otm1;
        POT=        POTtm1;
        Pr_liq=     Pr_liqtm1;
        Pr_sno=     Pr_snotm1;
        Q_channel=	Q_channel;
        Q_exit=     Q_exit;
        q_runon=    q_runon;
        QE=         QEtm1;
        QEV=        QEVtm1;
        Qfm=        Qfmtm1;
        Qi_in=      Qi_in;
        Qi_in_Ro=   Qi_out;
        Qi_out_Rout=Qi_out_Rout;
        Qi_out=     Qi_outtm1;
        Qsub_exit=	Qsub_exit;
        Qv=         Qvtm1;
        r_litter=	r_littertm1;
        r_soil=     r_soiltm1;
        ra=         ratm1;
        Rd=         Rdtm1;
        Rh=         Rhtm1;
        Rn=         Rntm1;
        ros=        rostm1;
        SE_rock=	SE_rocktm1;
        SE_urb=     SE_urbtm1;
        Slo_head=	Slo_head;
        Smelt=      Smelttm1;
        SND=        SNDtm1;
        snow_albedo=snow_albedotm1;
        soil_albedo=soil_albedotm1;
        SP_wc=      SP_wctm1;
        surface_albedo=surface_albedotm1;
        SWE=        SWEtm1;
        SWE_avalanched=   SWE_avalanchedtm1;
        t_sls=      t_slstm1;
        tau_sno=	tau_snotm1;
        Tdamp=      Tdamptm1;
        Tdeb=       Tdebtm1;
        Tdp=        Tdptm1;
        Tdp_snow = Tdp_snowtm1;
        Tice=       Ticetm1;
        Tstm0=      Tstm0;
        Ts=         Tstm1;
        Ts_under = Ts_undertm1;
        TsVEG=      TsVEGtm1;
        U_SWE=      U_SWEtm1;
        Vice=       Vicetm1;
        V=          Vtm1;
        WAT=        WATtm1;
        WIS=        WIStm1;
        WR_IP=      WR_IPtm1;
        WR_SP=      WR_SPtm1;
        Ws_under=	Ws_undertm1;
        ZWT=        ZWTtm1;
        
        %%%%%%%%%%%%%%%%% SPECIFICATIONS FOR HIGH & LOW VEGETATION
        AgeDL_H=	AgeDL_Htm1;
        AgeDL_L=	AgeDL_Ltm1;
        AgeL_H=     AgeL_Htm1;
        AgeL_L=     AgeL_Ltm1;
        AgePl_H=	AgePl_Htm1;
        AgePl_L=	AgePl_Ltm1;
        An_H=       An_Htm1;
        An_L=       An_Ltm1;
        ANPP_H=     ANPP_Htm1;
        ANPP_L=     ANPP_Ltm1;
        B_H=        B_Htm1;
        B_L=        B_Ltm1;
        BA_H=       BA_Htm1;
        BA_L=       BA_Ltm1;
        Bfac_dayH=	Bfac_dayHtm1;
        Bfac_dayL=	Bfac_dayLtm1;
        Bfac_weekH=	Bfac_weekHtm1;
        Bfac_weekL=	Bfac_weekLtm1;
        Ci_shdH=	Citm1_shdH;
        Ci_shdL=	Citm1_shdL;
        Ci_sunH=	Citm1_sunH;
        Ci_sunL=	Citm1_sunL;
        dflo_H=     dflo_Htm1;
        dflo_L=     dflo_Ltm1;
        Dr_H=       Dr_Htm1;
        Dr_L=       Dr_Ltm1;
        e_rel_H=	e_rel_Htm1;
        e_rel_L=	e_rel_Ltm1;
        e_relN_H=	e_relN_Htm1;
        e_relN_L=	e_relN_Ltm1;
        EIn_H=      EIn_Htm1;
        EIn_L=      EIn_Ltm1;
        fapar_H=	fapar_Htm1;
        fapar_L=	fapar_Ltm1;
        FNC_H=      FNC_Htm1;
        FNC_L=      FNC_Ltm1;
        gsr_H=      gsr_Htm1;
        gsr_L=      gsr_Ltm1;
        hc_H=       hc_Htm1;
        hc_L=       hc_Ltm1;
        In_H=       In_Htm1;
        In_L=       In_Ltm1;
        ISOIL_H=	ISOIL_Htm1;
        ISOIL_L=	ISOIL_Ltm1;
        Jsx_H=      Jsx_Htm1;
        Jsx_L=      Jsx_Ltm1;
        Jxl_H=      Jxl_Htm1;
        Jxl_L=      Jxl_Ltm1;
        Kleaf_H=	Kleaf_Htm1;
        Kleaf_L=	Kleaf_Ltm1;
        Kreserve_H=	Kreserve_Htm1;
        Kreserve_L=	Kreserve_Ltm1;
        Kuptake_H=	Kuptake_Htm1;
        Kuptake_L=	Kuptake_Ltm1;
        Kx_H=       Kx_Htm1;
        Kx_L=       Kx_Ltm1;
        LAI_H=      LAI_Htm1;
        LAI_L=      LAI_Ltm1;
        LAIdead_H=	LAIdead_Htm1;
        LAIdead_L=	LAIdead_Ltm1;
        ManIH=      ManIHtm1;
        ManIL=      ManILtm1;
        NBLeaf_H=	NBLeaf_Htm1;
        NBLeaf_L=	NBLeaf_Ltm1;
        NBLI_H=     NBLI_Htm1;
        NBLI_L=     NBLI_Ltm1;
        NPP_H=      NPP_Htm1;
        NPP_L=      NPP_Ltm1;
        NPPI_H=     NPPI_Htm1;
        NPPI_L=     NPPI_Ltm1;
        Nreserve_H=	Nreserve_Htm1;
        Nreserve_L=	Nreserve_Ltm1;
        NuLit_H=	NuLit_Htm1;
        NuLit_L=	NuLit_Ltm1;
        NupI_H=     NupI_Htm1;
        NupI_L=     NupI_Ltm1;
        Nuptake_H=	Nuptake_Htm1;
        Nuptake_L=	Nuptake_Ltm1;
        OH=         OHtm1;
        OL=         OLtm1;
        PARI_H=     PARI_Htm1;
        PARI_L=     PARI_Ltm1;
        PHE_S_H=	PHE_S_Htm1;
        PHE_S_L=	PHE_S_Ltm1;
        Preserve_H=	Preserve_Htm1;
        Preserve_L=	Preserve_Ltm1;
        Psi_l_H=	Psi_l_Htm1;
        Psi_l_L=	Psi_l_Ltm1;
        Psi_s_H=	Psi_s_Htm1;
        Psi_s_L=	Psi_s_Ltm1;
        Psi_x_H=	Psi_x_Htm1;
        Psi_x_L=	Psi_x_Ltm1;
        Puptake_H=	Puptake_Htm1;
        Puptake_L=	Puptake_Ltm1;
        RA_H=       RA_Htm1;
        RA_L=       RA_Ltm1;
        rap_H=      rap_Htm1;
        rap_L=      rap_Ltm1;
        RB_H=       RB_Htm1;
        rb_H=       rb_Htm1;
        RB_L=       RB_Ltm1;
        rb_L=       rb_Ltm1;
        Rdark_H=	Rdark_Htm1;
        Rdark_L=	Rdark_Ltm1;
        Rexmy_H=	Rexmy_Htm1;
        Rexmy_L=	Rexmy_Ltm1;
        Rg_H=       Rg_Htm1;
        Rg_L=       Rg_Ltm1;
        rKc_H=      rKc_Htm1;
        rKc_L=      rKc_Ltm1;
        Rmc_H=      Rmc_Htm1;
        Rmc_L=      Rmc_Ltm1;
        Rmr_H=      Rmr_Htm1;
        Rmr_L=      Rmr_Ltm1;
        Rms_H=      Rms_Htm1;
        Rms_L=      Rms_Ltm1;
        rNc_H=      rNc_Htm1;
        rNc_L=      rNc_Ltm1;
        rPc_H=      rPc_Htm1;
        rPc_L=      rPc_Ltm1;
        Rrootl_H=	Rrootl_Htm1;
        Rrootl_L=	Rrootl_Ltm1;
        rs_shdH=	rs_shdHtm1;
        rs_shdL=	rs_shdLtm1;
        rs_sunH=	rs_sunHtm1;
        rs_sunL=	rs_sunLtm1;
        SAI_H=      SAI_Htm1;
        SAI_L=      SAI_Ltm1;
        Sfr_H=      Sfr_Htm1;
        Sfr_L=      Sfr_Ltm1;
        SIF_H=      SIF_Htm1;
        SIF_L=      SIF_Ltm1;
        Slf_H=      Slf_Htm1;
        Slf_L=      Slf_Ltm1;
        Sll_H=      Sll_Htm1;
        Sll_L=      Sll_Ltm1;
        Sr_H=       Sr_Htm1;
        Sr_L=       Sr_Ltm1;
        SupK_H=     SupK_Htm1;
        SupK_L=     SupK_Ltm1;
        SupN_H=     SupN_Htm1;
        SupN_L=     SupN_Ltm1;
        SupP_H=     SupP_Htm1;
        SupP_L=     SupP_Ltm1;
        Swm_H=      Swm_Htm1;
        Swm_L=      Swm_Ltm1;
        T_H=        T_Htm1;
        T_L=        T_Ltm1;
        TBio_H=     TBio_Htm1;
        TBio_L=     TBio_Ltm1;
        Tden_H=     Tden_Htm1;
        Tden_L=     Tden_Ltm1;
        Tdp_H=      Tdp_Htm1;
        Tdp_L=      Tdp_Ltm1;
        TdpI_H=     TdpI_Htm1;
        TdpI_L=     TdpI_Ltm1;
        TexC_H=     TexC_Htm1;
        TexC_L=     TexC_Ltm1;
        TexK_H=     TexK_Htm1;
        TexK_L=     TexK_Ltm1;
        TexN_H=     TexN_Htm1;
        TexN_L=     TexN_Ltm1;
        TexP_H=     TexP_Htm1;
        TexP_L=     TexP_Ltm1;
        TNIT_H=     TNIT_Htm1;
        TNIT_L=     TNIT_Ltm1;
        TPHO_H=     TPHO_Htm1;
        TPHO_L=     TPHO_Ltm1;
        TPOT_H=     TPOT_Htm1;
        TPOT_L=     TPOT_Ltm1;
        Vl_H=       Vl_Htm1;
        Vl_L=       Vl_Ltm1;
        Vx_H=       Vx_Htm1;
        Vx_L=       Vx_Ltm1;
    end

    %% LOOP
    %======================================================================
    % LOOP OVER CELLS
    %======================================================================
    
    parfor ij=3550:4000   
        % Good practice to use a simple for loop for debugging/testing
        % ij is the index to go pixel by pixel through the mask
        % disp(strcat('in the loop', ij))

        if MASKn(ij)== 1
            Elev=DTM(ij);
            %[i,j] = ind2sub([m_cell,n_cell],ij);

            % BOUNDARY CONDITION  
            % INTRODUCED SOIL AND VEG. for ij
            [aR,Zs,EvL_Zs,Inf_Zs,Bio_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Ks_Zs,Dz,...
                ms,Kbot,Krock,zatm,Ccrown,Cbare,Crock,Curb,Cwat,...
                Color_Class,OM_H,OM_L,PFT_opt_H,PFT_opt_L,d_leaf_H,...
                d_leaf_L,SPAR,Phy,Soil_Param,Interc_Param,SnowIce_Param,...
                VegH_Param,VegL_Param,fpr,VegH_Param_Dyn,VegL_Param_Dyn,...
                Stoich_H,aSE_H,Stoich_L,aSE_L,fab_H,fbe_H,fab_L,fbe_L,...
                ZR95_H,ZR95_L,In_max_urb,In_max_rock,K_usle,Urb_Par,Deb_Par,...
                Zs_deb,Sllit,Kct,ExEM,ParEx_H,Mpar_H,ParEx_L,...
                Mpar_L,]=PARAMETERS_SOIL(ksv(ij),...
                PSAN(ij),PCLA(ij),PORG(ij),DEB_MAP(ij),md_max,...
                zatm_surface,Afirn(ij),SOIL_TH(ij),s);
            %%%
            
            % If hour (Datam_S(4)) is equal to 1
            if (Datam_S(4)==1)                                
                %% SOIL BIOGEOCHEMISTRY MODULE
                [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM]=Biogeo_environment([squeeze(Tdp_t(ij,:,:))]',[squeeze(O_t(ij,:,:))]',[squeeze(V_t(ij,:,:))]',...
                    Soil_Param,Phy,SPAR,Bio_Zs);%
                
                % Biogeochemistry Unit
                Nuptake_H(ij,:)= 0.0;
                Puptake_H(ij,:)= 0.0;
                Kuptake_H(ij,:)= 0.0; %% [gK/m^2 day]
                %%%
                Nuptake_L(ij,:)= 0.0; %% [gN/m^2 day]
                Puptake_L(ij,:)= 0.0;
                Kuptake_L(ij,:)= 0.0;
                %%%
                NavlI(ij,:)=[1 1 1];
                Bam(ij)=0; Bem(ij)=0;
                %%% 
                
                %% DEBUGGING: Create a structure and assign field names and values
                % For vegetation unit

                %{
                disp('here')
                if t == 3
                vg_Struct = struct(...
                    'cc_max', cc_max, ...
                    'Ccrown', Ccrown, ...
                    'ZR95_H', ZR95_H, ...
                    'ZR95_L', ZR95_L, ...
                    'B_Htm1', B_Htm1(ij,:,:), ...
                    'PHE_S_Htm1', PHE_S_Htm1(ij,:), ...
                    'dflo_Htm1', dflo_Htm1(ij,:), ...
                    'AgeL_Htm1', AgeL_Htm1(ij,:), ...
                    'AgeDL_Htm1', AgeDL_Htm1(ij,:), ...
                    'Ta_t', Ta_t(ij,:), ...
                    'PAR_t', PAR_t(ij,:), ...
                    'Tdp_H_t', Tdp_H_t(ij,:,:), ...
                    'Psi_x_H_t', Psi_x_H_t(ij,:,:), ...
                    'Psi_l_H_t', Psi_l_H_t(ij,:,:), ...
                    'An_H_t', An_H_t(ij,:,:), ...
                    'Rdark_H_t', Rdark_H_t(ij,:,:), ...
                    'NPP_Htm1', NPP_Htm1(ij,:), ...
                    'jDay', jDay, ...
                    'Datam_S', Datam_S, ...
                    'NPPI_Htm1', NPPI_Htm1(ij,:), ...
                    'TdpI_Htm1', TdpI_Htm1(ij,:), ...
                    'Bfac_weekHtm1', Bfac_weekHtm1(ij,:), ...
                    'Stoich_H', Stoich_H, ...
                    'aSE_H', aSE_H, ...
                    'VegH_Param_Dyn', VegH_Param_Dyn, ...
                    'Nreserve_Htm1', Nreserve_Htm1(ij,:), ...
                    'Preserve_Htm1', Preserve_Htm1(ij,:), ...
                    'Kreserve_Htm1', Kreserve_Htm1(ij,:), ...
                    'Nuptake_H', Nuptake_H(ij,:), ...
                    'Puptake_H', Puptake_H(ij,:), ...
                    'Kuptake_H', Kuptake_H(ij,:), ...
                    'FNC_Htm1', FNC_Htm1(ij,:), ...
                    'Tden_Htm1', Tden_Htm1(ij,:), ...
                    'AgePl_Htm1', AgePl_Htm1(ij,:), ...
                    'fab_H', fab_H, ...
                    'fbe_H', fbe_H, ...
                    'ParEx_H', ParEx_H, ...
                    'Mpar_H', Mpar_H, ...
                    'TBio_H', TBio_H(ij,:), ...
                    'SAI_Htm1', SAI_Htm1(ij,:), ...
                    'hc_Htm1', hc_Htm1(ij,:), ...
                    'B_Ltm1', B_Ltm1(ij,:,:), ...
                    'PHE_S_Ltm1', PHE_S_Ltm1(ij,:), ...
                    'dflo_Ltm1', dflo_Ltm1(ij,:), ...
                    'AgeL_Ltm1', AgeL_Ltm1(ij,:), ...
                    'AgeDL_Ltm1', AgeDL_Ltm1(ij,:), ...
                    'Tdp_L_t', Tdp_L_t(ij,:,:), ...
                    'Psi_x_L_t', Psi_x_L_t(ij,:,:), ...
                    'Psi_l_L_t', Psi_l_L_t(ij,:,:), ...
                    'An_L_t', An_L_t(ij,:,:), ...
                    'Rdark_L_t', Rdark_L_t(ij,:,:), ...
                    'NPP_Ltm1', NPP_Ltm1(ij,:), ...
                    'NPPI_Ltm1', NPPI_Ltm1(ij,:), ...
                    'TdpI_Ltm1', TdpI_Ltm1(ij,:), ...
                    'Bfac_weekLtm1', Bfac_weekLtm1(ij,:), ...
                    'NupI_Htm1', NupI_Htm1(ij,:,:), ...
                    'NupI_Ltm1', NupI_Ltm1(ij,:,:), ...
                    'NuLit_Htm1', NuLit_Htm1(ij,:,:), ...
                    'NuLit_Ltm1', NuLit_Ltm1(ij,:,:), ...
                    'NBLeaf_Htm1', NBLeaf_Htm1(ij,:), ...
                    'NBLeaf_Ltm1', NBLeaf_Ltm1(ij,:), ...
                    'PARI_Htm1', PARI_Htm1(ij,:,:), ...
                    'NBLI_Htm1', NBLI_Htm1(ij,:), ...
                    'PARI_Ltm1', PARI_Ltm1(ij,:,:), ...
                    'NBLI_Ltm1', NBLI_Ltm1(ij,:), ...
                    'Stoich_L', Stoich_L, ...
                    'aSE_L', aSE_L, ...
                    'VegL_Param_Dyn', VegL_Param_Dyn, ...
                    'Nreserve_Ltm1', Nreserve_Ltm1(ij,:), ...
                    'Preserve_Ltm1', Preserve_Ltm1(ij,:), ...
                    'Kreserve_Ltm1', Kreserve_Ltm1(ij,:), ...
                    'Nuptake_L', Nuptake_L(ij,:), ...
                    'Puptake_L', Puptake_L(ij,:), ...
                    'Kuptake_L', Kuptake_L(ij,:), ...
                    'FNC_Ltm1', FNC_Ltm1(ij,:), ...
                    'Tden_Ltm1', Tden_Ltm1(ij,:), ...
                    'AgePl_Ltm1', AgePl_Ltm1(ij,:), ...
                    'fab_L', fab_L, ...
                    'fbe_L', fbe_L, ...
                    'ParEx_L', ParEx_L, ...
                    'Mpar_L', Mpar_L, ...
                    'TBio_L', TBio_L(ij,:), ...
                    'SAI_Ltm1', SAI_Ltm1(ij,:), ...
                    'hc_Ltm1', hc_Ltm1(ij,:), ...
                    'ExEM', ExEM, ...
                    'Lat', Lat, ...
                    'Se_bio', Se_bio, ...
                    'Tdp_bio', Tdp_bio, ...
                    'OPT_EnvLimitGrowth', OPT_EnvLimitGrowth, ...
                    'OPT_VD', OPT_VD, ...
                    'OPT_VCA', OPT_VCA, ...
                    'OPT_ALLOME', OPT_ALLOME, ...
                    'OPT_SoilBiogeochemistry', OPT_SoilBiogeochemistry);
                save("M:\19_ISTA\1_TC\3_Model_Source\TeC_Source_Code\Model_comparison\Cinca_vg_param.mat", 'vg_Struct')
                end
                %}
               
                %% VEGETATION MODULE
                %==========================================================
                % FUNCTION: VEGETATION_MODULE_PAR
                %==========================================================

                [LAI_H(ij,:),B_H(ij,:,:),NPP_H(ij,:),ANPP_H(ij,:),Rg_H(ij,:),RA_H(ij,:),...
                    Rms_H(ij,:),Rmr_H(ij,:),Rmc_H(ij,:),PHE_S_H(ij,:),...
                    dflo_H(ij,:),AgeL_H(ij,:),e_rel_H(ij,:),e_relN_H(ij,:),...
                    LAI_L(ij,:),B_L(ij,:,:),NPP_L(ij,:),ANPP_L(ij,:),Rg_L(ij,:),RA_L(ij,:),...
                    Rms_L(ij,:),Rmr_L(ij,:),Rmc_L(ij,:),PHE_S_L(ij,:),...
                    dflo_L(ij,:),AgeL_L(ij,:),e_rel_L(ij,:),e_relN_L(ij,:),...
                    SAI_H(ij,:),hc_H(ij,:),SAI_L(ij,:),hc_L(ij,:),...
                    LAIdead_H(ij,:),NBLeaf_H(ij,:),Sr_H(ij,:),Slf_H(ij,:),Sfr_H(ij,:),Sll_H(ij,:),Swm_H(ij,:),Rexmy_H(ij,:,:),NupI_H(ij,:,:),NuLit_H(ij,:,:),...
                    LAIdead_L(ij,:),NBLeaf_L(ij,:),Sr_L(ij,:),Slf_L(ij,:),Sfr_L(ij,:),Sll_L(ij,:),Swm_L(ij,:),Rexmy_L(ij,:,:),NupI_L(ij,:,:),NuLit_L(ij,:,:),...
                    Rrootl_H(ij,:),AgeDL_H(ij,:),Bfac_dayH(ij,:),Bfac_weekH(ij,:),NPPI_H(ij,:),TdpI_H(ij,:),PARI_H(ij,:,:),NBLI_H(ij,:),RB_H(ij,:,:),FNC_H(ij,:),...
                    Nreserve_H(ij,:),Preserve_H(ij,:),Kreserve_H(ij,:),rNc_H(ij,:),rPc_H(ij,:),rKc_H(ij,:),ManIH(ij,:),...
                    Rrootl_L(ij,:),AgeDL_L(ij,:),Bfac_dayL(ij,:),Bfac_weekL(ij,:),NPPI_L(ij,:),TdpI_L(ij,:),PARI_L(ij,:,:),NBLI_L(ij,:),RB_L(ij,:,:),FNC_L(ij,:),...
                    Nreserve_L(ij,:),Preserve_L(ij,:),Kreserve_L(ij,:),rNc_L(ij,:),rPc_L(ij,:),rKc_L(ij,:),ManIL(ij,:),...
                    TexC_H(ij,:),TexN_H(ij,:),TexP_H(ij,:),TexK_H(ij,:),TNIT_H(ij,:),TPHO_H(ij,:),TPOT_H(ij,:),...
                    SupN_H(ij,:),SupP_H(ij,:),SupK_H(ij,:),ISOIL_H(ij,:,:),...
                    TexC_L(ij,:),TexN_L(ij,:),TexP_L(ij,:),TexK_L(ij,:),TNIT_L(ij,:),TPHO_L(ij,:),TPOT_L(ij,:),...
                    SupN_L(ij,:),SupP_L(ij,:),SupK_L(ij,:),ISOIL_L(ij,:,:),...
                    BA_H(ij,:),Tden_H(ij,:),AgePl_H(ij,:),BA_L(ij,:),Tden_L(ij,:),AgePl_L(ij,:),Ccrown_t(ij,:)]=VEGETATION_MODULE_PAR( ...
                    cc_max, ...
                    Ccrown, ...
                    ZR95_H, ...
                    ZR95_L, ...
                    B_Htm1(ij,:,:),...
                    PHE_S_Htm1(ij,:), ...
                    dflo_Htm1(ij,:), ...
                    AgeL_Htm1(ij,:), ...
                    AgeDL_Htm1(ij,:), ...
                    Ta_t(ij,:), ...
                    PAR_t(ij,:), ...
                    Tdp_H_t(ij,:,:), ...
                    Psi_x_H_t(ij,:,:), ...
                    Psi_l_H_t(ij,:,:), ...
                    An_H_t(ij,:,:), ...
                    Rdark_H_t(ij,:,:), ...
                    NPP_Htm1(ij,:), ...
                    jDay, ...
                    Datam_S,...
                    NPPI_Htm1(ij,:), ...
                    TdpI_Htm1(ij,:), ...
                    Bfac_weekHtm1(ij,:), ...
                    Stoich_H, ...
                    aSE_H, ...
                    VegH_Param_Dyn,...
                    Nreserve_Htm1(ij,:), ...
                    Preserve_Htm1(ij,:), ...
                    Kreserve_Htm1(ij,:), ...
                    Nuptake_H(ij,:), ...
                    Puptake_H(ij,:), ...
                    Kuptake_H(ij,:), ...
                    FNC_Htm1(ij,:), ...
                    Tden_Htm1(ij,:), ...
                    AgePl_Htm1(ij,:),...
                    fab_H, ...
                    fbe_H, ...
                    ParEx_H, ...
                    Mpar_H, ...
                    TBio_H(ij,:), ...
                    SAI_Htm1(ij,:), ...
                    hc_Htm1(ij,:),...
                    B_Ltm1(ij,:,:), ...
                    PHE_S_Ltm1(ij,:), ...
                    dflo_Ltm1(ij,:), ...
                    AgeL_Ltm1(ij,:), ...
                    AgeDL_Ltm1(ij,:),...
                    Tdp_L_t(ij,:,:), ...
                    Psi_x_L_t(ij,:,:), ...
                    Psi_l_L_t(ij,:,:), ...
                    An_L_t(ij,:,:), ...
                    Rdark_L_t(ij,:,:), ...
                    NPP_Ltm1(ij,:),...
                    NPPI_Ltm1(ij,:), ...
                    TdpI_Ltm1(ij,:), ...
                    Bfac_weekLtm1(ij,:),...
                    NupI_Htm1(ij,:,:), ...
                    NupI_Ltm1(ij,:,:), ...
                    NuLit_Htm1(ij,:,:), ...
                    NuLit_Ltm1(ij,:,:), ...
                    NBLeaf_Htm1(ij,:), ...
                    NBLeaf_Ltm1(ij,:),...
                    PARI_Htm1(ij,:,:), ...
                    NBLI_Htm1(ij,:), ...
                    PARI_Ltm1(ij,:,:), ...
                    NBLI_Ltm1(ij,:),...
                    Stoich_L, ...
                    aSE_L, ...
                    VegL_Param_Dyn,...
                    NavlI(ij,:), ...
                    Bam(ij), ...
                    Bem(ij), ...
                    Ccrown_t_tm1(ij,:),...
                    Nreserve_Ltm1(ij,:), ...
                    Preserve_Ltm1(ij,:), ...
                    Kreserve_Ltm1(ij,:), ...
                    Nuptake_L(ij,:), ...
                    Puptake_L(ij,:), ...
                    Kuptake_L(ij,:), ...
                    FNC_Ltm1(ij,:), ...
                    Tden_Ltm1(ij,:), ...
                    AgePl_Ltm1(ij,:),...
                    fab_L, ...
                    fbe_L, ...
                    ParEx_L, ...
                    Mpar_L, ...
                    TBio_L(ij,:), ...
                    SAI_Ltm1(ij,:), ...
                    hc_Ltm1(ij,:),...
                    ExEM, ...
                    Lmax_day, ...
                    L_day, ...                    
                    Se_bio, ...
                    Tdp_bio, ...
                    OPT_EnvLimitGrowth, ...
                    OPT_VD, ...
                    OPT_VCA, ...
                    OPT_ALLOME, ...
                    OPT_SoilBiogeochemistry);
                
                BLit(ij,:)= 0.0 ; % %% %%[kg DM / m2]
            end        

            %% HYDROLOGY MODULE
            %==============================================================
            % FUNTION: HYDROLOGY_MODULE_PAR
            %==============================================================

            [V(ij,:),O(ij,:),Vice(ij,:),Oice(ij,:),ZWT(ij),OF(ij),OS(ij),OH(ij,:),OL(ij,:),Psi_s_H(ij,:),Psi_s_L(ij,:),Rd(ij),Qi_out(ij,:),...
                Rh(ij),Lk(ij),f(ij),WIS(ij),Ts(ij),Csno(ij),Cice(ij),NDVI(ij),...
                Pr_sno(ij),Pr_liq(ij),rb_H(ij,:),rb_L(ij,:),rs_sunH(ij,:),rs_sunL(ij,:),rs_shdH(ij,:),rs_shdL(ij,:),...
                rap_H(ij,:),rap_L(ij,:),r_soil(ij),b_soil(ij),alp_soil(ij),ra(ij),r_litter(ij,:),...
                WR_SP(ij),U_SWE(ij),NIn_SWE(ij),dQ_S(ij),DQ_S(ij),DT_S(ij),...
                WAT(ij),ICE(ij),ICE_D(ij),IP_wc(ij),WR_IP(ij),NIce(ij),Cicew(ij),Csnow(ij),FROCK(ij),...
                Dr_H(ij,:),Dr_L(ij,:),SE_rock(ij),SE_urb(ij),Lk_wat(ij),Lk_rock(ij),...
                An_L(ij,:),An_H(ij,:),Rdark_L(ij,:),Rdark_H(ij,:),Ci_sunH(ij,:),Ci_sunL(ij,:),Ci_shdH(ij,:),Ci_shdL(ij,:),Rn(ij),...
                H(ij),QE(ij),Qv(ij),Lpho(ij),T_H(ij,:),T_L(ij,:),EIn_H(ij,:),EIn_L(ij,:),EG(ij),ELitter(ij),ESN(ij),ESN_In(ij),...
                EWAT(ij),EICE(ij),EIn_urb(ij),EIn_rock(ij),dw_SNO(ij),Imelt(ij),Smelt(ij),...
                G(ij),Gfin(ij),Tdp(ij,:),Tdp_snow(ij,:),Tdeb(ij,:),Tice(ij,:),Tdamp(ij),Tdp_H(ij,:),Tdp_L(ij,:),SWE(ij),SND(ij),ros(ij),In_SWE(ij),SP_wc(ij),Qfm(ij),t_sls(ij),...
                In_H(ij,:),In_L(ij,:),In_Litter(ij),In_urb(ij),In_rock(ij),...
                gsr_H(ij,:),Psi_x_H(ij,:),Psi_l_H(ij,:),Jsx_H(ij,:),Jxl_H(ij,:),Kleaf_H(ij,:),Kx_H(ij,:),Vx_H(ij,:),Vl_H(ij,:),...
                gsr_L(ij,:),Psi_x_L(ij,:),Psi_l_L(ij,:),Jsx_L(ij,:),Jxl_L(ij,:),Kleaf_L(ij,:),Kx_L(ij,:),Vx_L(ij,:),Vl_L(ij,:),...
                fapar_H(ij,:),fapar_L(ij,:),SIF_H(ij,:),SIF_L(ij,:),...
                Ws_under(ij),er(ij),snow_albedo(ij,:),tau_sno(ij),e_sno(ij),...
                HV(ij),QEV(ij),dQVEG(ij),TsVEG(ij),Ts_under(ij),EK(ij),POT(ij,:),CK1(ij)]=HYDROLOGY_MODULE_PAR(Vtm1(ij,:),Oicetm1(ij,:),...
                aR,Zs,...
                EvL_Zs,Inf_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Dz,ms,Kbot,Pr_S(ij),Ta_S(ij),Ds_S(ij),Ws_S(ij),zatm,Tstm1(ij),dt,dth,ea_S(ij),N_S(ij),Pre_S(ij),Tstm0(ij),...
                LAI_H(ij,:),SAI_H(ij,:),LAI_L(ij,:),SAI_L(ij,:),LAIdead_H(ij,:),LAIdead_L(ij,:),...
                Rrootl_H(ij,:),Rrootl_L(ij,:),BLit(ij,:),Sllit,Kct,...
                Datam_S,DeltaGMT,Lon,Lat,t_bef,t_aft,...
                Ccrown,Cbare,Crock,Curb,Cwat,...
                SAB1_S(ij),SAB2_S(ij),SAD1_S(ij),SAD2_S(ij),PARB_S(ij),PARD_S(ij),SvF(ij),SNDtm1(ij),snow_albedotm1(ij,:),Color_Class,OM_H,OM_L, ...
                PFT_opt_H,PFT_opt_L,hc_H(ij,:),hc_L(ij,:),d_leaf_H,d_leaf_L,...
                Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,...
                Ca_S(ij),Oa,Citm1_sunH(ij,:),Citm1_shdH(ij,:),Citm1_sunL(ij,:),Citm1_shdL(ij,:),...
                e_rel_H(ij,:),e_relN_H(ij,:),e_rel_L(ij,:),e_relN_L(ij,:),...
                e_snotm1(ij),In_Htm1(ij,:),In_Ltm1(ij,:),In_Littertm1(ij),In_urbtm1(ij),In_rocktm1(ij),SWEtm1(ij),In_SWEtm1(ij),....
                Tdebtm1(ij,:),Ticetm1(ij),Tdptm1(ij,:),Tdp_snowtm1(ij,:),Tdamptm1(ij),Ts_undertm1(ij),...
                WATtm1(ij),ICEtm1(ij),IP_wctm1(ij),ICE_Dtm1(ij),Cicewtm1(ij),...
                Vx_Htm1(ij,:),Vl_Htm1(ij,:),Vx_Ltm1(ij,:),Vl_Ltm1(ij,:),Psi_x_Htm1(ij,:),Psi_l_Htm1(ij,:),Psi_x_Ltm1(ij,:),Psi_l_Ltm1(ij,:),...
                ZR95_H,ZR95_L,...
                FROCKtm1(ij),Krock,...
                Urb_Par,Deb_Par,Zs_deb,...
                Tdew_S(ij),t_slstm1(ij),rostm1(ij),SP_wctm1(ij),fpr,IrD_S(ij),...
                In_max_urb,In_max_rock,K_usle,tau_snotm1(ij),Ta_day(ij,:),...
                Slo_top(ij),Slo_head(ij,:),Asur(ij),Ared(ij),aTop(ij),EKtm1(ij),q_runon(ij),Qi_in(ij,:),...
                Ws_undertm1(ij),Pr_sno_t(ij,:),...
                pow_dis,a_dis,Salt_S(ij),...
                SPAR,SNn(ij),OPT_min_SPD,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST,Opt_ST2,OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH,parameterize_phase, hSTL, Albsno_method);

        end
    end
    

    %%post-compute sublimation from ESN, do this inside hydrology module
    %%once it turns out to be useful
    SSN = ESN.*(Ts<0);
    SSN_In = ESN_In.*(Ts<0);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii=1:ms_max
        Qi_out_Rout(:,:,ii)=   reshape(Qi_out(:,ii),m_cell,n_cell); %[mm/h]
    end
    Rd= reshape(Rd,m_cell,n_cell); %%[mm]
    Rh= reshape(Rh,m_cell,n_cell); %%[mm]
    
    %% ROUTING MODULE
    %======================================================================
    % Funtion: ROUTING_MODULE
    %======================================================================
    
    [q_runon,Q_channel,Qi_in_Rout,Slo_pot,Q_exit,Qsub_exit,T_pot,...
        QpointH,QpointC,UpointH,UpointC]= ROUTING_MODULE(dt,dth,Rd,Rh,...
        Qi_out_Rout,Q_channel,cellsize,Area,DTM,NMAN_H,NMAN_C,MRough,...
        WC,SN,T_flow,T_pot,Slo_top,ms_max,POT,ZWT,OPT_HEAD,Xout,Yout);
 
    
    for ii=1:ms_max
        Qi_in(:,ii)=	reshape(Qi_in_Rout(:,:,ii),num_cell,1); %[mm]
        Slo_head(:,ii)=	reshape(Slo_pot(:,:,ii),num_cell,1);
    end
    Rd = reshape(Rd,num_cell,1); %%[mm]
    Rh = reshape(Rh,num_cell,1); %%[mm]
    q_runon = reshape(q_runon,num_cell,1); %%[mm]
    Qi_in=Qi_in/dth;%% [mm/h]
    q_runon = q_runon/dth; %%% [mm/h]
    %%% Q_exit Qsub_exit [mm] over the entire domain
    if not(isreal(sum(sum(q_runon))))
        disp('The Program fails because of runon numerical instability')
        break
    end

    
    %% Glacier volume redistribution 

if (Datam_S(2)==1) && (Datam_S(3)==1) && (Datam_S(4)==1) && (Idyn > 0)
        
        GLA = (GLH > 0) & (MASK == 1); % Glacier mask

        % sum up total water equivalent in [m w.e.], store ratios
        WEbef = ICE + SWE; %total w.e. before redistribution
        WEym1 = ICEym1 + SWEym1; %total w.e. from the year before

        raICE = ICE./WEbef; raICE(isnan(raICE))=0; %ice ratio of total
        raSWE = SWE./WEbef; raSWE(isnan(raSWE))=0; %snow ratio of total mm w.e.

        THbef = reshape(WEbef,m_cell,n_cell)./916; %conversion to thickness in m, rhi parameter needed, and reshape
        THym1 = reshape(WEym1,m_cell,n_cell)./916; %conversion to thickness in m, rhi parameter needed, and reshape

        THbef(GLA ~= 1) = 0;
        THym1(GLA ~= 1) = 0;
        
        SMB = (WEbef - WEym1)./916; % Distributed surface mass balance of the previous year in m/y ice-eq (IGM)
        SMB(GLH == 0 | MASK~=1 | GLA ~= 1) = -10;
        SMB(isnan(SMB)) = -10;
        SMB = reshape(SMB,m_cell,n_cell);
 
        %%%% Create or update the .nc file %%%%%%%%%%%%

        x(2) = x(1) + (y(2)-y(1)); % Temporary fix
        IGM_netcdf([outlocation 'igm_inputs3.nc'],x,y,rot90(flipud(double(GLA)),3),rot90(flipud(THym1),3),rot90(flipud(DTM_Bedrock),3)...
            , rot90(flipud(DTM_orig),3),rot90(flipud(SMB),3)) % Create .nc file

        %%%% Run IGM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        pyenv('Version',path_igmEnv, 'ExecutionMode','OutOfProcess');   % setup python virtual env
        igm_run_path =  [path_igm '/igm/igm_run.py']; % path of IGM .py script
        addpath(genpath([path_igm '/igm/']))
 
        igm_arg = [' --lncd_input_file ' [outlocation 'igm_inputs3.nc']...
            ' --working_dir ' outlocation...
            ' --time_start ' num2str(1)...
            ' --time_end ' num2str(2) ... 
            ' --time_save ' num2str(1)];

        pyrunfile([igm_run_path igm_arg]) % Run IGM

        thk = ncread([outlocation 'output.nc'],'thk'); % Load updated glacier thickness from IGM

        THnew = double(rot90(thk(:,:,2))); % New thickness after IGM

        %Check IGM mass conservation
        dMB_igm = nansum(THnew,'all') - nansum(THbef,'all');

        % Add the icethickness which leaks outside of glacier mask back on existing glacier
        IGM_smb = nansum(THnew,'all') - nansum(THym1,'all');
        TC_smb = nansum(SMB(GLA>0),'all');

        THnew(flipud(GLA == 1) & THnew > 0) = THnew(flipud(GLA == 1) & THnew > 0) + (TC_smb-IGM_smb)./nansum(GLA>0,'all');
        THnew(flipud(GLA ~= 1)) = 0; %No glacier growing in area

        % reshape total water equivalent back to array map
        WEnew = reshape(flipud(THnew),num_cell,1).*916; % conversion back to mm w.e., rhi parameter needed, and reshape
        
        % recalculate snow/ice w.e. maps based on old ratios
        ICE = raICE.*WEnew; % re-calculate ice w.e.
        SWE(WEnew>0) = raSWE(WEnew>0).*WEnew(WEnew>0); % re-calculate SWE, only where glacier

        %Mass balance check
        dMass = nansum(ICE(WEnew >0) + SWE(WEnew >0)) - nansum(WEbef(WEnew >0));

        disp(['Total SMB: ' num2str(TC_smb)])
        disp(['SMB after IGM redistribution: ' num2str(dMass/916)])
        
        % store for next year redistribution
        ICEym1 = ICE;     
        SWEym1 = SWE;

        % calculate new snow/ice depths
        ICE_D = ICE/916; ICE_D(isnan(ICE_D))=0;
        SND = SWE./ros; SND(isnan(SND))=0;
end


    %% AVALANCHES COMPONENT

    SND=    reshape(SND,m_cell,n_cell); SND(isnan(SND)) = 0;
    SWE=    reshape(SWE,m_cell,n_cell); SWE(isnan(SWE)) = 0;
    ros=    reshape(ros,m_cell,n_cell); ros(isnan(ros)) = 0;

    if Aval == 1
        SWEpreava = SWE; 

        [SND,SWE,ros,Swe_exit]= AVALANCHES(DTM,cellsize,Area,...
        reshape(Asur,m_cell,n_cell),Slo_top,SND,SWE,ros,a_aval,C_aval);
        SWE_avalanched = SWE-SWEpreava;   
        clear SWEpreava;
    else 
        SWE_avalanched = SWE.*0;
        Swe_exit = SWE.*0;
    end 

    %% Catchment average snowline    
    SND=	reshape(SND,num_cell,1); SND(isnan(SND)) = 0;
    SWE=    reshape(SWE,num_cell,1); SWE(isnan(SWE)) = 0;
    ros=    reshape(ros,num_cell,1); ros(isnan(ros)) = 0;
    SWE_avalanched = reshape(SWE_avalanched,num_cell,1);

    %% FRACTURED ROCK COMPONENT   
    [Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK(Q_channel,FROCK,...
        SPRINGn,dth,m_cell,n_cell,num_cell,Kres_Rock);

    if t==2
        V_tgtm1=    sum(Ared.*Asur.*sum(Vtm1,2))*(cellsize^2)/Area;
        Vice_tgtm1= sum(Ared.*Asur.*sum(Vicetm1,2))*(cellsize^2)/Area;
        SWE_tgtm1=  sum(SWEtm1)*(cellsize^2)/Area;
        In_tgtm1=   (sum(sum(In_Htm1)) + sum(sum(In_Ltm1)) + ...
            sum(sum(In_Littertm1)) + sum(SP_wctm1) + ...
            sum(In_SWEtm1) + sum(In_urbtm1) + ...
            sum(In_rocktm1) + sum(IP_wctm1) ) * (cellsize^2)/Area;
        ICE_tgtm1=  sum(ICEtm1)*(cellsize^2)/Area; %%
        WAT_tgtm1=  sum(WATtm1)*(cellsize^2)/Area; %%
        FROCK_tgtm1=sum(FROCKtm1)*(cellsize^2)/Area; %%
    end
       
    %% ALBEDO MAP

    if t==2 
        snow_albedo_out = snow_albedotm1;
        surface_albedo_out = surface_albedotm1;
    elseif Datam_S(4)==12
    snow_albedo_out = snow_albedo;
    surface_albedo_out = surface_albedotm1;
    surface_albedo_out(ICE_D>0)=0.28;
    surface_albedo_out(DEB_MAP>0)=0.13;
    surface_albedo_out(SND>0)=snow_albedo(SND>0);
    end

    %% INITIAL CONDITION FOR THE NEXT STEP
    %======================================================================
    %
    %======================================================================
    Cicewtm1=   Cicew;
    Citm1_shdH= Ci_shdH;
    Citm1_shdL= Ci_shdL;
    Citm1_sunH= Ci_sunH;
    Citm1_sunL= Ci_sunL;
    e_snotm1=   e_sno;
    EKtm1=      EK;
    FROCKtm1=   FROCK;
    ICE_Dtm1=   ICE_D;
    ICEtm1=     ICE;
    In_Htm1=    In_H ;
    In_Littertm1=In_Litter;
    In_Ltm1=    In_L;
    In_rocktm1= In_rock;
    In_SWEtm1=  In_SWE;
    In_urbtm1=  In_urb;
    IP_wctm1=   IP_wc;
    Oicetm1=    Oice;
    Psi_l_Htm1=	Psi_l_H;
    Psi_l_Ltm1=	Psi_l_L;
    Psi_x_Htm1=	Psi_x_H;
    Psi_x_Ltm1=	Psi_x_L;
    rostm1=     ros;
    SNDtm1=     SND ;
    snow_albedotm1= snow_albedo ;
    soil_albedotm1= soil_albedo ;
    SP_wctm1=   SP_wc ;
    surface_albedotm1= surface_albedo ;
    SWEtm1=     SWE ;
    SWE_avalanchedtm1= SWE_avalanched ;
    t_slstm1=   t_sls;
    tau_snotm1= tau_sno;
    Tdamptm1=   Tdamp;
    Tdebtm1=    Tdeb;
    Tdp_snowtm1 = Tdp_snow;
    Tdptm1=     Tdp;
    Ticetm1=    Tice;
    %%% order!
    Tstm1=      Ts;
    Tstm0=      2*Ts-Tstm1;
    Ts_undertm1= Ts_under;
    %%%
    Vl_Htm1=    Vl_H;
    Vl_Ltm1=    Vl_L;
    Vtm1=       V;
    Vx_Htm1=    Vx_H;
    Vx_Ltm1=    Vx_L;
    WATtm1=     WAT;
    Ws_undertm1=Ws_under;
    
    if (Datam_S(4)==1)
        AgeDL_Htm1=     AgeDL_H;
        AgeDL_Ltm1=     AgeDL_L;
        AgeL_Htm1=      AgeL_H;
        AgeL_Ltm1=      AgeL_L ;
        AgePl_Htm1=     AgePl_H;
        AgePl_Ltm1=     AgePl_L;
        B_Htm1=         B_H;
        B_Ltm1=         B_L;
        Bfac_weekHtm1=  Bfac_weekH;
        Bfac_weekLtm1=  Bfac_weekL;
        Ccrown_t_tm1=   Ccrown_t;
        dflo_Htm1=      dflo_H;
        dflo_Ltm1=      dflo_L ;
        FNC_Htm1=       FNC_H;
        FNC_Ltm1=       FNC_L;
        hc_Htm1=        hc_H;
        hc_Ltm1=        hc_L;
        Kreserve_Htm1=  Kreserve_H;
        Kreserve_Ltm1=  Kreserve_L;
        NBLeaf_Htm1=    NBLeaf_H;
        NBLeaf_Ltm1=    NBLeaf_L;
        NBLI_Htm1=      NBLI_H;
        NBLI_Ltm1=      NBLI_L;
        NPP_Htm1=       NPP_H;
        NPP_Ltm1=       NPP_L;
        NPPI_Htm1=      NPPI_H;
        NPPI_Ltm1=      NPPI_L;
        Nreserve_Htm1=  Nreserve_H;
        Nreserve_Ltm1=  Nreserve_L;
        NuLit_Htm1=     NuLit_H;
        NuLit_Ltm1=     NuLit_L;
        NupI_Htm1=      NupI_H;
        NupI_Ltm1=      NupI_L;
        PARI_Htm1=      PARI_H;
        PARI_Ltm1=      PARI_L;
        PHE_S_Htm1=     PHE_S_H;
        PHE_S_Ltm1=     PHE_S_L;
        Preserve_Htm1=  Preserve_H;
        Preserve_Ltm1=  Preserve_L;
        SAI_Htm1=       SAI_H;
        SAI_Ltm1=       SAI_L;
        Tden_Htm1=      Tden_H;
        Tden_Ltm1=      Tden_L;
        TdpI_Htm1=      TdpI_H;
        TdpI_Ltm1=      TdpI_L;
    end
    
    
    %% MEMORY CONDITION FOR VEGETATION  MODEL     
    %%% 1 day
    if t > 24
        An_H_t(:,:,1:23)=       An_H_t(:,:,2:24);
        An_H_t(:,:,24)=         An_H;
        An_L_t(:,:,1:23)=       An_L_t(:,:,2:24);
        An_L_t(:,:,24)=         An_L;
        O_t(:,:,1:23)=          O_t(:,:,2:24);
        O_t(:,:,24)=            O;
        PAR_t(:,1:23)=          PAR_t(:,2:24);
        PAR_t(:,24)=            PARB_S + PARD_S;
        Pr_sno_t(:,1:23)=       Pr_sno_t(:,2:24);
        Pr_sno_t(:,24)=         Pr_sno;
        Psi_l_H_t(:,:,1:23)=    Psi_l_H_t(:,:,2:24);
        Psi_l_H_t(:,:,24)=      Psi_l_H;
        Psi_l_L_t(:,:,1:23)=    Psi_l_L_t(:,:,2:24);
        Psi_l_L_t(:,:,24)=      Psi_l_L;
        Psi_x_H_t(:,:,1:23)=    Psi_x_H_t(:,:,2:24);
        Psi_x_H_t(:,:,24)=      Psi_x_H;
        Psi_x_L_t(:,:,1:23)=    Psi_x_L_t(:,:,2:24);
        Psi_x_L_t(:,:,24)=      Psi_x_L;
        Rdark_H_t(:,:,1:23)=    Rdark_H_t(:,:,2:24);
        Rdark_H_t(:,:,24)=      Rdark_H;
        Rdark_L_t(:,:,1:23)=    Rdark_L_t(:,:,2:24);
        Rdark_L_t(:,:,24)=      Rdark_L;
        Ta_t(:,1:23)=           Ta_t(:,2:24);
        Ta_t(:,24)=             Ta_S;
        Tdp_H_t(:,:,1:23)=      Tdp_H_t(:,:,2:24);
        Tdp_H_t(:,:,24)=        Tdp_H;
        Tdp_L_t(:,:,1:23)=      Tdp_L_t(:,:,2:24);
        Tdp_L_t(:,:,24)=        Tdp_L;
        Tdp_t(:,:,1:23)=        Tdp_t(:,:,2:24);
        Tdp_t(:,:,24)=          Tdp;
        V_t(:,:,1:23)=          V_t(:,:,2:24);
        V_t(:,:,24)=            V;
    else
        An_H_t(:,:,t)=          An_H;
        An_L_t(:,:,t)=          An_L;
        O_t(:,:,t)=             O;
        PAR_t(:,24)=            PARB_S + PARD_S;
        Pr_sno_t(:,t)=          Pr_sno;
        Psi_l_H_t(:,:,t)=       Psi_l_H;
        Psi_l_L_t(:,:,t)=       Psi_l_L;
        Psi_x_H_t(:,:,t)=       Psi_x_H;
        Psi_x_L_t(:,:,t)=       Psi_x_L;
        Rdark_H_t(:,:,t)=       Rdark_H;
        Rdark_L_t(:,:,t)=       Rdark_L;
        Ta_t(:,t)=              Ta_S;
        Tdp_H_t(:,:,t)=         Tdp_H;
        Tdp_L_t(:,:,t)=         Tdp_L;
        Tdp_t(:,:,t)=           Tdp;
        V_t(:,:,t)=             V;
    end
    

    
    %% OUTPUT WRITING
    %======================================================================
    %
    %======================================================================
       
    run(['OUTPUT_MANAGER_DIST_LABEL.m']);

    % Save the workspace at frequent interval. Very useful in case it crashes 
    if  mod(t,25)==0    
        save([outlocation, Fstep], '-regexp', '^(?!(FF_BC|LWIN_BC|PARB|PARD|PP_BC|PRESS_BC|RH_BC|SAB1|SAB2|SAD1|SAD2|TA_BC|WS)$).');
    end   

    if  mod(t,8760)==0  ||  t==N_time_step
        Fstep2= strcat(Fstep,'_',num2str(t));
        save([outlocation, Fstep2], '-regexp', '^(?!(FF_BC|LWIN_BC|PARB|PARD|PP_BC|PRESS_BC|RH_BC|SAB1|SAB2|SAD1|SAD2|TA_BC|WS)$).');
    end


end
%close(bau)

%% Display
Computational_Time =toc;
disp('COMPUTATIONAL TIME [h] ')
disp(Computational_Time/3600)
disp(' COMPUTATIONAL TIME [s/cycle] ')
disp(Computational_Time/N_time_step)
