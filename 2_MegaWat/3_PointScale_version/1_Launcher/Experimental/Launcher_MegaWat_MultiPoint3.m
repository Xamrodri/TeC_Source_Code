%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%              POINT SCALE MODEL                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
% Created on Nov 25, 2024
% Author: MAXIMILIANO RODRIGUEZ
% Code originally from: ACHILLE JOUBERTON
% Area of Study: Apennines - Tiber River
% Region: Apennines
% Code explanation: This code launches the Point Scale version of TC model.
% Output variables are managed by the file OUTPUT_MANAGER_DIST_LABEL.
% Some code and names still refers to previous versions of the code.
%==========================================================================

%% Clear all
clc; clear all;

%% Directories
folder_path = 'C:/Users/mrodrigu/Desktop/19_ISTA/1_TC/3_Model_Source/2_MegaWat/'; % Put here the path of where you downloaded the repository
forc_path = 'C:/Users/mrodrigu/Desktop/19_ISTA/7_Forcing/2_Extractions_Point/3_Radiation_Partition/'; % Put here the path of where you downloaded the repository

%% Catchment selection
sel_basin = 2; % ID for catchment selection for variable SITE

%% Pixel selection and time
%==========================================================================
% Depends on the point to be modelled
% Select for which pixel to run the point-scale version of T&C
% 1: AWS_OnGlacier
% 2: Pluvio
%==========================================================================
sel_forc = 1;  % selection of forcing

%% Study site details
%==========================================================================
%{
Sites
1) "Cinca_Mid"
2) "Tiber"
3) "Velino"

%}
%==========================================================================
SITE = 'Velino';
FORCING = "ERA5Land";
UTM_zone = 33; % for Italy
DeltaGMT= 1; % for Italy

%% Modelling period
%For some reason  "01-Jan-2008 00:00:00" does not work. Only  "01-Jan-2008 01:00:00"
date_start =  "01-Jan-2008 01:00:00"; % Starting point of the simulation
date_end  =  "30-Dec-2008 23:00:00"; % Last timestep of the simulation

%% MODEL PARAMETERS
%study_name = ["Pyrenees_pointscale", "Apennine_pointscale", "Apennine_pointscale", "Apennine_pointscale"];
% Time step for the model
dt=3600; % [s]
dth=1; % [h]

% Integration interval for Solar variables
% Hours or fraction before and after. Values obtained from the
% Automatic_Radiation_Partition_I
t_bef = 1.5; t_aft = -0.5;

%% Data storing  
%==========================================================================
% How to store outputs
% Recommended for long-term simulations (> 10-20 years, otherwise data volume is too important)
%==========================================================================
output_daily = 1; 

%% Load DEM  

%Tmod = 0; % temperature modification above clean ice [°C];
%Pmod = 0; % factor Pmod preciptation modification, e.g. 0.3 means 30% more precipitation at highest elevation
%Z_min = 570; % lowest elevation for linear precipitation modification (min factor -> 0)
%Z_max = 5000; % highest elevation for linear precipitation modification (max factor -> Pmod)

%% Precipitation phase parametrization
% 1 = 2-threshold, 
% 2 = Ding 2017, 
% 3 = single-threshold, 
% 4 = Pomeroy 2013,
% 5 = Wang 2019, 
% 6 = Jennings 2018

parameterize_phase.OPT_Pr_Part = 2; % Choice of the precipitation phase scheme
parameterize_phase.Tmax = 2; % Upper air temperature for dual temperature threshold
parameterize_phase.Tmin = 0; % Lower air temperature for dual temperature threshold
parameterize_phase.Tconst = 2; % Air temperature for constant thresholds

% Skin layer thickness:
hSTL = 0.003; %m

% Albedo scheme choice
Albsno_method = 5; % 3 doesn't work, 4 is Brock 2000, 5 is Ding 2017

%% Create the directory where model outputs will be stored
outlocation = [folder_path,'3_PointScale_version/4_Outputs/'];
if ~exist(outlocation, 'dir'); mkdir(outlocation); addpath(genpath(outlocation)); end

% Saving initial conditions of the model
out = strcat(outlocation,'INIT_COND_', SITE ,'_MultiPoint.mat'); % file path initial conditions

%% Dependencies
addpath(genpath([folder_path,'1_Functions'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path,'5_Common_inputs'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path,'3_PointScale_version/2_Forcing'])); % Where is located the meteorological forcing and Shading matrix 
addpath(genpath([folder_path,'3_PointScale_version/3_Inputs'])); % Add path to Ca_Data

%% Load DEM and geographical information
%==========================================================================
%{
There are different dtm
        1) dtm_Cinca_Mid_250m.mat
        2) dtm_Tiber_250m.mat
        3) dtm_Velino_250m.mat
%}
%==========================================================================
dtm_file = "dtm_Velino_250m.mat"; 
res = 250; % simulation resolution [m]
disp(strcat('Model resolution: ',num2str(res)))

dtm_file_op = strcat(folder_path,'5_Common_inputs/',SITE,'/',dtm_file);
load(dtm_file_op); % Distributed maps pre-processing. Useful here to get the DTM and initial snow depth
DTM = DTM_orig; % Use the full DEM in case running POI outside of mask
DTM(isnan(DTM)) == 0; %%??

[m_cell,n_cell]=size(DTM);

% Precipitation vertical gradient
% Pmod_S = MASK;
% rate = Pmod/(Z_max-Z_min);
% Pmod_S(DTM>Z_min) = 1+rate.*(DTM(DTM>Z_min)-Z_min);

%% Load CO2 data
load('Ca_Data.mat');
Ca_all = Ca;
topo = 1;

%% FAKE SNOW (TRIAL)
%load(strcat(folder_path,'4_Preparation_files/Apennine_preparation_Achille_Method_Distributed/','7_SnowCover/4_Snow_Depth/Apennine_Init_Snow_Depth_virtual.mat'));
%load(strcat(folder_path,'4_Preparation_files/Apennine_preparation_Achille_Method_Distributed/','7_SnowCover/5_Snow_Albedo/Apennine_Init_Snow_Albedo_virtual.mat'));

%% Impose measured albedo on glacier areas
fn_alb_elev = strcat(SITE, '_Albedo_vs_elev.mat');

if exist(fn_alb_elev,'file')>0
disp('Using measured glacier albedo')
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
else 
    Afirn = DTM.*0 + 0.28;   
end

%% Topography for parfor
% load is not possible inside the parfor

Topo_data = load([folder_path,'4_Preparation_files/4_GeoTerrain_MultiPoint/4_Results/1_Velino/',SITE,'_ShF.mat']); % ShF matrix created during pre-processing step

Par_time = Topo_data.Par_time;
Par_points = Topo_data.Par_points;

%DEM
num_cell=numel(DTM);

%MASK = MASK.*0+1;
MASKn=reshape(MASK,num_cell,1);

%% Forcing (for parfor)
% load is not possible inside the parfor
forc_file = strcat(forc_path,'Forcing_ERA5_Land_',SITE,'_2008_corr_all.mat'); % Put here the path of where you downloaded the repository
load(forc_file); % Load forcing table for the current POI

%% Soil parameters (It must be outside the loop)
PSAN=reshape(PSAN,num_cell,1)/100;
PCLA=reshape(PCLA,num_cell,1)/100; 
PORG=reshape(PORG,num_cell,1)/100;

ms=10 ; %% 11 ; 
SOIL_TH=reshape(SOIL_TH,num_cell,1);
ms_max = 10; %% Number of soil layers

%% GLACIERS
%%% DEBRIS
% INIT_COND_v2 has md_max parameter
% Do not delete this code
DEB_MAP=reshape(DEB_MAP,num_cell,1);
md_max = 10; %% % Number of debris layers

%% SNOW
%%% Initial snow depth
SNOWD=reshape(SNOWD,num_cell,1);


% Initial snow albedo
if ~exist('SNOWALB','var')
    SNOWALB = SNOWD;
    SNOWALB(SNOWD>0) = 0.6;
else 
    SNOWALB=reshape(SNOWALB,num_cell,1);
end 

% Debugger
% disp(strcat('Before INIT_COND_V2 caller',num2str(size(Ca,2))))

%% POIs
POI = readtable(strcat(folder_path,'3_PointScale_version/3_Inputs/2_Apennine/Velino_MultiPoints.txt')); %import table with points info
[POI.LAT, POI.LON] = utm2ll(POI.UTM_X, POI.UTM_Y, UTM_zone);

%% Get location for POI
%==========================================================================
%{
yllcorner and xllcorner represent the bottom corner of the original DEM
before flipup it which was set in the preparation file. 
DTM is flipped from the preparation file
Check: imagesc(DTM)
Consider that in matlab the cell (1,1) in the the upper left corner and
the (n,m) cell is the the bottom right corner. 
%}
%==========================================================================


%k = 10;
for k = 1:height(POI)
id_location = char(string(POI.Name(k))); %id

y_coord = POI.UTM_Y(k);
x_coord = POI.UTM_X(k);

pixelX = floor((x_coord - xllcorner) / cellsize) + 1;
pixelY = floor((y_coord - yllcorner) / cellsize) + 1;

%ij = POI.idx(loc);
ij = sub2ind(size(DTM),pixelY, pixelX); % Location
[j, i] = ind2sub(size(DTM), ij); % Location

POI.ij(k) = ij;
POI.i(k) = i;
POI.j(k) = j;

POI.Zbas(k) = DTM(j,i); % Altitude
end

%% categories    [fir     larch    grass  shrub  BLever    BLdec   ]  
zatm_surface = [18      18       2      2      18        18      ]; %Depend on vegetation
zatm_hourly_on = 0;

%% LAND COVER
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

ksv=reshape(VEG_CODE,num_cell,1);

%% LAND COVER PARTITION
%How corine classification behaves
k=98
for k = 1:height(POI)
disp(k)
    switch ksv(POI.ij(k))

    case 1 % Decidious Broad-leaved forest %
        %   Case 1 includes from CORINE: 
        %       1) Broad-leaved forest (32.5%)
        %       2) Mixed forest (1.2%)
        %
        %   From CORINE website:
        %       1) deciduous and evergreen broad-leaved tree species listed under 
        %          the “applicable for” section with >75% cover
        %       2) sporadically occurring <25 ha patches of
        %          shrubs and dwarf shrubs;
        %          herbaceous vegetation (grasses and herbs);
        %          mosses and lichens;
        %          denuded spots.
        %       3) optionally sporadically occurring patches of coniferous trees
        %          not exceeding 25 % share of the tree covered area;
        %       4) palm trees;
        %
        %   From Simone:
        %       1) Broad-leaved forest assumed to be mostly decidious. Evergreen 
        %          broad-leaved are not very common and in the Appennine
        
        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {[0.1 0.1 0.8]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub    BLever    BLdec]
        POI.II(k) =          {[0       0        0      1        1         1    ]>0}; 

 
    case 2 % Grassland/pasture %
        %    Case 2 includes from CORINE:              
        %       1) Pastures (1.3%)       
        %       2) Green urban areas (0.1%)
        %       3) Inland marshes (>0.05%)

        %     It also includes, but not in Tiber basin:
        %       1) Peat bogs
        %       2) Salt marshes
        %       3) Salines
        %       4) Intertaidal flats

        %
        %     From Simone
        %       1) A classification for Grassland/pasture

        POI.Cwat(k) = 0.1; POI.Curb(k) = 0.1; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {[0.7 0.1]};         
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec]
        POI.II(k) =  {[0       0        1      1      0         0    ]>0};  
    

    case 3 % Crops %
        %   Case 3 includes from CORINE:
        %       1)  Non-irrigated arable land (25.4%)
        %       2)  Land principally occupied by agriculture with significant
        %           areas of natural vegetation. (8.8%)
        %       3)  Complex cultivation patterns (7%)
        %       4)  Fruit trees and berry plantations (1.2%)
        %       5)  Permanently irrigated land (0.8%)
        %       6)  Vineyards (0.5%) 
        %       7)  Annual crops associated with permanent crops (0.1%)
        %       8)  Rice fields
        %       9)  Agro-forestry areas (0%)
        %
        %       From Simone: 
        %       1)  Crops (Choose one crop, wheat and sunflowers 
        %           are good choices for the region)
        % 
        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {[0.5 0.5]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec ]
        POI.II(k) =  {[0       0        1      1      0         0     ]>0};  
    

    case 4 % Evergreen needleaves %   
        %    Case 4 includes from CORINE:
        %      1) Coniferous forest (1%)
        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {[1.0]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]
        POI.II(k) =          {[0       0        0      0      1         0      ]>0};  
    


    case 5 % Mediterranean shrublands %
        %    Case 5 includes from CORINE:
        %       1)  Transitional woodland-shrub (5%)
        %       2)  Natural grasslands (4.3%)
        %       3)  Sclerophyllous vegetation (0.4%)
        %       4)  Moors and heathland (>0.05%)
        
        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.1; POI.Cbare(k) = 0.1;
        POI.Ccrown(k) = {[0.6 0.1 0.1]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]
        POI.II(k) =  {[0       0        0      1      1         1      ]>0};  
    

    case 6 % Olives %
        %    Case 6 includes from CORINE:
        %       1) Olive groves (4%)

        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {[1.0]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]
        POI.II(k) =  {[0       0        0      1      0         0      ]>0};  
    

    case 7 % Urban %
        %   Case 7 includes from CORINE:
        %         1) Discontinuous urban fabric (3%)
        %         2) Industrial or commercial units (0.7%)
        %         3) Continuous urban fabric (0.6%)
        %         4) Sport and leisure facilities (0.1%)
        %         5) Road and rail networks and associated land (0.1%)
        %         6) Construction sites (0.1%)
        %         7) Port areas (>0.05%)   
        %         8) Airports (0.1%)

        POI.Cwat(k) = 0; POI.Curb(k) = 0.8 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.1;
        POI.Ccrown(k) = {[0.1]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]
        POI.II(k) =  {[0       0        1      0      0         0      ]>0};  
        

    case 8 % Rock %
        %    Case 8 includes from CORINE:
        %        1) Bare rocks (0.2%)
        %        2) Glaciers and perpetual snow (0%)

        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 1.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {[0.0]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]  
        POI.II(k) =  {[0       0        0      0      0         0      ]>0}; 


    case 9 % Water %
        %    Case 9 includes from CORINE:
        %        1) Water bodies (0.3%)
        %        2) Water courses (0.2%)

        %     It also includes, but not in Tiber basin:
        %       1) Coastal lagoons
        %       2) Estuaries
        %       3) Sea and Ocean


        POI.Cwat(k) = 1.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {[0.0]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]  
        POI.II(k) =  {[0       0        0      0      0         0      ]>0};
  

    case 10 % Bare soils %
        %    Case 10 includes from CORINE:
        %        1) Mineral extraction sites (0.2%)
        %        2) Burnt areas (0.1%)
        %        3) Beaches - dunes - sands (>0.05%)
        %        4) Sparsely vegetated areas (0.7%)
        
        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.9;
        POI.Ccrown(k) = {[0.1]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]  
        POI.II(k) =  {[0       0        0      1      0         0      ]>0};
  

    otherwise
        disp('INDEX FOR SOIL VEGETATION PARAMETER INCONSISTENT')
        %return
end

% Defining zatm in each point
if ~isempty(zatm_surface(cell2mat(POI.II(k))>0))
POI.zatm(k) = max(zatm_surface(cell2mat(POI.II(k))>0)); %choose correct atmospheric reference height
else
POI.zatm(k) = 2;
end



end

%% Initial conditions
cc_max = 3;

INIT_COND_v3(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Ca,SNOWD,SNOWALB,out, cell2mat(POI.II), POI.ij);

load(out);

%e_relN_Htm1(POI.ij(35),:)


%% Creation of variables for parfor
timeDifference = hours(datetime(date_end)-datetime(date_start))+1;
Datam = zeros(timeDifference,4);

%% Fetch time and do date handling
Date = (datetime(date_start):hours(1):datetime(date_end))';
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(1:timeDifference,:) = [YE MO DA HO]; 

%clear YE MO DA HO MI SE

%% Pre allocation
% Pre-allocate output arrays




%% FOR LOOP for locations
loc=1;
parfor loc = 1:height(POI)  

%% Preallocation
Deb_Par = struct( ...
    'alb', [], ...    
    'e_sur', [], ...
    'lan', [], ...
    'rho', [], ...  
    'cs', [], ...
    'zom', []);

Stoich_H = struct( ...
           'Nl', [], ...
           'Ns', [], ...
           'Nr', [], ...
           'Nf', [], ...
           'Nh', [], ...
         'Phol', [], ...
         'Phos', [], ...
         'Phor', [], ...
         'Phof', [], ...
         'Phoh', [], ...
        'Kpotl', [], ...
        'Kpots', [], ...
        'Kpotr', [], ...
        'Kpotf', [], ...
        'Kpoth', [], ...
      'ftransR', [], ...
      'ftransL', [], ...
          'FiS', [], ...
     'Lig_fr_l', [], ...
    'Lig_fr_fr', [], ...
     'Lig_fr_h', [], ...
     'Lig_fr_r', []); 

Stoich_L = struct( ...
           'Nl', [], ...
           'Ns', [], ...
           'Nr', [], ...
           'Nf', [], ...
           'Nh', [], ...
         'Phol', [], ...
         'Phos', [], ...
         'Phor', [], ...
         'Phof', [], ...
         'Phoh', [], ...
        'Kpotl', [], ...
        'Kpots', [], ...
        'Kpotr', [], ...
        'Kpotf', [], ...
        'Kpoth', [], ...
      'ftransR', [], ...
      'ftransL', [], ...
          'FiS', [], ...
     'Lig_fr_l', [], ...
    'Lig_fr_fr', [], ...
     'Lig_fr_h', [], ...
     'Lig_fr_r', []); 

ParEx_H = struct( ...
    'dexmy', [], ...
    'OPT_EXU', [], ...
       'bfix', [] , ...
        'kc1', [], ...
        'kc2', [], ...
        'kc3', [] , ...
        'kn1', [], ...
        'kn2', [],...
        'kn3', [], ...
        'kp1', [], ...
        'kp2', [], ...
        'kp3', [], ....
        'kk1', [], ...
        'kk2', [], ...
        'kk3', []);

ParEx_L = struct( ...
    'dexmy', [], ...
    'OPT_EXU', [], ...
       'bfix', [] , ...
        'kc1', [], ...
        'kc2', [], ...
        'kc3', [] , ...
        'kn1', [], ...
        'kn2', [],...
        'kn3', [], ...
        'kp1', [], ...
        'kp2', [], ...
        'kp3', [], ....
        'kk1', [], ...
        'kk2', [], ...
        'kk3', []);

    Mpar_H = struct( ...  
           'jDay_cut', [], ...
            'LAI_cut', [], ...
          'jDay_harv', [], ...
             'B_harv', [], ...
           'Date_log', [], ...
          'fract_log', [], ...
          'Date_fire', [], ...
           'fire_eff', [], ...
           'funb_nit', [], ...
      'Date_girdling', [], ...
     'fract_girdling', [], ...
        'Date_sowing', [], ...
    'Date_harvesting', [], ...
             'Crop_B', [], ...
         'Crop_crown', [], ...
     'fract_resprout', [], ...
         'fract_left', [], ... 
      'fract_left_fr', [], ...
      'fract_left_AB', [], ...
      'fract_left_BG', []);


    Mpar_L = struct( ...  
           'jDay_cut', [], ...
            'LAI_cut', [], ...
          'jDay_harv', [], ...
             'B_harv', [], ...
           'Date_log', [], ...
          'fract_log', [], ...
          'Date_fire', [], ...
           'fire_eff', [], ...
           'funb_nit', [], ...
      'Date_girdling', [], ...
     'fract_girdling', [], ...
        'Date_sowing', [], ...
    'Date_harvesting', [], ...
             'Crop_B', [], ...
         'Crop_crown', [], ...
     'fract_resprout', [], ...
         'fract_left', [], ... 
      'fract_left_fr', [], ...
      'fract_left_AB', [], ...
      'fract_left_BG', []);

% Or appropriate size and data type

PFT_opt_H = struct( ...
    'chiL', [], ...
    'alf_lf_vis', [], ...
    'alf_lf_nir', [], ...
    'alf_st_vis', [], ...
    'alf_st_nir',  [], ...
    'alf_ld_vis', [], ...
    'alf_ld_nir', [], ...
    'tau_lf_vis', [], ...
    'tau_lf_nir', [], ...
    'tau_st_vis', [], ...
    'tau_st_nir', [], ...
    'tau_ld_vis', [], ...
    'tau_ld_nir', [])

PFT_opt_L = struct( ...
    'chiL', [], ...
    'alf_lf_vis', [], ...
    'alf_lf_nir', [], ...
    'alf_st_vis', [], ...
    'alf_st_nir',  [], ...
    'alf_ld_vis', [], ...
    'alf_ld_nir', [], ...
    'tau_lf_vis', [], ...
    'tau_lf_nir', [], ...
    'tau_st_vis', [], ...
    'tau_st_nir', [], ...
    'tau_ld_vis', [], ...
    'tau_ld_nir', [])




%% Crowns
cc_aux = POI.NCrown(loc); 
II = cell2mat(POI.II(loc));
Cwat = POI.Cwat(loc); 
Curb = POI.Curb(loc); 
Crock = POI.Crock(loc);
Cbare = POI.Cbare(loc);
Ccrown = cell2mat(POI.Ccrown(loc));
zatm = POI.zatm(loc);
id_location = string(POI.Name(loc));

%% Locations
Zbas = POI.Zbas(loc)
Lat = POI.LAT(loc);
Lon = POI.LON(loc);
ij = POI.ij(loc);

%% FORCING
%==========================================================================
% ERA5
%==========================================================================

fieldNames = fieldnames(forc_f);  % Get all field names as cell array
forc = forc_f.(string(POI.Name(loc))); % Table of 

Date_all=forc.Date; 

%define period and time zone info
x1=find(date_start == Date_all,1);
x2=find(date_end == Date_all,1);


%% Displaying modelling parameters
disp(strcat("Site selected: ", SITE))
disp(['Forcing selected: ' char(FORCING)])
disp(['Running T&C for pixel: ' char(id_location)])
disp(['Simulation period: ' datestr(date_start) ' to ' datestr(date_end)])



%% Carbon data
%==========================================================================
% Load carbon data and narrow down to period
%==========================================================================

formattedDate_CO2 = datetime(Date_CO2, 'ConvertFrom', 'datenum');

d1 = find(formattedDate_CO2 == Date(1));
d2 = find(formattedDate_CO2 == Date(end));

%d1 = find(abs(Date_CO2-datenum(Date(1)))<1/36);
%d2 = find(abs(Date_CO2-datenum(Date(end)))<1/36);
Ca=Ca_all(d1:d2); 
%clear d1 d2

Oa= 210000; % Intercellular Partial Pressure Oxygen [umolO2/mol]

% Debugging
% formattedDate = datetime(Date_CO2, 'ConvertFrom', 'datenum');
% disp(strcat('Launcher',num2str(size(Ca,2))))


%% Meteorological input
%==========================================================================
%{
FORCINGS: 
    Precipitation [mm]
    Air Pressure [Pa]
    Temperature [C]
    Wind Speed [m s-1]    
    Radiation [W m-2]    
    Dew Point temperature [C]

CALCULATIONS:
    Relative Humidity [-]

NOTE:
Forcings must be put in the model as double. If Ta is as single then it
crashes. 
%}
%==========================================================================

% Period of forcing data
forcing = forc(x1:x2,:);
NN= height(forcing);%%% time Step

% Height of virtual station
zatm_hourly = repmat(2.00,height(forcing),1);



%% Precipitation
% Precipitation from ERA5Land
Pr=double(forcing.Total_Precipitation_HH);
Pr(isnan(Pr))=0;
Pr(Pr<0.01)=0;

%% Air pressure
% Air Pressure in mbar based on the PDF for variables and parameters of TC
% Pressure comes in Pa from ERA5
Pre=double(forcing.Pressure/100);    

%% Temperature
% 2m air temperature
Ta=double(forcing.Temperature);

%% Wind Speed
Ws=double(forcing.Wind_Speed);
Ws(Ws < 0.01) = 0.01;

%% Relative humidity
% Divided by 100 to set the number in the range 0-1
U=double(forcing.RH/100);

%% Longwave radiation
% N can be cloud cover [-] or longwave incoming radiation [W m-2]
% Latm=forcing.LW_rad_downward_HH; % Latm:Incoming long wave radiation [W m-2]
% N=ones(NN,1); % cloud cover [-]
% N = forcing.LW_rad_downward_HH;
N = forcing.N;

%% Radiation partition
SAD1=double(forcing.SAD1); SAD2=double(forcing.SAD2); 
SAB1=double(forcing.SAB1); SAB2=double(forcing.SAB2);
PARB=double(forcing.PARB); PARD=double(forcing.PARD);

%% Albedo parameters
%Ameas = ones(NN,1);
alpha = 0; % switch for albedo
%Ameas_t=0; % albedo
%Aice_meas_on_hourly = ones(height(forcing),1)/2; % albedo
%Asno_meas_on_hourly = ones(height(forcing),1)/2; % albedo

%% Vapor pressure - Dew Point temperature
%esat/ea/Ds/Tdew
esat=double(forc.es);   % Vapour pressure at saturation (Pa)
ea=double(forc.ea);     % Vapour pressure (Pa)
Ds= esat - ea;  % Vapor Pressure Deficit (Pa)
Ds(Ds<0)=0; 
Tdew= double(forc.Dew_Point_Temp);

%a=17.27; b=237.3;
%clear a b xr;
%xr=a*Ta./(b+Ta)+log10(U);
%Tdew=b*xr./(a-xr);          % Presumed dewpoint temperature (°C)

%% DING PARAMETERIZATION
% Initial daily mean values for Ding parametrization
Ta_Ding_d = nanmean(Ta(1:24));
Pre_Ding_d = nanmean(Pre(1:24));
ea_Ding_d = nanmean(ea(1:24));

%% RADIATION AND TOPOGRAPHY

if topo == 1
    
    %% Topography

    Par_time_table=Par_time.(string(POI.Name(loc))); %Table to use
    Par_time_period = Par_time_table(Par_time_table.Time>=date_start & Par_time_table.Time<=date_end, :); %variables per period

    %ShF
    ShF = Par_time_period.ShF_S;
    
    rho_g = 0.35; %%% Spatial Albedo
    zeta_S = Par_time_period.zeta_Sts;
    h_S = Par_time_period.h_Sts;
    %clear Par_time

    Par_points_data = Par_points(Par_points.Name == string(POI.Name(loc)),:);
    SvF = Par_points_data.SvF_S; % Sky view factor at pixel ij
    Ct = Par_points_data.Ct_S;
    Slo_top = Par_points_data.Slo_top_S; % Slope at pixel ij
    Aspect = Par_points_data.Aspect; % Aspect at pixel ij
    %clear Points

    cos_fst = cos(atan(Slo_top)).*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;

    %clear zeta_S 
    
    %% Radiative parameters
    SAB1(sin(h_S) <= 0.10)  =  0;
    SAB2(sin(h_S) <= 0.10)  =  0;
    SAD1(sin(h_S) <= 0.10)  =  0;
    SAD2(sin(h_S) <= 0.10)  =  0;
    PARB(sin(h_S) <= 0.10)  =  0;
    PARD(sin(h_S) <= 0.10)  =  0;

    SAD1 = SAD1.*SvF + Ct.*rho_g.*((SAB1./sin(h_S)).*cos_fst + (1-SvF).*SAD1);
    SAD2 = SAD2.*SvF + Ct.*rho_g.*((SAB2./sin(h_S)).*cos_fst + (1-SvF).*SAD2);
    PARD = PARD.*SvF + Ct.*rho_g.*((PARB./sin(h_S)).*cos_fst + (1-SvF).*PARD);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAB1 =(SAB1./sin(h_S)).*cos_fst.*ShF;
    SAB2 =(SAB2./sin(h_S)).*cos_fst.*ShF;
    PARB = (PARB./sin(h_S)).*cos_fst.*ShF;

    %correSITEions, temporary
    SAB1(SAB1<0)=0; SAB2(SAB2<0)=0;
    PARB(PARB<0)=0; PARD(PARD<0)=0;
    SAD1(SAD1<0)=0; SAD2(SAD2<0)=0;

    SAB1(isnan(SAB1))= 0; SAB2(isnan(SAB2)) = 0;
    SAD1(isnan(SAD1))= 0; SAD2(isnan(SAD2)) = 0;
    PARB(isnan(PARB))= 0; PARD(isnan(PARD)) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   SvF=1; %% Sky View factor = 1 if topography is not considered
   Slo_top_S = MASKn*0;
   Slo_top=Slo_top_S(ij);
end


%% SOIL 
%==========================================================================
% Vector is divided by 100 to put the numbers within 0-1
% Original PSAN, PCLA and PORG come with values between 0-100, g/100g
%==========================================================================
Psan = PSAN(ij); % Soil sand content at pixel ij
Pcla = PCLA(ij); % Soil clay content at pixel ij
Porg= PORG(ij); % Soil organic content at pixel ij

%% SAVING INITIAL CONDITIONS AND PARAMETERS 
%==========================================================================
% In MultiPoint analysis INIT_COND_Tiber_MultiPoint.mat is created.
% This can cause problems with the initial conditions. Specially with Ca.
% INIT_COND_v2 depends on the Land Cover. Here it changes the initial
% condition 
%==========================================================================
% (run this only once in MultiPoint model!)
%if exist(out, 'file') == 2
%load(out);
%else



%% RUN MODEL
%==========================================================================
% PARAM_IC: Define parameter file
% MAIN_FRAME: Contains the model
%==========================================================================

PARAM_IC = strcat(folder_path,'3_PointScale_version/3_Inputs/MOD_PARAM_Multipoint.m');


%% MAIN_FRAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% MAIN_FRAME OF HBM-VEG %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INITIZIALIZATION VARIABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Debugger
%disp(strcat('Initial MAIN_FRAME:',num2str(size(Ca,2))))

%%% j time dt = 1 day 
%%% i time dt = 1h
%%% ms soil layer
%%% cc Crown Area number
%%% NN time step
dtd = 1; %%Time step [day]
dth = dt/3600; %% Time step [hour]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=zeros(NN,ms); 
O=zeros(NN,ms);
Qi_out=zeros(NN,ms); % Outgoing Lateral subsurface flow
WTR=zeros(NN,ms); % Water flow due to water table rising
POT=zeros(NN,ms); % Soil water potential
Tdp=zeros(NN,ms); % Soil Temperature of the layer
%Sdp=zeros(NN,ms);
Vice=zeros(NN,ms);
Oice=zeros(NN,ms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
OF=zeros(NN,1); OS=zeros(NN,1);
ZWT=zeros(NN,1);
Pr_sno=zeros(NN,1); Pr_liq=zeros(NN,1);
er=zeros(NN,1);
Ts=zeros(NN,1);ra=zeros(NN,1);r_soil=zeros(NN,1);
b_soil=zeros(NN,1); Tice=zeros(NN,1);
Rn=zeros(NN,1);H=zeros(NN,1);
QE=zeros(NN,1);Qv=zeros(NN,1); Lpho=zeros(NN,1);
EG=zeros(NN,1); ESN=zeros(NN,1);ESN_In=zeros(NN,1);
EWAT=zeros(NN,1);EIn_urb=zeros(NN,1);
EIn_rock=zeros(NN,1); dw_SNO=zeros(NN,1);
G=zeros(NN,1); SWE=zeros(NN,1);
SND=zeros(NN,1);
%snow_alb=ones(NN,1);
ros=zeros(NN,1);In_SWE=zeros(NN,1);SP_wc=zeros(NN,1);
WR_SP=zeros(NN,1);U_SWE=zeros(NN,1);NIn_SWE=zeros(NN,1);
dQ=zeros(NN,1);Qfm=zeros(NN,1);t_sls=zeros(NN,1);
DQ=zeros(NN,1);DT=zeros(NN,1);
In_urb=zeros(NN,1);In_rock=zeros(NN,1);
SE_rock=zeros(NN,1);SE_urb=zeros(NN,1);
Csno=zeros(NN,1);WIS=zeros(NN,1);
Lk=zeros(NN,1);f=zeros(NN,1);
Rh=zeros(NN,1);Rd=zeros(NN,1);
NDVI=zeros(NN,1);alp_soil=zeros(NN,1);
tau_sno=zeros(NN,1);e_sno=zeros(NN,1);ALB=zeros(NN,1);
EK=zeros(NN,1);
dQVEG=zeros(NN,1);TsV=zeros(NN,1);
HV=zeros(NN,1);QEV=zeros(NN,1);
Cice=zeros(NN,1);
Lk_wat=zeros(NN,1); Lk_rock=zeros(NN,1);
EICE=zeros(NN,1);
WAT=zeros(NN,1);ICE=zeros(NN,1);ICE_D=zeros(NN,1);
WR_IP=zeros(NN,1);NIce=zeros(NN,1);Cicew=zeros(NN,1);
IP_wc=zeros(NN,1);
Csnow=zeros(NN,1);FROCK=zeros(NN,1);
Imelt=zeros(NN,1);Smelt=zeros(NN,1);
Tdamp=zeros(NN,1); Gfin=zeros(NN,1);
In_Litter=zeros(NN,1); ELitter=zeros(NN,1);
Ws_under=zeros(NN,1);
Ts_under=NaN*ones(NN,1);
Tdpsnow = zeros(NN,5);



IrD=zeros(NN,1);
Salt=zeros(NN,1); %%% Salt = Salt Concentration [mol Salt/ m-3] 

NetWatWet=zeros(NN,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_litter=zeros(NN,cc_aux);
%%%
In_H=zeros(NN,cc_aux);In_L=zeros(NN,cc_aux);
OH=zeros(NN,cc_aux);OL=zeros(NN,cc_aux);
Psi_s_H=zeros(NN,cc_aux); Psi_s_L=zeros(NN,cc_aux);
Tdp_H=zeros(NN,cc_aux);Tdp_L=zeros(NN,cc_aux);
rap_H=zeros(NN,cc_aux);rap_L=zeros(NN,cc_aux);
Dr_H=zeros(NN,cc_aux);Dr_L=zeros(NN,cc_aux);
T_H=zeros(NN,cc_aux);T_L=zeros(NN,cc_aux);EIn_H=zeros(NN,cc_aux);EIn_L=zeros(NN,cc_aux);
rb_H=zeros(NN,cc_aux);rb_L=zeros(NN,cc_aux);
An_H=zeros(NN,cc_aux);Rdark_H=zeros(NN,cc_aux);
An_L=zeros(NN,cc_aux);Rdark_L=zeros(NN,cc_aux);
Ci_sunH=zeros(NN,cc_aux);Ci_sunL=zeros(NN,cc_aux);
Ci_shdH=zeros(NN,cc_aux);Ci_shdL=zeros(NN,cc_aux);
rs_sunH=zeros(NN,cc_aux);rs_sunL=zeros(NN,cc_aux);
rs_shdH=zeros(NN,cc_aux);rs_shdL=zeros(NN,cc_aux);
%%%%%
Vx_H=zeros(NN,cc_aux);Vl_H=zeros(NN,cc_aux);
Vx_L=zeros(NN,cc_aux);Vl_L=zeros(NN,cc_aux);
Psi_x_H=zeros(NN,cc_aux);Psi_l_H=zeros(NN,cc_aux);
Psi_x_L=zeros(NN,cc_aux);Psi_l_L=zeros(NN,cc_aux);
gsr_H=zeros(NN,cc_aux);
Jsx_H=zeros(NN,cc_aux); Jxl_H=zeros(NN,cc_aux);
Kleaf_H=zeros(NN,cc_aux);Kx_H=zeros(NN,cc_aux);
gsr_L=zeros(NN,cc_aux);
Jsx_L=zeros(NN,cc_aux);Jxl_L=zeros(NN,cc_aux);
Kleaf_L=zeros(NN,cc_aux);Kx_L=zeros(NN,cc_aux);
fapar_H=zeros(NN,cc_aux);fapar_L=zeros(NN,cc_aux);
SIF_H=zeros(NN,cc_aux);SIF_L=zeros(NN,cc_aux);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Variables for Carbon processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NCP = 8;               % Number of Carbon Pool
NNd = ceil(NN/24)+1;   % Daily time step number

LAI_L=zeros(NNd,cc_aux);   % Leaf area index - low vegetation
B_L=zeros(NNd,cc_aux,NCP); % Carbon pool biomass
NPP_L=zeros(NNd,cc_aux);   % Net Primary Production - low vegetation
Rg_L=zeros(NNd,cc_aux);    % Growth respiration - low vegetation
RA_L=zeros(NNd,cc_aux);    % Autotrophic respiration - low vegetation
Rmc_L=zeros(NNd,cc_aux);   % Maitenance respiration - low vegetation
ANPP_L=zeros(NNd,cc_aux);  % Above ground net primary production - low vegetation
Rms_L=zeros(NNd,cc_aux);   % Maitenance respiration sapwood - low vegetation
Rmr_L=zeros(NNd,cc_aux);   % Maitenance respiration roots - low vegetation
PHE_S_L=zeros(NNd,cc_aux); % Phenology state - low vegetation 
dflo_L=zeros(NNd,cc_aux);  % Days from leaf onset - low vegetation
AgeL_L=zeros(NNd,cc_aux);  % Leaf Age - low vegetation
e_rel_L=ones(NNd,cc_aux);  % Relative efficiency of the photosynthesis apparatus
SAI_L=zeros(NNd,cc_aux);   % Steam area index - low vegetation
hc_L=zeros(NNd,cc_aux);    % Vegetation height - low vegetation
e_relN_L=ones(NNd,cc_aux); % Relative Efficiency of the photosynthesis apparatus due to N limitations - low vegetation
BA_L=zeros(NNd,cc_aux);    %

%%%
LAI_H=zeros(NNd,cc_aux); B_H=zeros(NNd,cc_aux,NCP);NPP_H=zeros(NNd,cc_aux);Rg_H=zeros(NNd,cc_aux);
RA_H=zeros(NNd,cc_aux);Rms_H=zeros(NNd,cc_aux);Rmr_H=zeros(NNd,cc_aux); ANPP_H=zeros(NNd,cc_aux);
PHE_S_H=zeros(NNd,cc_aux); dflo_H=zeros(NNd,cc_aux);Rmc_H=zeros(NNd,cc_aux);
AgeL_H=zeros(NNd,cc_aux);e_rel_H=ones(NNd,cc_aux);SAI_H=zeros(NNd,cc_aux); hc_H=zeros(NNd,cc_aux); e_relN_H=ones(NNd,cc_aux); BA_H=zeros(NNd,cc_aux);
%%%
Rrootl_H=zeros(NNd,cc_aux); Rrootl_L=zeros(NNd,cc_aux);
Bfac_dayH=ones(NNd,cc_aux); Bfac_weekH=ones(NNd,cc_aux); NPPI_H=zeros(NNd,cc_aux); TdpI_H=zeros(NNd,cc_aux);
Bfac_dayL=ones(NNd,cc_aux); Bfac_weekL=ones(NNd,cc_aux); NPPI_L=zeros(NNd,cc_aux); TdpI_L=zeros(NNd,cc_aux);
%%%
Sr_H=zeros(NNd,cc_aux);  Slf_H=zeros(NNd,cc_aux);
Sfr_H=zeros(NNd,cc_aux); Swm_H=zeros(NNd,cc_aux);  Sll_H=zeros(NNd,cc_aux);
Sr_L=zeros(NNd,cc_aux); Slf_L=zeros(NNd,cc_aux);
Sfr_L=zeros(NNd,cc_aux); Swm_L=zeros(NNd,cc_aux); Sll_L=zeros(NNd,cc_aux);
LAIdead_H=zeros(NNd,cc_aux); LAIdead_L=zeros(NNd,cc_aux);
Rexmy_H=zeros(NNd,cc_aux,3); AgeDL_H=zeros(NNd,cc_aux); Nreserve_H=zeros(NNd,cc_aux);
Preserve_H=zeros(NNd,cc_aux); Kreserve_H=zeros(NNd,cc_aux);
rNc_H=zeros(NNd,cc_aux);rPc_H=zeros(NNd,cc_aux);rKc_H=zeros(NNd,cc_aux);
Rexmy_L=zeros(NNd,cc_aux,3); AgeDL_L=zeros(NNd,cc_aux); Nreserve_L=zeros(NNd,cc_aux);
Preserve_L=zeros(NNd,cc_aux); Kreserve_L=zeros(NNd,cc_aux);
rNc_L=zeros(NNd,cc_aux);rPc_L=zeros(NNd,cc_aux);rKc_L=zeros(NNd,cc_aux);
NBLeaf_H =zeros(NNd,cc_aux); PARI_H=zeros(NNd,cc_aux,3) ; NBLI_H=zeros(NNd,cc_aux);
NBLeaf_L =zeros(NNd,cc_aux);PARI_L=zeros(NNd,cc_aux,3) ; NBLI_L=zeros(NNd,cc_aux);
%%%%%
NupI_H=zeros(NNd,cc_aux,3);
NupI_L=zeros(NNd,cc_aux,3);
NavlI=zeros(NNd,3);
%%%
RB_L=zeros(NNd,cc_aux,7);
RB_H=zeros(NNd,cc_aux,7);
NuLit_H =zeros(NNd,cc_aux,3);
NuLit_L =zeros(NNd,cc_aux,3);
%%%%%
BLit=zeros(NNd,cc_aux);
%%%%%
Nuptake_H=zeros(NNd,cc_aux);
Puptake_H=zeros(NNd,cc_aux);
Kuptake_H=zeros(NNd,cc_aux);
FNC_H=zeros(NNd,cc_aux);
TexC_H=zeros(NNd,cc_aux);TexN_H=zeros(NNd,cc_aux);TexP_H=zeros(NNd,cc_aux);TexK_H=zeros(NNd,cc_aux);
TNIT_H=zeros(NNd,cc_aux);TPHO_H=zeros(NNd,cc_aux);TPOT_H=zeros(NNd,cc_aux);
SupN_H=zeros(NNd,cc_aux);SupP_H=zeros(NNd,cc_aux);SupK_H=zeros(NNd,cc_aux);
%%%%
Nuptake_L=zeros(NNd,cc_aux);
Puptake_L=zeros(NNd,cc_aux);
Kuptake_L=zeros(NNd,cc_aux);
FNC_L=zeros(NNd,cc_aux);
TexC_L=zeros(NNd,cc_aux);TexN_L=zeros(NNd,cc_aux);TexP_L=zeros(NNd,cc_aux);TexK_L=zeros(NNd,cc_aux);
TNIT_L=zeros(NNd,cc_aux);TPHO_L=zeros(NNd,cc_aux);TPOT_L=zeros(NNd,cc_aux);
SupN_L=zeros(NNd,cc_aux);SupP_L=zeros(NNd,cc_aux);SupK_L=zeros(NNd,cc_aux);
%%%%
ISOIL_H=zeros(NNd,cc_aux,18);
ISOIL_L=zeros(NNd,cc_aux,18);
ManIH = zeros(cc_aux,1);
ManIL = zeros(cc_aux,1);



AgrHarNut =  zeros(NNd,3);

%%%%%%%%%%%%%%%%%

jDay=zeros(NNd,1);L_day=zeros(NNd,1);
%%%%
Ccrown_t =ones(NNd,cc_aux);
AgePl_H =zeros(NNd,cc_aux); AgePl_L =zeros(NNd,cc_aux);
Tden_H =zeros(NNd,cc_aux); Tden_L =zeros(NNd,cc_aux);
TBio_Ht =zeros(NNd,cc_aux); TBio_Lt =zeros(NNd,cc_aux);
ZR95_Ht =zeros(NNd,cc_aux); ZR95_Lt =zeros(NNd,cc_aux);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OPT_SoilBiogeochemistry = 0;

%{
if  not(exist('OPT_BG','var'))
    OPT_SoilBiogeochemistry = 0;
else
    if  OPT_BG == 0
        OPT_SoilBiogeochemistry = 0;
    else
        OPT_SoilBiogeochemistry = 1;
    end 
end
%}
%%%%%%%%%%%%%
if OPT_SoilBiogeochemistry == 1
    P=zeros(NNd,55);
    R_litter=zeros(NNd,1);
    R_microbe=zeros(NNd,1);
    R_litter_sur=zeros(NNd,1);
    R_ew=zeros(NNd,1);
    VOL=zeros(NNd,1);
    N2flx = zeros(NNd,1);
    Min_N = zeros(NNd,1);
    Min_P = zeros(NNd,1);
    R_bacteria= zeros(NNd,1);
    RmycAM = zeros(NNd,1);
    RmycEM = zeros(NNd,1);
    Prod_B = zeros(NNd,1);
    Prod_F = zeros(NNd,1);
    BfixN = zeros(NNd,1);
    LitFirEmi =  zeros(NNd,2);
    LEAK_NH4 = zeros(NNd,1);
    LEAK_NO3 = zeros(NNd,1);
    LEAK_P = zeros(NNd,1);
    LEAK_K = zeros(NNd,1);
    LEAK_DOC = zeros(NNd,1);
    LEAK_DON = zeros(NNd,1);
    LEAK_DOP = zeros(NNd,1);
    R_NH4 = zeros(NNd,1);
    R_NO3 = zeros(NNd,1);
    R_P = zeros(NNd,1);
    R_K = zeros(NNd,1);
    R_DOC = zeros(NNd,1);
    R_DON = zeros(NNd,1);
    R_DOP = zeros(NNd,1);
    RexmyI=zeros(NNd,3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CALL PARAMETERS AND INITIAL CONDITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%categories  [fir    larch     grass    shrub  BLever   BLdec  ]  
Kbot    =    [0.0    0.0       0.0      0.0    0.0      0.0    ]; % Conductivity at the bedrock layer [mm/h] 
Krock   =    [NaN    NaN       NaN      NaN    NaN      NaN    ]; % Conductivity of Fractured Rock [mm/h] 

Kbot= Kbot(ksv(ij)); % Conductivity of the bedrock [mm/h] 
Krock=Krock(ksv(ij)); % Hydraulic conductivity fractured rock [mm/h]

% Soil layers
Zs= [0    10    20    50   100   150   200   300   400    700   1000]; % Depth of top of the soil layer [mm],  ms+1
Zdes = 10; % Depth of evaporation layer [mm]
Zinf=  10; % Depth of infiltration layer (=first layer) [mm]
Zbio = 250; % Depth of the active Biogeochemistry zone [mm]

if  not(length(Zs)==ms+1)
    disp('SOIL LAYER MESH INCONSISTENT')
    %return
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
ix = [];


%% Soil parameter function

% Matric Potential
Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]

% Function
[Ofc,Oss_Lp,Owp_Lp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,1,1,Phy);
%clear Oss_Lp Owp_Lp
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

%if exist('Afirn','var')
Aice=Afirn(ij); %Use landsat distributed measured albedo
%end 

%% Debris Cover Glacier
%==========================================================================
% Debris Cover Glacier
% MAKE SURE DEBRIS THICKNESS 0 FOR CLEAN GLACIER!!
%==========================================================================

albs   =     [0.153  0.115     0.13     0.13   0.13     0.13 ]; % Ask Achille??
lans   =     [1.65   0.985     1.45     0.94   0.94     0.94 ];
zoms   =     [0.38   0.081     0.15     0.016  0.016    0.016];

Deb_Par.alb= albs(sel_basin);
Deb_Par.e_sur =  0.94;
Deb_Par.lan = lans(sel_basin);
Deb_Par.rho = 1496;  % [kg/m^3]
Deb_Par.cs = 948;   % [J/kg K]
Deb_Par.zom = zoms(sel_basin);

dbThick=DEB_MAP(ij);%% [mm]

%%  ROOT PARAMETER 

%BLever: broadleaf evergreen vegetation
%BLdec: broadleaf deciduous vegetation

%categories  [fir    larch     grass    shrub  BLever   BLdec  ]  
ExEM   =     [0.0    0.0       0.0      0.0    0.0      0.0    ]; % Fraction of EM mycorrhizal on the toal mycorrhizal [-]
% Selection based on II
ExEM = ExEM(1); %%CHECK! FOR simplicity

CASE_ROOT= 1;  % Type of Root Profile
%categories  [fir    larch     grass    shrub  BLever   BLdec ]  
ZR95_H   =   [800    800       0        0      1000     800   ]; % Root depth 95 percentile, high vegetation [mm]
ZR95_L   =   [0      0         200      600    0        0     ]; % Root depth 95 percentile, low vegetation [mm]
ZR50_H   =   [NaN    NaN       NaN      NaN    NaN      NaN   ]; % Root depth 50 percentile, high vegetation [mm]
ZR50_L   =   [NaN    NaN       NaN      NaN    NaN      NaN   ]; % Root depth 50 percentile, low vegetation [mm]
ZRmax_H  =   [NaN    NaN       NaN      NaN    NaN      NaN   ]; % Maximum root depth, high vegetation [mm]
ZRmax_L  =   [NaN    NaN       NaN      NaN    NaN      NaN   ]; % Maximum root depth, low vegetation [mm]

% Selection based on II
ZR95_H =ZR95_H(II); ZR50_H =ZR50_H(II); ZRmax_H =ZRmax_H(II);
ZR95_L =ZR95_L(II); ZR50_L =ZR50_L(II); ZRmax_L =ZRmax_L(II);

%debugger
%disp(ZRmax_H)


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

In_max_urb = 5;    % Maximum interception capacity in urban [mm]
In_max_rock = 0.1; % Maximum interception capacity in rocks [mm]
Kct=0.75;         % Foliage cover decay factor for throughfall [-]

%% Interception Parameter
gcI=3.7;       % Interception parameter [1/mm]
KcI=0.06;      % Interception drainage rate coefficient [mm] - Mahfouf and Jacquemin 1989
Sp_SN_In= 5.9; % Specific interception of rainfall for unit leaf area. Average of high vegetation [mm/LAI]

%categories  [fir    larch     grass    shrub  BLever   BLdec ]  
Sp_LAI_H_In= [0.1     0.1      0.2      0.2    0.2      0.2   ]; % Specific interception of rainfall for unit leaf area, high vegetation [mm/LAI]
Sp_LAI_L_In= [0.2     0.2      0.2      0.1    0.2      0.2   ]; % Specific interception of rainfall for unit leaf area, low vegetation [mm/LAI]

% Selection based on II
Sp_LAI_H_In =Sp_LAI_H_In(II);
Sp_LAI_L_In =Sp_LAI_L_In(II);

%% Leaf Dimension

%categories  [fir    larch     grass    shrub  BLever   BLdec  ]  
d_leaf_H =   [0.25   0.8       2        2      5        4      ]; % Leaf characteristic dimension, high vegetation [cm]
d_leaf_L =   [2      2         0.8      3      2        2      ]; % Leaf characteristic dimension, low vegetation [cm]

% Selection based on II
d_leaf_H = d_leaf_H(II);
d_leaf_L =d_leaf_L(II);

%% Veg Biochemical parameter

%categories  [fir    larch     grass    shrub  BLever   BLdec  ]  
KnitH      = [0.35   0.2       NaN      NaN    0.35     0.30   ]; % Canopy Nitrogen Decay coefficient, high vegetation [-]
KnitL      = [NaN    NaN       0.15     0.25   NaN      NaN    ]; % Canopy Nitrogen Decay coefficient, low vegetation [-]

mSl_H      = [0      0         NaN      NaN    0        0      ]; % Linear coefficient of increasing specific leaf area with LAI, high vegetation [m2 PFT /gC] 
mSl_L      = [NaN    NaN       0        0      NaN      NaN    ]; % Linear coefficient of increasing specific leaf area with LAI, low vegetation  [m2 LAI/gC * m2 PFT / m2 LAI]  0.0 - 0.004  Brod. Dec. Tree

% Selection based on II
KnitH =KnitH(II); 
KnitL =KnitL(II);
mSl_H =mSl_H(II); 
mSl_L =mSl_L(II);

%%  Photosynthesis Parameter

% High vegetation
%categories  [fir    larch     grass    shrub  BLever   BLdec ] 
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
%categories  [fir    larch     grass    shrub  BLever   BLdec ] 
FI_L   =     [NaN    NaN       0.081    0.081  NaN      NaN   ]; % Intrinsec quantum Efficiency [umolCO2/umolPhotons]
Do_L   =     [NaN    NaN       1000     1000   NaN      NaN   ]; % [Pa] 
a1_L   =     [NaN    NaN       8        7      NaN      NaN   ]; % WUE parameter [-]
go_L   =     [NaN    NaN       0.01     0.01   NaN      NaN   ]; % Minimum Stomatal Conductance [mol / s m^2] 
CT_L   =     [NaN    NaN       3        3      NaN      NaN   ]; % Photosyntesis Typology for Plants 'CT' == 3  'CT' ==  4  
DSE_L  =     [NaN    NaN       0.649    0.649  NaN      NaN   ]; % Activation Energy - Plant Dependent [kJ/mol] 
Ha_L   =     [NaN    NaN       62       72     NaN      NaN   ]; % Entropy factor - Plant Dependent [kJ / mol K]  
gmes_L =     [NaN    NaN       Inf      Inf    NaN      NaN   ]; % Mesophyll conductance [mol CO2 / s m^2 ];  
rjv_L  =     [NaN    NaN       1.9      2.2    NaN      NaN   ]; % Ratio Jmax - Vmax  [umol electrons / umolCO2 ]

% Selection based on II
FI_L=FI_L(II); Do_L=Do_L(II); a1_L=a1_L(II); go_L=go_L(II);
CT_L=CT_L(II); DSE_L=DSE_L(II); Ha_L=Ha_L(II); gmes_L=gmes_L(II);
rjv_L=rjv_L(II);

%categories  [fir    larch     grass    shrub  BLever   BLdec  ] 
Vmax_H  =    [45     64        0        0      32       48     ]; % Maximum Rubisco Capacity [umol CO2 /m2 s]  
Vmax_L  =    [0      0         50       46     0        0      ]; % Maximum Rubisco Capacity [umol CO2 /m2 s]

% Selection based on II
Vmax_H  =   Vmax_H(II); Vmax_L =Vmax_L(II);


%% Hydraulic Parameters
%categories    [fir    larch   grass    shrub  BLever   BLdec   ] 
Psi_sto_00_H = [-0.8   -0.8    NaN      NaN    -1.0     -0.8    ]; % Water Potential at 2% loss conductivity [MPa] 
Psi_sto_50_H = [-2.5   -2.5    NaN      NaN    -2.8     -2.5    ]; % Water Potential at 50% loss conductivity [MPa]  

% Leaf
PsiL00_H     = [-1     -1      NaN      NaN    -1.2     -1.0    ]; % Water Potential at 2% loss conductivity [MPa] 
PsiL50_H     = [-3.2   -3.2    NaN      NaN    -4.0     -3.0    ]; % Water Potential at 50% loss conductivity [MPa]  
Kleaf_max_H  = [10     10      NaN      NaN    10       10      ]; % Leaf maximum hydraulic conductivity [mmolH20 m^2 leaf s /MPa]
Cl_H         = [1200   1200    NaN      NaN    1200     1200    ]; % Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]

% Xylem
Axyl_H       = [15     15      NaN      NaN    15       15      ]; % [cm^2 stem /m^2 PFT]
Kx_max_H     = [80000  80000   NaN      NaN    80000    80000   ]; % Xylem Conductivity specific for water. 5550-555550 [mmolH20 /m s MPa]  
PsiX50_H     = [-5     -5      NaN      NaN    -6       -4.5    ]; % Water Potential at 50% loss conductivity [MPa]
Cx_H         = [150    150     NaN      NaN    150      150     ]; % [kg / m^3 sapwood MPa]

% Stomata
Psi_sto_00_L = [NaN    NaN    -0.5      -1     NaN      NaN     ]; % Water Potential at 2% loss conductivity  [MPa]  
Psi_sto_50_L = [NaN    NaN    -2.8      -3     NaN      NaN     ]; % Water Potential at 50% loss conductivity [MPa]  

% Leaf
PsiL00_L     = [NaN    NaN    -1        -2.5   NaN      NaN     ]; % Water Potential at 2% loss conductivity [MPa]  
PsiL50_L     = [NaN    NaN    -3.5      -4.5   NaN      NaN     ]; % Water Potential at 50% loss conductivity [MPa]  
Kleaf_max_L  = [NaN    NaN    5         5      NaN      NaN     ]; % [mmolH20 m^2 leaf s /MPa]
Cl_L         = [NaN    NaN    1200      1200   NaN      NaN     ]; % Leaf capacitance [mmolH20 / m^2 leaf MPa] [500 - 3000]

% Xylem
Axyl_L       = [NaN    NaN    0         0      NaN      NaN     ]; % Xylem area over PFT area [cm^2 stem /m^2 PFT]
Kx_max_L     = [NaN    NaN    80000     80000  NaN      NaN     ]; % Xylem Conductivity specific for water;5550-555550 [mmolH20 /m s MPa]  
PsiX50_L     = [NaN    NaN    -4.5      -9     NaN      NaN     ]; % Water Potential at 50% loss conductivity[MPa] 
Cx_L         = [NaN    NaN    150       150    NaN      NaN     ]; % Steam capacitance low vegetation [kg / m^3 sapwood MPa]

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
%categories    [fir    larch   grass    shrub  BLever   BLdec  ] 
PsiG50_H   =   [-0.5   -0.8    NaN      NaN    NaN      NaN    ]; % Water potential at 50% impairment of growth and allocation control [MPa]
PsiG99_H   =   [-2.5   -2.5    NaN      NaN    NaN      NaN    ]; % Water potential at 90% impairment of growth and allocation control [MPa]
gcoef_H    =   [3.5    3.5     NaN      NaN    NaN      NaN    ]; % Parameter for maximum growth in perfect conditions, related to Env. controls of growth [gC/m2 day]

% Low vegetation
PsiG50_L   =   [NaN    NaN     -2.8     -3     NaN      NaN    ]; % Water potential at 50% impairment of growth and allocation control [MPa]
PsiG99_L   =   [NaN    NaN     -4       -4.5   NaN      NaN    ]; % Water potential at 90% impairment of growth and allocation control [MPa]
gcoef_L    =   [NaN    NaN     3.5      3.5    NaN      NaN    ]; % Parameter for maximum growth in perfect conditions, related to Env. controls of growth[gC/m2 day]

% Selection based on II
PsiG50_H=PsiG50_H(II); PsiG99_H=PsiG99_H(II); gcoef_H=gcoef_H(II);
PsiG50_L=PsiG50_L(II); PsiG99_L=PsiG99_L(II); gcoef_L=gcoef_L(II);

%% Vegetation Optical Parameter
%categories    [fir    larch   grass    shrub  BLever   BLdec  ] 
OPT_PROP_H  =  [2      3       0        0      5        7      ];   % =PFT_Class for "Veg_Optical_Parameter"-function
OPT_PROP_L  =  [0      0       13       2      0        0      ];

% Selection based on II
OPT_PROP_H = OPT_PROP_H(II);
OPT_PROP_L = OPT_PROP_L(II);

%debugger
%disp(num2str(OPT_PROP_L))

for i=1:cc_aux
    %%%%%%%% Vegetation Optical Parameter
    [PFT_opt_H(i)]=Veg_Optical_Parameter(OPT_PROP_H(i));
    [PFT_opt_L(i)]=Veg_Optical_Parameter(OPT_PROP_L(i));
end

%categories    [fir    larch   grass    shrub  BLever   BLdec ] 
OM_H  =        [1      1       NaN      NaN    1        1     ]; % Within canopy clumping factor [-]
OM_L  =        [NaN    NaN     1        1      NaN      NaN   ]; % Within canopy clumping factor [-]

% Selection based on II
OM_H=OM_H(II); OM_L=OM_L(II);

%% Specific leaf area of litter
Sllit = 2 ; % Litter Specific Leaf area [m2 Litter / kg DM]

%% High Vegetation

%categories   [fir   larch      grass  shrub  BLever       BLdec  ]  
aSE_H    =    [0        1       NaN    NaN    0            1      ]; % Allocation to reserve carbohydrate Values: 1 for Seasonal Plant and 0 for Evergreen
Sl_H     =    [0.010    0.025   NaN    NaN    0.016        0.020  ]; % Specific leaf area of  biomass [m^2 /gC]. Values: 0.05 -0.005
Nl_H     =    [42       26      NaN    NaN    40           28     ]; % Leaf Nitrogen Concentration [kgC/kgN ] 
r_H      =    [0.058    0.055   NaN    NaN    0.045        0.035  ]; % respiration rate at 10° [gC/gN d ]. Values: [0.066 -0.011]
gR_H     =    [0.25     0.25    NaN    NaN    0.25         0.25   ]; % Growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
Tcold_H  =    [-40      -3      NaN    NaN    3            3.5    ]; % Cold Leaf Shed [°C]
age_cr_H =    [950      180     NaN    NaN    365          180    ]; % Critical Leaf Age [day]
Trr_H    =    [0.25     3       NaN    NaN    0.5          3.5    ]; % Translocation rate [gC /m^2 d]
LtR_H    =    [0.8      0.8     NaN    NaN    1.0          0.9    ]; % Leaf to Root ratio maximum
eps_ac_H =    [0.25     1       NaN    NaN    0.5          1.0    ]; % Allocation to reserve parameter [0-1]
fab_H    =    [0.74     0.8     NaN    NaN    0.74         0.74   ]; % fraction above-ground sapwood and reserve
ff_r_H   =    [0.1      0.1     NaN    NaN    0.1          0.1    ]; % Reference allocation to Fruit and reproduction
Wm_H     =    [0        0       NaN    NaN    0            0      ]; % wood turnover coefficient [1/d]

dd_max_H = 1./[150      200     NaN    NaN   200           100    ]; % Death maximum for drought [1/d] 
dc_C_H   = 1./[5        10      NaN    NaN   365           10     ]; % Factor of increasing mortality for cold
drn_H    = 1./[900      1100    NaN    NaN   550           800    ]; % Turnover root  [1/d]
dsn_H    = 1./[1100     750     NaN    NaN   800           700    ]; % Normal transfer rate sapwood [1/d]
Mf_H     = 1./[80       50      NaN    NaN   80            80     ]; % Fruit maturation turnover [1/d]
Klf_H    = 1./[40       30      NaN    NaN   30            28     ]; % Dead Leaves fall turnover [1/d]

fbe_H = 1-fab_H; %% fraction below-ground sapwood and reserve
% [Stoich_H(1)]=Veg_Stoichiometric_Parameter(Nl_H(1));
% [Stoich_H(2)]=Veg_Stoichiometric_Parameter(Nl_H(2));
% [Stoich_H(3)]=Veg_Stoichiometric_Parameter(Nl_H(3));
% [Stoich_H(4)]=Veg_Stoichiometric_Parameter(Nl_H(4));
% [Stoich_H(5)]=Veg_Stoichiometric_Parameter(Nl_H(5));
% [Stoich_H(6)]=Veg_Stoichiometric_Parameter(Nl_H(6));

%% Phenology 
%categories   [fir     larch    grass  shrub  BLever    BLdec  ]  
Bfac_lo_H =   [0.99    0.99     NaN    NaN    0.95      0.95   ]; % Leaf Onset Water Stress
Bfac_ls_H =   [NaN     NaN      NaN    NaN    NaN       NaN    ]; % Not-used 
Tlo_H     =   [4.5     3.5      NaN    NaN    5.5       2.8    ]; % Mean Temperature for Leaf onset
Tls_H     =   [NaN     NaN      NaN    NaN    NaN       NaN    ]; % Not-used 
PAR_th_H  =   [NaN     NaN      NaN    NaN    NaN       NaN    ]; % Light Phenology Threshold 
dmg_H     =   [30      30       NaN    NaN    45        30     ]; % Day of Max Growth
LAI_min_H =   [0.001   0.01     NaN    NaN    0.001     0.01   ];
mjDay_H   =   [220     250      NaN    NaN    250       250    ]; % Maximum Julian day for leaf onset
LDay_min_H =  [12.8    12.7     NaN    NaN    12.1      11.7   ]; % Minimum Day duration for leaf onset
LDay_cr_H =   [11.8    11.6     NaN    NaN    11.8      12.0   ]; % Threshold for senescence day light [h]

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

%i = 3
for i=1:cc_aux
    [Stoich_H(i)]=Veg_Stoichiometric_Parameter(Nl_H(i));
    [ParEx_H(i)]=Exudation_Parameter(0);
    [Mpar_H(i)]=Vegetation_Management_Parameter;
end

%% Low Vegetation 
%categories   [fir     larch    grass  shrub    BLever    BLdec  ]  
aSE_L   =     [NaN     NaN      2      1        NaN       NaN    ]; % Allocation to reserve carbohydrate -- 1 Seasonal Plant --  0 Evergreen
Sl_L    =     [NaN     NaN      0.016  0.018    NaN       NaN    ]; % Specific leaf area of  biomass [m^2 /gC] 0.05 -0.005
Nl_L    =     [NaN     NaN      23     35       NaN       NaN    ]; % Leaf Nitrogen Concentration [kgC/kgN ]
r_L     =     [NaN     NaN      0.025  0.025    NaN       NaN    ]; % Respiration rate at 10° [gC/gN d ]  [0.066 -0.011]
gR_L    =     [NaN     NaN      0.25   0.25     NaN       NaN    ]; % Growth respiration  [] -- [Rg/(GPP-Rm)] [0.22 - 0.28]
Tcold_L =     [NaN     NaN      3      2        NaN       NaN    ]; % Cold Leaf Shed [°C] 


dd_max_L = 1./[NaN     NaN      45     100      NaN       NaN    ]; % Death maximum for drought [1/d]
dc_C_L   = 1./[NaN     NaN      32     15       NaN       NaN    ]; % Factor of increasing mortality for cold
drn_L    = 1./[NaN     NaN      950    1600     NaN       NaN    ]; % turnover root  [1/d]
dsn_L    = 1./[NaN     NaN      365    1800     NaN       NaN    ]; % normal transfer rate sapwood [1/d]
Mf_L     = 1./[NaN     NaN      50     50       NaN       NaN    ]; % Fruit maturation turnover [1/d]
Klf_L    = 1./[NaN     NaN      20     40       NaN       NaN    ]; % Dead Leaves fall turnover [1/d]

age_cr_L =    [NaN     NaN      180    180      NaN       NaN    ]; % [day] Critical Leaf Age
Trr_L    =    [NaN     NaN      2      0.6      NaN       NaN    ]; % Translocation rate [gC /m^2 d]
LtR_L    =    [NaN     NaN      0.7    1        NaN       NaN    ]; % Leaf to Root ratio maximum
Wm_L     =    [NaN     NaN      0      0        NaN       NaN    ] ;% wood turnover coefficient [1/d]
eps_ac_L =    [NaN     NaN      1      1        NaN       NaN    ]; % Allocation to reserve parameter [0-1]
fab_L    =    [NaN     NaN      0      0.75     NaN       NaN    ]; % fraction above-ground sapwood and reserve
fbe_L    = 1-fab_L;                                                           % fraction below-ground sapwood and reserve
ff_r_L   =    [NaN     NaN      0.1    0.1      NaN       NaN    ]; % Reference allocation to Fruit and reproduction

% Phenology 
Bfac_lo_L  =  [NaN     NaN      0.99   0.99     NaN       NaN    ]; % Leaf Onset Water Stress
Bfac_ls_L  =  [NaN     NaN      0.15   NaN      NaN       NaN    ]; % 
Tlo_L      =  [NaN     NaN      2.5    2.5      NaN       NaN    ]; % Mean Temperature for Leaf onset
Tls_L      =  [NaN     NaN      NaN    NaN      NaN       NaN    ]; % Not-used 
PAR_th_L   =  [NaN     NaN      NaN    NaN      NaN       NaN    ]; % Light Phenology Threshold 
dmg_L      =  [NaN     NaN      20     25       NaN       NaN    ]; % Day of Max Growth
LAI_min_L  =  [NaN     NaN      0.05   0.001    NaN       NaN    ];
mjDay_L    =  [NaN     NaN      250    180      NaN       NaN    ]; % Maximum Julian day for leaf onset
LDay_min_L =  [NaN     NaN      12     12.2     NaN       NaN    ]; % Minimum Day duration for leaf onset
LDay_cr_L  =  [NaN     NaN      12     12       NaN       NaN    ]; % Threshold for senescence day light [h]

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

for i=1:cc_aux
    [Stoich_L(i)]=Veg_Stoichiometric_Parameter(Nl_L(i));
    [ParEx_L(i)]=Exudation_Parameter(0);
    [Mpar_L(i)]=Vegetation_Management_Parameter;
end

%% INITIAL CONDITIONS
%==========================================================================
%
%==========================================================================

L_day=zeros(NNd,1);
for j_aux=2:24:NN
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j_aux)]= SetSunVariables(Datam(j_aux,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end
Lmax_day = max(L_day);
%clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')

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

%%%%%%%%%%%%%%%%%

%% CARBON POOLS
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
LAI_H(1,:) = LAI_Htm1(ij,:);  Rrootl_H(1,:)= Rrootl_Htm1(ij); 
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
Nreserve_H(1,:)= Nreserve_Htm1(ij); Preserve_H(1,:)= Preserve_Htm1(ij); Kreserve_H(1,:)= Kreserve_Htm1(ij);
FNC_H(1,:)=FNC_Htm1(ij); NupI_H(1,:,:)= NupI_Htm1(ij,:,:); Nreserve_L(1,:)= Nreserve_Ltm1(ij);
Preserve_L(1,:)=Preserve_Ltm1(ij); Kreserve_L(1,:)=Kreserve_Ltm1(ij); FNC_L(1,:)=FNC_Ltm1(ij); NupI_L(1,:,:)= NupI_Ltm1(ij,:,:);
%%%%%%%%%%%%%%%%%%
TdpI_H(1,:)=TdpI_Htm1(ij); TdpI_L(1,:)=TdpI_Ltm1(ij);






%%

%run(PARAM_IC)
Restating_parameters2; %uses cc_aux instead of cc

if length(Oice)==1
    Oice=zeros(NN,ms);
end

%%%%%%%%%%
Tdeb =zeros(NN,max(1,length(Zs_deb)-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Check_Land_Cover_Fractions(Crock,Curb,Cwat,Cbare,Ccrown);
CcrownFIX = Ccrown;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Lateral Contribution
q_runon=zeros(NN,1); %%[mm/h]
Qi_in=zeros(NN,ms); %%[mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%
tic ;
CK1=zeros(NN,1);  CK2=zeros(NN,1);
%profile on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=1; %%% time dt = 1 day  %%%%%%%%%%%%%%%%
Tstm0= Ts(1);
%%% i time dt = 1h
%%% ms soil layer
%%% cc Crown Area present

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% NUMERICAL METHODS OPTIONS (Numerical tolerance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Opt_CR = optimset('TolX',3);
%Opt_ST = optimset('TolX',0.1);
%Opt_ST = optimset('TolX',0.1,'Display','iter');
Opt_CR = optimset('TolFun',1);%,'UseParallel','always');   % for internal CO2 computation
Opt_ST = optimset('TolFun',0.1);%,'UseParallel','always'); % for surface temperature
Opt_ST2 = optimset('TolFun',0.1,'Display','off');
OPT_SM=  odeset('AbsTol',0.05,'MaxStep',dth);              % for soil moisture
OPT_VD=  odeset('AbsTol',0.05);                            % for carbon budget
OPT_PH= odeset('AbsTol',0.01);                             % for internal plant hydraulic
OPT_STh = odeset('AbsTol',0.02);                           % for heat transfer
OPT_VegSnow = 1;                                           % Option for computing energy budget of vegetation when there is snow at the ground   
OPT_SoilTemp = 1;                                          % Option for computing soil temperature or not. 
OPT_FR_SOIL = 1;      % whether you want soil freezing on or off. Mike found that it was causing problems for Maipo streamflow, so I turned it off
OPT_min_SPD = 0.006;  % [m] minimum snow pack depth to have a multilayer snow 

%%%%
OPT_VCA = 0;
OPT_ALLOME = 0;
OPT_DROOT = 0;
OPT_PlantHydr = 0;
OPT_EnvLimitGrowth = 0;
OPT_WET = 0;
Wlev=zeros(NN,1)
Wlevm1 = Rd(1)+Rh(1);
%{
if  not(exist('OPT_VCA','var'))
    OPT_VCA = 0;
    OPT_ALLOME = 0;
end
if  not(exist('OPT_DROOT','var'))
    OPT_DROOT = 0;
end
if  not(exist('OPT_PlantHydr','var'))
    OPT_PlantHydr = 0;
end
if  not(exist('OPT_EnvLimitGrowth','var'))
    OPT_EnvLimitGrowth = 0;
end
if  not(exist('OPT_WET','var'))
    OPT_WET = 0;
else
    if  not(exist('Wlev','var'))
        Wlev=zeros(NN,1);
    end
    Wlevm1 = Rd(1)+Rh(1);
end
%}

%% Debugger
%disp(strcat('MAIN_FRAME, Before iter:',num2str(size(Ca,2))))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Iterations for each time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:NN    
    %%%%Displaying iteration in the command window    
    if  (mod(i,1000) == 0) || (i == 2)
        disp('Iter:'); disp(i);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pdind = [max(1,i-24):i-1]; %% previous day indexes
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% Model calculations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (Datam(i,4)==1) % Condition on date
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j=j+1; [jDay(j)]= JULIAN_DAY(Datam(i,:));
        [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)] = SetSunVariables(Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
        % PARFOR does not likethis
        %clear h_S delta_S zeta_S T_sunrise T_sunset 
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Biogeochemistry submodule
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM]=Biogeo_environment(Tdp(pdind,:),O(pdind,:),V(pdind,:),Soil_Param,Phy,SPAR,Bio_Zs); % sum(V(i,:))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if OPT_SoilBiogeochemistry == 1
            IS= Ccrown*squeeze(ISOIL_L(j-1,:,:)) + Ccrown*squeeze(ISOIL_H(j-1,:,:));
            REXMY= Ccrown*squeeze(Rexmy_L(j-1,:,:)) + Ccrown*squeeze(Rexmy_H(j-1,:,:));
            FireA = 1*((sum(ManIH==-5) + sum(ManIL==-5)) > 0);
            
            %% BIOGEO_UNIT
            %==============================================================
            %
            %==============================================================
            [P(j,:),LEAK_NH4(j),LEAK_NO3(j),LEAK_P(j),LEAK_K(j),LEAK_DOC(j),LEAK_DON(j),LEAK_DOP(j),...
                R_NH4(j),R_NO3(j),R_P(j),R_K(j),R_DOC(j),R_DON(j),R_DOP(j),...
                Nuptake_H(j,:),Puptake_H(j,:),Kuptake_H(j,:),Nuptake_L(j,:),Puptake_L(j,:),Kuptake_L(j,:),RexmyI(j,:),...
                R_litter(j),R_microbe(j),R_litter_sur(j),R_ew(j),VOL(j),N2flx(j),Min_N(j),Min_P(j),...
                R_bacteria(j),RmycAM(j),RmycEM(j),Prod_B(j),Prod_F(j),BfixN(j),NavlI(j,:),LitFirEmi(j,:)]=BIOGEO_UNIT(P(j-1,:),IS,Zbio,sum(Bio_Zs.*rsd),PHs,Tdp_bio,mean(Ta(pdind)),Psi_bio,Se_bio,Se_fc,VSUM,VTSUM,...
                Ccrown,Bio_Zs,RfH_Zs,RfL_Zs,sum(Lk(pdind)),sum(Rd(pdind)),sum(Rh(pdind)),sum(Pr(pdind)),sum(T_H(pdind,:),1),sum(T_L(pdind,:),1),B_H(j-1,:,3),B_L(j-1,:,3),LAI_H(j-1,:),LAI_L(j-1,:),...
                SupN_H(j-1,:),SupP_H(j-1,:),SupK_H(j-1,:),SupN_L(j-1,:),SupP_L(j-1,:),SupK_L(j-1,:),...
                REXMY,RexmyI(j-1,:),ExEM,NavlI(j-1,:),Pcla,Psan,...
                B_IO,jDay(j),FireA,0);
            %%%%%%%%%%%%% End BIOGEO_UNIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%
            BLit(j,:)= 0.002*sum(P(j,1:5))*Ccrown; % Total litter content [kg DM / m2]
            Bam =  P(j,20); %%[gC/m2]
            Bem =  P(j,21); %%[gC/m2]
        else % if Biogeochemical module off, then paramters are set zero            
            BLit(j,:)= 0.0 ; % Total litter content  (0.100) [kg DM / m2]
            %%% Plant uptake
            Nuptake_H(j,:)= 0.0; % of nitrogen - high vegetation
            Puptake_H(j,:)= 0.0; % of phosphorus - high vegetation
            Kuptake_H(j,:)= 0.0; % of potassium - high vegetation [gK/m^2 day]
            %%% Pland uptake
            Nuptake_L(j,:)= 0.0;  % of nitrogen [gN/m^2 day]
            Puptake_L(j,:)= 0.0;  % of phosphorus [gN/m^2 day]
            Kuptake_L(j,:)= 0.0;  % of potassium [gN/m^2 day]
            %%%%
            NavlI(j,:)=[0 0 0];
            %P(j,:)=0.0;
            %LEAK_NH4(j)= 0.0; LEAK_NO3(j)= 0.0; LEAK_P(j)= 0.0; LEAK_K(j) = 0.0; LEAK_DOC(j)= 0.0;
            %R_litter(j)= 0.0; R_microbe(j)= 0.0; R_ew(j)=0; R_litter_sur(j)= 0.0; VOL(j)= 0.0; N2flx(j)= 0.0; BfixN(j)=0; LitFirEmi(j,:)=0;
            %LEAK_DOP(j)= 0.0; LEAK_DON(j)= 0.0; Min_N(j)=0; Min_P(j) =0; R_bacteria(j)=0; RmycAM(j)=0; RmycEM(j)=0;
            %RexmyI(j,:)=0.0;
            Bam=NaN; Bem=NaN;
        end %%% End Biogeochemical module 
        
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%% Calculations per Crown areas 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%% projected area n-coordinate
        for cc=1:length(Ccrown)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%% High vegetation for Crown areas
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if (ZR95_H(cc) > 0) || (ZRmax_H(cc) > 0)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%% VEGGIE UNIT
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [LAI_H(j,cc),B_H(j,cc,:),NPP_H(j,cc),ANPP_H(j,cc),Rg_H(j,cc),RA_H(j,cc),Rms_H(j,cc),Rmr_H(j,cc),Rmc_H(j,cc),PHE_S_H(j,cc),...
                    dflo_H(j,cc),AgeL_H(j,cc),e_rel_H(j,cc),e_relN_H(j,cc),LAIdead_H(j,cc),NBLeaf_H(j,cc),Sr_H(j,cc),Slf_H(j,cc),Sfr_H(j,cc),Swm_H(j,cc),Sll_H(j,cc),Rexmy_H(j,cc,:),Rrootl_H(j,cc),...
                    AgeDL_H(j,cc),Bfac_dayH(j,cc),Bfac_weekH(j,cc),NPPI_H(j,cc),TdpI_H(j,cc),NupI_H(j,cc,:),PARI_H(j,cc,:),NBLI_H(j,cc),RB_H(j,cc,:),FNC_H(j,cc),Nreserve_H(j,cc),Preserve_H(j,cc),Kreserve_H(j,cc),...
                    rNc_H(j,cc),rPc_H(j,cc),rKc_H(j,cc),ManIH(cc)]= VEGGIE_UNIT(B_H(j-1,cc,:),PHE_S_H(j-1,cc),dflo_H(j-1,cc),AgeL_H(j-1,cc),AgeDL_H(j-1,cc),...
                    Ta(pdind),Tdp_H(pdind,cc),PARB(pdind)+PARD(pdind),Psi_x_H(pdind,cc),Psi_l_H(pdind,cc),An_H(pdind,cc),Rdark_H(pdind,cc),NPP_H(j-1,cc),jDay(j),Datam(i,:),...
                    NPPI_H(j-1,cc),TdpI_H(j-1,cc),Bfac_weekH(j-1,cc),NupI_H(j-1,cc,:),NavlI(j,:),PARI_H(j-1,cc,:),NBLI_H(j-1,cc),NBLeaf_H(j-1,cc),...
                    L_day(j),Lmax_day,VegH_Param_Dyn,cc,...
                    Nreserve_H(j-1,cc),Preserve_H(j-1,cc),Kreserve_H(j-1,cc),Nuptake_H(j,cc),Puptake_H(j,cc),Kuptake_H(j,cc),FNC_H(j-1,cc),Se_bio,Tdp_bio,...
                    ParEx_H(cc),ExEM,Bam,Bem,Mpar_H(cc),TBio_Ht(j-1,cc),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);
                %%%%%%%%%%%%%%%%%% End VEGGIE_UNIT %%%%%%%%%%%%%%%%%%%%%%%%
                
                %% PLANT EXPORT
                %==========================================================
                % Function
                %==========================================================
                [TexC_H(j,cc),TexN_H(j,cc),TexP_H(j,cc),TexK_H(j,cc),TNIT_H(j,cc),TPHO_H(j,cc),TPOT_H(j,cc),NuLit_H(j,cc,:),Nreserve_H(j,cc),Preserve_H(j,cc),Kreserve_H(j,cc),...
                    SupN_H(j,cc),SupP_H(j,cc),SupK_H(j,cc),ISOIL_H(j,cc,:)]= Plant_Exports(B_H(j,cc,:),B_H(j-1,cc,:),NuLit_H(j-1,cc,:),...
                    Slf_H(j,cc),Sfr_H(j,cc),Swm_H(j,cc),Sll_H(j,cc),Sr_H(j,cc),Rexmy_H(j,cc,:),Stoich_H(cc),Mpar_H(cc),fab_H(cc),fbe_H(cc),RB_H(j,cc,:),...
                    Nreserve_H(j,cc),Preserve_H(j,cc),Kreserve_H(j,cc),rNc_H(j,cc),rPc_H(j,cc),rKc_H(j,cc),ManIH(cc),OPT_SoilBiogeochemistry);
                %%%%%%%%%%%%%%%% End Plant_Export %%%%%%%%%%%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%% Change in height SAI and Ccrown
                SAI_H(j,cc) =SAI_H(j-1,cc);
                if aSE_H(cc) == 2
                    [hc_H(j,cc)] = GrassHeight(LAI_H(j,cc),LAIdead_H(j,cc));
                elseif aSE_H(cc) == 5
                    [hc_H(j,cc),SAI_H(j,cc),B_H(j,:,:),Ccrown,Nreserve_H(j-1:j,:),Preserve_H(j-1:j,:),Kreserve_H(j-1:j,:),AgrHarNut(j,:)] = CropHeightType(LAI_H(j,cc),LAIdead_H(j,cc),cc,B_H(j,:,:),...
                        Ccrown,Nreserve_H(j,:),Preserve_H(j,:),Kreserve_H(j,:),ManIH,Mpar_H,VegH_Param_Dyn,OPT_SoilBiogeochemistry);
                    %%%%
                else
                    hc_H(j,cc)= hc_H(j-1,cc); %%%[m]
                    if OPT_VCA == 1
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%%%%%% Vegetation_Structural_Attributes
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        [Ccrown_t(j,cc),hc_H(j,cc),SAI_H(j,cc),BA_H(j,cc),Tden_H(j,cc),AgePl_H(j,cc),TBio_Ht(j,cc)]=Vegetation_Structural_Attributes(dtd,...
                            Ccrown_t(j-1,cc),B_H(j,cc,:),fab_H(cc),Tden_H(j-1,cc),AgePl_H(j-1,cc),OPT_ALLOME);
                        %%%% End of Vegetation_Structural_Attributes

                        Ccrown(cc) = CcrownFIX(cc)*Ccrown_t(j,cc);
                        B_H(j,cc,:)= B_H(j,cc,:)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Nreserve_H(j,cc)=Nreserve_H(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Preserve_H(j,cc)=Preserve_H(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Kreserve_H(j,cc)=Kreserve_H(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                    end
                end
                if OPT_DROOT == 1 
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%% Root_Depth_Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [ZR95_H(cc),RfH_Zs(cc,:)]=Root_Depth_Dynamics(CASE_ROOT,B_H(j,cc,:),B_H(j-1,cc,:),Rrootl_H(j,cc),Zs,ZR95_H(cc),ZR50_H(cc),ZRmax_H(cc),...
                        Bfac_dayH(j,cc),Psan,Tdp_H(pdind,cc),O(pdind,:),Soil_Param,a_root_H(cc));
                    %%%%%%%%%%%% End Root_Depth_Dynamics
                end
                ZR95_Ht(j,cc)=ZR95_H(cc); 
                
            else
                %%% variables are set zero
                LAI_H(j,cc)=0;B_H(j,cc,:)=0;NPP_H(j,cc)=0;ANPP_H(j,cc)=0;Rg_H(j,cc)=0;RA_H(j,cc)=0;Rms_H(j,cc)=0;
                Rmr_H(j,cc)=0;PHE_S_H(j,cc)=1;dflo_H(j,cc)=0;AgeL_H(j,cc)=0;e_rel_H(j,cc)=0; e_relN_H(j,cc)=0;
                SAI_H(j,cc)=0; hc_H(j,cc)=0;
                LAIdead_H(j,cc) =LAIdead_H(j-1,cc);
                Sr_H(j,cc)=0;Slf_H(j,cc)=0;Sfr_H(j,cc)=0;Swm_H(j,cc)=0; Rrootl_H(j,cc)=0;
                Rexmy_H(j,cc,:)=0;AgeDL_H(j,cc)=0;FNC_H(j,cc)=1; Nreserve_H(j,cc)=0;Preserve_H(j,cc)=0;Kreserve_H(j,cc)=0;
                rNc_H(j,cc)=1;rPc_H(j,cc)=1;rKc_H(j,cc)=1;
                Bfac_dayH(j,cc)=0; Bfac_weekH(j,cc)=0; NPPI_H(j,cc)=0; TdpI_H(j,cc)=0; NupI_H(j,cc,:)=0; RB_H(j,cc,:)=0;
                TexC_H(j,cc)=0;TexN_H(j,cc)=0;TexP_H(j,cc)=0;TexK_H(j,cc)=0;TNIT_H(j,cc)=0;TPHO_H(j,cc)=0;TPOT_H(j,cc)=0;
                ISOIL_H(j,cc,:)=0; NuLit_H(j,cc,:)=0;
            
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%% Low vegetation for Crown areas
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%% ZR95_L: Root depth 95 percentile - low vegetation
            %%% ZR_max_L: Maximum root depth - low vegetation
            if (ZR95_L(cc) > 0) || (ZRmax_L(cc) > 0)
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%% Veggie Unit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [LAI_L(j,cc),B_L(j,cc,:),NPP_L(j,cc),ANPP_L(j,cc),Rg_L(j,cc),RA_L(j,cc),Rms_L(j,cc),Rmr_L(j,cc),Rmc_L(j,cc),PHE_S_L(j,cc),...
                    dflo_L(j,cc),AgeL_L(j,cc),e_rel_L(j,cc),e_relN_L(j,cc),LAIdead_L(j,cc),NBLeaf_L(j,cc),Sr_L(j,cc),Slf_L(j,cc),Sfr_L(j,cc),Swm_L(j,cc),Sll_L(j,cc),Rexmy_L(j,cc,:),Rrootl_L(j,cc),...
                    AgeDL_L(j,cc),Bfac_dayL(j,cc),Bfac_weekL(j,cc),NPPI_L(j,cc),TdpI_L(j,cc),NupI_L(j,cc,:),PARI_L(j,cc,:),NBLI_L(j,cc),RB_L(j,cc,:),FNC_L(j,cc),Nreserve_L(j,cc),Preserve_L(j,cc),Kreserve_L(j,cc),...
                    rNc_L(j,cc),rPc_L(j,cc),rKc_L(j,cc),ManIL(cc)]= VEGGIE_UNIT(B_L(j-1,cc,:),PHE_S_L(j-1,cc),dflo_L(j-1,cc),AgeL_L(j-1,cc),AgeDL_L(j-1,cc),...
                    Ta(pdind),Tdp_L(pdind,cc),PARB(pdind)+PARD(pdind),Psi_x_L(pdind,cc),Psi_l_L(pdind,cc),An_L(pdind,cc),Rdark_L(pdind,cc),NPP_L(j-1,cc),jDay(j),Datam(i,:),...
                    NPPI_L(j-1,cc),TdpI_L(j-1,cc),Bfac_weekL(j-1,cc),NupI_L(j-1,cc,:),NavlI(j,:),PARI_L(j-1,cc,:),NBLI_L(j-1,cc),NBLeaf_L(j-1,cc),...
                    L_day(j),Lmax_day,VegL_Param_Dyn,cc,...
                    Nreserve_L(j-1,cc),Preserve_L(j-1,cc),Kreserve_L(j-1,cc),Nuptake_L(j,cc),Puptake_L(j,cc),Kuptake_L(j,cc),FNC_L(j-1,cc),Se_bio,Tdp_bio,...
                    ParEx_L(cc),ExEM,Bam,Bem,Mpar_L(cc),TBio_Lt(j-1,cc),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);
                %%%%%%%%%%%%%%%%%%% End of veggie unit
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%% Plant export
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [TexC_L(j,cc),TexN_L(j,cc),TexP_L(j,cc),TexK_L(j,cc),TNIT_L(j,cc),TPHO_L(j,cc),TPOT_L(j,cc),NuLit_L(j,cc,:),Nreserve_L(j,cc),Preserve_L(j,cc),Kreserve_L(j,cc),...
                    SupN_L(j,cc),SupP_L(j,cc),SupK_L(j,cc),ISOIL_L(j,cc,:)]= Plant_Exports(B_L(j,cc,:),B_L(j-1,cc,:),NuLit_L(j-1,cc,:),...
                    Slf_L(j,cc),Sfr_L(j,cc),Swm_L(j,cc),Sll_L(j,cc),Sr_L(j,cc),Rexmy_L(j,cc,:),Stoich_L(cc),Mpar_L(cc),fab_L(cc),fbe_L(cc),RB_L(j,cc,:),...
                    Nreserve_L(j,cc),Preserve_L(j,cc),Kreserve_L(j,cc),rNc_L(j,cc),rPc_L(j,cc),rKc_L(j,cc),ManIL(cc),OPT_SoilBiogeochemistry);
                %%%%%%%%%%% End of Plant exports

                %%%%%%%%%%%%%%%%%%%%%% Change in height SAI and Ccrown
                SAI_L(j,cc) =SAI_L(j-1,cc);
                if aSE_L(cc) == 2
                    [hc_L(j,cc)] = GrassHeight(LAI_L(j,cc),LAIdead_L(j,cc));
                elseif aSE_L(cc) == 5
                    [hc_L(j,cc),SAI_L(j,cc),B_L(j,:,:),Ccrown,Nreserve_L(j,:),Preserve_L(j,:),Kreserve_L(j,:),AgrHarNut(j,:)] = CropHeightType(LAI_L(j,cc),LAIdead_L(j,cc),cc,B_L(j,:,:),...
                        Ccrown,Nreserve_L(j,:),Preserve_L(j,:),Kreserve_L(j,:),ManIL,Mpar_L,VegL_Param_Dyn,OPT_SoilBiogeochemistry);
                else
                    hc_L(j,cc)= hc_L(j-1,cc); %%%[m]
                    if OPT_VCA == 1
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %%%%%% Vegetation_Structural_Attributes
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        [Ccrown_t(j,cc),hc_L(j,cc),SAI_L(j,cc),BA_L(j,cc),Tden_L(j,cc),AgePl_L(j,cc),TBio_Lt(j,cc)]=Vegetation_Structural_Attributes(dtd,...
                            Ccrown_t(j-1,cc),B_L(j,cc,:),fab_L(cc),Tden_L(j-1,cc),AgePl_L(j-1,cc),OPT_ALLOME);
                        %%%%%%%%%% End of Vegetation_Structural_Attributes

                        Ccrown(cc) = CcrownFIX(cc)*Ccrown_t(j,cc);
                        B_L(j,cc,:)= B_L(j,cc,:)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Nreserve_L(j,cc)=Nreserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Preserve_L(j,cc)=Preserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Kreserve_L(j,cc)=Kreserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                    end
                end

                if OPT_DROOT == 1 
                    %%%%%%%%% Root_Depth_Dynamics %%%%%%%%%%%%%%%%%%%%%%%%%
                    [ZR95_L(cc),RfL_Zs(cc,:)]=Root_Depth_Dynamics(CASE_ROOT,B_L(j,cc,:),B_L(j-1,cc,:),Rrootl_L(j,cc),Zs,ZR95_L(cc),ZR50_L(cc),ZRmax_L(cc),...
                        Bfac_dayL(j,cc),Psan,Tdp_L(pdind,cc),O(pdind,:),Soil_Param,a_root_L(cc));
                    %%%%%%%%%%%%%% End Root_Depth_Dynamics
                end
                ZR95_Lt(j,cc)=ZR95_L(cc); 

            else
                LAI_L(j,cc)=0;B_L(j,cc,:)=0;NPP_L(j,cc)=0;ANPP_L(j,cc)=0;Rg_L(j,cc)=0;RA_L(j,cc)=0;Rms_L(j,cc)=0;
                Rmr_L(j,cc)=0;PHE_S_L(j,cc)=1;dflo_L(j,cc)=0;AgeL_L(j,cc)=0;e_rel_L(j,cc)=0; e_relN_L(j,cc)=0;
                SAI_L(j,cc)=0; hc_L(j,cc)=0;
                LAIdead_L(j,cc) =LAIdead_L(j-1,cc);
                Sr_L(j,cc)=0;Slf_L(j,cc)=0;Sfr_L(j,cc)=0;Swm_L(j,cc)=0;  Rrootl_L(j,cc)=0;
                Rexmy_L(j,cc,:)=0;AgeDL_L(j,cc)=0;FNC_L(j,cc)=1;Nreserve_L(j,cc)=0;Preserve_L(j,cc)=0;Kreserve_L(j,cc)=0;
                rNc_L(j,cc)=1;rPc_L(j,cc)=1;rKc_L(j,cc)=1;
                Bfac_dayL(j,cc)=0; Bfac_weekL(j,cc)=0; NPPI_L(j,cc)=0; TdpI_L(j,cc)=0;  NupI_L(j,cc,:)=0; RB_L(j,cc,:)=0;
                TexC_L(j,cc)=0;TexN_L(j,cc)=0;TexP_L(j,cc)=0;TexK_L(j,cc)=0;TNIT_L(j,cc)=0;TPHO_L(j,cc)=0;TPOT_L(j,cc)=0;
                ISOIL_L(j,cc,:)=0; NuLit_L(j,cc,:)=0;
            end %%% End of conditions for root depth - low vegetation
            %%%%%%%%%%%%%%%%%
        end %%% End of calculations for Crowns
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%% Check on land cover fractions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if OPT_VCA == 0  %%% Solutions for Crops aSE==5 
            Ccrown_t(j,:) = Ccrown;
            Cbare = 1 - sum(Ccrown) - Cwat - Curb - Crock;
            Check_Land_Cover_Fractions(Crock,Curb,Cwat,Cbare,Ccrown);
        end
    end %End of iteration  

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Oil Palm option
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OPT_VCA == 1
        if((OPT_ALLOME == 1)&& (ZR95_L(2))>0) %%% Ad hoc solution for oil palm
            Ccrown_t(j,2) = 1-Ccrown_t(j,1);  % understory of Oil Palm
            Ccrown(2) = Ccrown_t(j,2);
            %%% Theoretically to preserve mass // but this is considered
            %%% lost understory in the  oil palm 
            %B_L(j,cc,:)= B_L(j,cc,:)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
            %Nreserve_L(j,cc)=Nreserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
            %Preserve_L(j,cc)=Preserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
            %Kreserve_L(j,cc)=Kreserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
        end
        Cbare = 1 - sum(Ccrown);
        if zatm < (max(max(hc_H),max(hc_L))+2)
            zatm = zatm + 2;
        end
    end
    
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Wetland option
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Qi_in(i,:)=Qi_out(i-1,:); % Incoming and outcomming lateral subsurface flow
    q_runon(i)=0; %Rd(i-1)+Rh(i-1); Runon
    if OPT_WET == 1 %%% Wetland option with water standing in the surface
        q_runon(i)=Wlev(i);
        if isnan(Wlev(i))
            q_runon(i)=Wlevm1;
        end
        NetWatWet(i) = q_runon(i) - Wlevm1 ;
    end
    
    % Debugging
    % disp(size(Ca))
    
   
    %% HYDROLOGIC UNIT
    %======================================================================
    %
    %======================================================================
    [V(i,:),Vice(i,:),O(i,:),Oice(i,:),ZWT(i),OF(i),OS(i),OH(i,:),OL(i,:),Psi_s_H(i,:),Psi_s_L(i,:),Rd(i),Qi_out(i,:),WTR(i,:),...
        Rh(i),Lk(i),f(i),WIS(i),Ts(i),Pr_sno(i),Pr_liq(i),Csno(i),Cice(i),NDVI(i),rb_H(i,:),rb_L(i,:),rs_sunH(i,:),...
        rs_sunL(i,:),rs_shdH(i,:),rs_shdL(i,:),r_litter(i,:),...
        An_L(i,:),An_H(i,:),Rdark_L(i,:),Rdark_H(i,:),Ci_sunH(i,:),Ci_sunL(i,:),Ci_shdH(i,:),Ci_shdL(i,:),...
        rap_H(i,:),rap_L(i,:),r_soil(i),b_soil(i),alp_soil(i),ra(i),Rn(i),...
        H(i),QE(i),Qv(i),Lpho(i),T_H(i,:),T_L(i,:),EIn_H(i,:),EIn_L(i,:),EG(i),ESN(i),ESN_In(i),ELitter(i),EWAT(i),EICE(i),EIn_urb(i),EIn_rock(i),dw_SNO(i),...
        G(i),Gfin(i),Tdp(i,:),Tdpsnow(i,:),Tdeb(i,:),Tdamp(i),Tice(i),Tdp_H(i,:),Tdp_L(i,:),SWE(i),SND(i),ros(i),In_SWE(i),SP_wc(i),WR_SP(i),U_SWE(i),NIn_SWE(i),dQ(i),Qfm(i),t_sls(i),DQ(i),DT(i),...
        WAT(i),ICE(i),ICE_D(i),IP_wc(i),WR_IP(i),NIce(i),Cicew(i),Csnow(i),FROCK(i),Imelt(i),Smelt(i),...
        In_H(i,:),In_L(i,:),In_Litter(i),In_urb(i),In_rock(i),Dr_H(i,:),Dr_L(i,:),SE_rock(i),SE_urb(i),Lk_wat(i),Lk_rock(i),er(i),...
        gsr_H(i,:),Psi_x_H(i,:),Psi_l_H(i,:),Jsx_H(i,:),Jxl_H(i,:),Kleaf_H(i,:),Kx_H(i,:),Vx_H(i,:),Vl_H(i,:),...
        gsr_L(i,:),Psi_x_L(i,:),Psi_l_L(i,:),Jsx_L(i,:),Jxl_L(i,:),Kleaf_L(i,:),Kx_L(i,:),Vx_L(i,:),Vl_L(i,:),...
        fapar_H(i,:),fapar_L(i,:),SIF_H(i,:),SIF_L(i,:),...
        snow_alb,tau_sno(i),e_sno(i),Ws_under(i),dQVEG(i),HV(i),QEV(i),TsV(i),Ts_under(i),EK(i),...
        POT(i,:)]=HYDROLOGIC_UNIT(V(i-1,:), ...
        Oice(i-1,:),aR,Zs,...
        EvL_Zs,Inf_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Dz,ms,Kbot,Pr(i),Ta(i),Ds(i),Ws(i),zatm,Ts(i-1),Ts_under(i-1),IrD(i),dt,dth,ea(i),N(i),Pre(i),Tstm0,...
        LAI_H(j,:),SAI_H(j,:),LAI_L(j,:),SAI_L(j,:),LAIdead_H(j,:),LAIdead_L(j,:),Rrootl_H(j,:),Rrootl_L(j,:),BLit(j,:),Sllit,Kct,...
        Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft,...
        Ccrown,Cbare,Crock,Curb,Cwat,...
        Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,...
        Zs_deb,Deb_Par,...
        ZR95_H,ZR95_L,...
        SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i),SvF,...
        SND(i-1),snow_alb,Color_Class,OM_H,OM_L,...
        PFT_opt_H,PFT_opt_L,hc_H(j,:),hc_L(j,:),d_leaf_H,d_leaf_L,...
        Ca(i),Oa,Ci_sunH(i-1,:),Ci_shdH(i-1,:),Ci_sunL(i-1,:),Ci_shdL(i-1,:),...
        e_rel_H(j,:),e_relN_H(j,:),e_rel_L(j,:),e_relN_L(j,:),...
        e_sno(i-1),In_H(i-1,:),In_L(i-1,:),In_Litter(i-1),In_urb(i-1),In_rock(i-1),SWE(i-1),In_SWE(i-1),....
        Tdeb(i-1,:),Tdp(i-1,:),Tdpsnow(i-1,:),Tdamp(i-1),Tice(i-1),...
        WAT(i-1),ICE(i-1),IP_wc(i-1),ICE_D(i-1),Cicew(i-1),...
        Vx_H(i-1,:),Vl_H(i-1,:),Vx_L(i-1,:),Vl_L(i-1,:),Psi_x_H(i-1,:),Psi_l_H(i-1,:),Psi_x_L(i-1,:),Psi_l_L(i-1,:),...
        FROCK(i-1),Krock,Ws_under(i-1),...
        Tdew(i),t_sls(i-1),ros(i-1),SP_wc(i-1),fpr,Pr_sno(pdind),...
        Urb_Par,In_max_urb,In_max_rock,K_usle,tau_sno(i-1),Ta(pdind),Slo_top,Slo_pot,Asur,Ared,aTop,EK(i-1),q_runon(i),Qi_in(i,:),...
        pow_dis,a_dis,Salt(i),...
        SPAR,SN,OPT_min_SPD,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST,Opt_ST2,OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH, ...
        parameterize_phase, hSTL, Albsno_method);    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% End of HYDROLOGIC UNIT %%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%% Prognostic temperature
    %Ts: Soil/Snow prognostic temperature for the energy balance
    Tstm0 =2*Ts(i)-Ts(i-1);

    %%%%%%%%% Radiation and albedo
    STOT= SAB1(i)+SAB2(i)+SAD1(i)+SAD2(i);
    ALB(i)= SAB1(i)/STOT*snow_alb.dir_vis + SAD1(i)/STOT*snow_alb.dif_vis + ... % Albedo
        SAB2(i)/STOT*snow_alb.dir_nir + SAD2(i)/STOT*snow_alb.dif_nir;
    
    %%%%%%%%% Numerical optimization
    if OPT_WET == 1 
        Wlevm1 = Rd(i)+Rh(i); 
    end 
    
    %%%%%%%%% Mass balance
    % v-coordinate
    % CK1: Check on mass balance
    CK1(i) = f(i)*dth*Asur*Ared + sum(V(i-1,:) - V(i,:))*Asur*Ared + sum(Vice(i-1,:) - Vice(i,:))*Asur*Ared - EG(i)*dth - Lk(i)*dth ...
        - sum(Qi_out(i,:))*dth -Rd(i) -sum(Jsx_L(i,:)).*dth -sum(Jsx_H(i,:)).*dth  + sum(Qi_in(i,:))*dth  ;
    
end  %End of iteration 



%% CHECKS
%==========================================================================
% CHECKS
%==========================================================================
%close(bau)
Computational_Time =toc;
%profile off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('COMPUTATIONAL TIME [h] ')
disp(Computational_Time/3600)
disp(' COMPUTATIONAL TIME [ms/cycle] ')
disp(1000*Computational_Time/NN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROF1 = profile('info');
%profile('status')
%profile viewer
if OPT_SoilBiogeochemistry == 1
    NEE = -(NPP_H+NPP_L)*Ccrown' + R_litter + R_microbe + R_ew; %% [gC/m2 day]
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% MASS BALANCE CHECK
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dV= (V(1,:)-V(end,:))*Asur*Ared +(Vice(1,:)-Vice(end,:))*Asur*Ared;
Ck = sum((Pr_liq(2:end)+Pr_sno(2:end))*dth)+ sum(IrD*dth) + sum(dV)...
    -sum(EG*dth)  -sum(ELitter*dth) -sum(sum(T_L)*dth)  -sum(sum(T_H)*dth) -sum(sum(EIn_L)*dth)  -sum(sum(EIn_H)*dth) -sum(Rh)...
    -sum(sum(Qi_out*dth)) -sum(Rd) + sum(sum(Qi_in*dth)) + sum(q_runon*dth)  - sum(EWAT*dth) ...
    - sum(ESN*dth) -sum(ESN_In*dth) -sum(EIn_urb*dth)- sum(EIn_rock*dth) + ...
    (SWE(1)-SWE(end)) + (In_SWE(1) - In_SWE(end)) + (SP_wc(1) -SP_wc(end))  + ...
    sum(In_H(1,:) - In_H(end,:)) + sum((In_L(1,:) - In_L(end,:)))+ (In_Litter(1) - In_Litter(end)) + (In_urb(1) - In_urb(end)) + (In_rock(1) - In_rock(end)) + ...
    - sum(EICE*dth) + (IP_wc(1) -IP_wc(end)) + (ICE(1) -ICE(end)) + (WAT(1) -WAT(end)) + (FROCK(1) -FROCK(end)) + ...
    Asur*( (Vx_H(1,:)-Vx_H(end,:)) + (Vl_H(1,:)-Vl_H(end,:))  + (Vx_L(1,:)-Vx_L(end,:)) + (Vl_L(1,:)-Vl_L(end,:)) )*Ccrown' ; %%%[mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%disp(Ck);
%disp(mean(DQ));
% T_H, T_L  EG, EIn_urb, EIn_rock, [mm/h]

%% Evaporation check
%==========================================================================
% QE: Latent heat [W m-2]
% Laten: Latent heat of vaporization in [J kg-1]
% ETen: [mm h-1]
% EG: Evaporation from Bare soil [mm h-1]
% ELitter: Evaporation from the Litter [mm h-1]
% ESN: Evaporation from the snowpack at the ground [mm h-1]
% ESN_In: Evaporation from intercepted snow [mm h-1]
% ET: Evapotranspiration [mm h-1]
% EWAT: Evaporation from water and ponds [mm h-1]
% EICE: Evaporation/sublimation from Ice [mm h-1]
% EIn_urb: Evaporation from Urban [mm h-1] 
% EIn_rock: Evaporation from Rocks [mm h-1] 
% EIn_H: Evaporation from intercepted water High Vegetation [mm h-1]
% EIn_L: Evaporation from intercepted water Low Vegetation [mm h-1]
% T_L: Transpiration Low Vegetation [mm h-1]
% T_H: Transpiration High Vegetation [mm h-1]
%==========================================================================

% Evaporation from latent heat using temperature
Laten= 1000*(2501.3 - 2.361*(Ta));
ETen = (QE)*1000*3600./(reshape(Laten,size(QE))*1000);

% Evaporation from model components
ET =  sum(T_H+EIn_H,2) + sum(T_L+EIn_L,2) +  EG +  ELitter + ESN + ESN_In + EWAT +  EICE+ EIn_urb + EIn_rock ;

%% Carbon balance check
%==========================================================================
% CARBON BALANCE CHECK
%==========================================================================
%{

for kj=1:cc
    if RA_H(end,kj) == 0
        dB_H= squeeze((B_H(1,kj,:)-B_H(end-1,kj,:)));
    else
        dB_H= squeeze((B_H(1,kj,:)-B_H(end,kj,:)));
    end
    if RA_L(end,kj) == 0
        dB_L= squeeze((B_L(1,kj,:)-B_L(end-1,kj,:)));
    else
        dB_L= squeeze((B_L(1,kj,:)-B_L(end,kj,:)));
    end
    %%%%%%%%
    %CkC_H(kj) =  sum(dB_H)+ sum(NPP_H(:,j)) - sum(Swm_H(:,j))-sum(Sfr_H(:,j))-sum(Sr_H(:,j))-sum(Slf_H(:,j)) - sum(Rexmy_H(:,j,:));%
    CkC_H(kj) =  sum(dB_H)+ sum(1.0368*An_H(:,kj)/24)- sum(Rmr_H(:,kj)) -sum(Rmc_H(:,kj)) -sum(Rms_H(:,kj)) -sum(Rg_H(:,kj))...
        -sum(TexC_H(:,kj));
    %- sum(Swm_H(:,kj))-sum(Sfr_H(:,kj))-sum(Sr_H(:,kj))-sum(Slf_H(:,kj)) - sum(Rexmy_H(:,kj,:));%
    % CkC_H(kj) =  sum(dB_H)+ sum(1.0368*An_H(:,kj)/24)- sum(Rmr_H(:,kj)) -sum(Rmc_H(:,kj)) -sum(Rms_H(:,kj)) -sum(Rg_H(:,kj))...
    %     -sum(RB_H(:,kj,:),[1,3])- sum(Swm_H(:,kj))-sum(Sfr_H(:,kj))-sum(Sr_H(:,kj))-sum(Slf_H(:,kj))-sum(Rexmy_H(:,kj,:),[1,3]);
    CkC_L(kj) =  sum(dB_L)+ sum(1.0368*An_L(:,kj)/24)- sum(Rmr_L(:,kj)) -sum(Rmc_L(:,kj)) -sum(Rms_L(:,kj)) -sum(Rg_L(:,kj))...
        -sum(TexC_L(:,kj));
    %- sum(Swm_L(:,kj))-sum(Sfr_L(:,kj))-sum(Sr_L(:,kj))-sum(Slf_L(:,kj)) - sum(Rexmy_L(:,kj,:));%
    %    CkC_L(kj) =  sum(dB_L)+ sum(1.0368*An_L(:,kj)/24)- sum(Rmr_L(:,kj)) -sum(Rmc_L(:,kj)) -sum(Rms_L(:,kj)) -sum(Rg_L(:,kj))...
    %    -sum(RB_L(:,kj,:),[1,3])- sum(Swm_L(:,kj))-sum(Sfr_L(:,kj))-sum(Sr_L(:,kj))-sum(Slf_L(:,kj))-sum(Rexmy_L(:,kj,:),[1,3]);
end
CkC_ALL = sum(CkC_H)+sum(CkC_L);


for kj=1:cc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% NUTRIENT BALANCE CHECK    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if RA_H(end,kj) == 0
        ed=length(Nreserve_H)-1;
    else
        ed=length(Nreserve_H);
    end

    dNres_H=Nreserve_H(1,kj) - Nreserve_H(ed,kj);
    CkN_H(kj)= dNres_H -  sum(TexN_H(:,kj))  + (TNIT_H(2,kj) - TNIT_H(ed,kj)) + sum(Nuptake_H(:,kj));
    dPres_H=Preserve_H(1,kj) - Preserve_H(ed,kj);
    CkP_H(kj)= dPres_H -  sum(TexP_H(:,kj))  + (TPHO_H(2,kj) - TPHO_H(ed,kj)) + sum(Puptake_H(:,kj));
    dKres_H=Kreserve_H(1,kj) - Kreserve_H(ed,kj);
    CkK_H(kj)= dKres_H -  sum(TexK_H(:,kj))  + (TPOT_H(2,kj) - TPOT_H(ed,kj)) + sum(Kuptake_H(:,kj));
    %%%%
    dNres_L=Nreserve_L(1,kj) - Nreserve_L(ed,kj);
    CkN_L(kj)= dNres_L -  sum(TexN_L(:,kj))  + (TNIT_L(2,kj) - TNIT_L(ed,kj)) + sum(Nuptake_L(:,kj));
    dPres_L=Preserve_L(1,kj) - Preserve_L(ed,kj);
    CkP_L(kj)= dPres_L -  sum(TexP_L(:,kj))  + (TPHO_L(2,kj) - TPHO_L(ed,kj)) + sum(Puptake_L(:,kj));
    dKres_L=Kreserve_L(1,kj) - Kreserve_L(ed,kj);
    CkK_L(kj)= dKres_L -  sum(TexK_L(:,kj))  + (TPOT_L(2,kj) - TPOT_L(ed,kj)) + sum(Kuptake_L(:,kj));
end
%}

%%%% Difference between carbon and nutrient ISOIL and Tex (0 except for
%%%% management cases)
%%% To add LitFirEmi to the litter budget
IS = zeros(18, 1);
Tex = 0;
for jj = 1:length(Ccrown_t)
    IS = IS + squeeze(ISOIL_L(jj, :, :))'*Ccrown_t(jj, :)' + squeeze(ISOIL_H(jj, :, :))'*Ccrown_t(jj, :)';
    Tex = Tex + TexC_L(jj,:)*Ccrown_t(jj, :)' + TexC_H(jj,:)*Ccrown_t(jj, :)';
end

CkExC = Tex - sum(IS(1:9));

%PARFOR does not like this
%clear dB_H  dB_L dNres_H  dNres_L  dPres_H dKres_H  dPres_L dKres_L ed
%clear Se_bio Psi_bio Tdp_bio VSUM VTSUM  IS Tex 
%clear  Tstm0 snow_alb  kj m p  pdind STOT  PARAM_IC
%clear Bam Bem FireA Se_fc REXMY









%% Post-compute calculations
%==========================================================================
% Negative evapotranspiration (ET) is dew
%==========================================================================
% post-compute sublimation from ESN
SSN = ESN.*(Ts<0);
%ET(ET < 0) = 0; 

%% Output manager
%==========================================================================
%{
As outputs:
    Hourly outputs
    Daily outputs
    Parameters
%}
%==========================================================================

% Make date usable for R
Date_R = char(Date); 

%% Parameter output
Param_t = table(Lat,Lon,Zbas,dbThick,'VariableNames',{'Lat','Lon','Zbas','dbThick'});
Param_t = [Param_t, struct2table(SnowIce_Param), struct2table(Deb_Par)];
Param_t = rows2vars(Param_t);
Param_t = renamevars(Param_t,{'OriginalVariableNames','Var1'},{'Parameter','Value'});

% Exporting as .txt
writetable(Param_t, strcat(outlocation,id_location,'_param.txt') )

% Here I manually choose the T&C outputs I want to save at each point.
Outputs_t = table(Date,EICE,ESN,SND,SWE,...
Ta,Ws,U,N,SAD1+SAD2+SAB1+SAB2,Pre,Pr,Pr_sno,ALB,Smelt,Imelt,SSN,ICE,ET, ETen ,QE,ros,'VariableNames',{ ...
'Date','EICE','ESN','SND','SWE',...
'Ta','Ws','U','N','Rsw',...
'Pre','Pr','Pr_sno','Albedo','Smelt','Imelt','SSN','ICE','ET','ETen','QE','ros'});

%% Hourly output
writetable(Outputs_t, strcat(outlocation,id_location,'_hourly_results.txt'))

%% Daily outputs
% If daily outputs are activated

if output_daily == 1

Outputs_tt = table2timetable(Outputs_t); 

Outputs_ds = retime(Outputs_tt,'daily',@nansum);
Outputs_dm = retime(Outputs_tt,'daily',@nanmean);

Outputs_d = Outputs_dm; 

Outputs_d.Pr = Outputs_ds.Pr;
Outputs_d.Pr_sno = Outputs_ds.Pr_sno;
Outputs_d.ET = Outputs_ds.ET;
Outputs_d.ETen = Outputs_ds.ETen;

Outputs_t = timetable2table(Outputs_d);

% Exporting as .txt
writetable(Outputs_t, strcat(outlocation,id_location,'_daily_results.txt'))
end 


end 

%% OTHER CALCULATIONS OUT OF THE MODEL
%==========================================================================
%
%==========================================================================
%{
% Resample the DEM to a new resolution (e.g., 30 meters)
DEM_resampled_out = resample(DEM, 500); 

% Display the original and resampled DEMs
imagesc(DEM);
figure;
imagesc(DEM_resampled);
%}