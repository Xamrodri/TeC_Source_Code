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
sel_point = 3; % selected position of the point for modelling in modelling_pixel
sel_time = 2;  % selected time for the variable x1s and x2s
sel_forc = 1;  % selection of forcing

%% Study site details
SITE = ["Cinca_Mid", "Tiber"]; 
FORCING = "ERA5Land";
UTM_zone = 33; % for Italy
DeltaGMT= 1; % for Italy
modelling_pixel = ["Cinca_Mid_OUT", "Apennine_out", "Monte_Terminillo", "Lago_di_Corbara", "Point_1"];
forcing_point = ["Monte_Terminillo"]; %Change in the future to modelling_pixel variable
%outlet_name = char(outlet_names);

%% Modelling period
%For some reason  "01-Jan-2008 00:00:00" does not work. Only  "01-Jan-2008 01:00:00"
x1s =  ["01-Nov-2022 00:00:00", "01-Jan-2008 01:00:00"]; % Starting point of the simulation
x2s =  ["01-Jun-2023 23:00:00", "30-Dec-2008 23:00:00"]; % Last timestep of the simulation

date_start = datetime(x1s(sel_time));
date_end = datetime(x2s(sel_time));

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
out = strcat(outlocation,'INIT_COND_', SITE(sel_basin) ,'_MultiPoint.mat'); % file path initial conditions

%% Dependencies
addpath(genpath([folder_path,'1_Functions'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path,'5_Common_inputs'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path,'3_Pyrenees_PointScale/2_Forcing'])); % Where is located the meteorological forcing and Shading matrix 
addpath(genpath([folder_path,'3_Pyrenees_PointScale/3_Inputs'])); % Add path to Ca_Data

%% Load DEM and geographical information
dtm_file = ["dtm_Cinca_Mid_250m.mat" "dtm_Tiber_250m.mat" "dtm_Tiber_250m.mat" "dtm_Tiber_250m.mat"]; 
res = 250; % simulation resolution [m]
disp(strcat('Model resolution: ',num2str(res)))

dtm_file_op = strcat(folder_path,'5_Common_inputs/',SITE(sel_basin),'/',dtm_file(sel_basin));
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
fn_alb_elev = strcat(SITE(sel_basin), '_Albedo_vs_elev.mat');

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

%% POIs
POI = readtable(strcat(folder_path,'3_PointScale_version/3_Inputs/2_Apennine/Apennine_MultiPoints.txt')); %import table with points info
[POI.LAT, POI.LON] = utm2ll(POI.UTM_X, POI.UTM_Y, UTM_zone);

LOC = find(strcmp(POI.Name, modelling_pixel(sel_point))); %Location within POI table

%% FOR LOOP for locations
Loc=30;
for loc = LOC  

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
id_location = char(string(POI.Name(loc))); %id
y_coord = POI.UTM_Y(loc);
x_coord = POI.UTM_X(loc);

pixelX = floor((x_coord - xllcorner) / cellsize) + 1;
pixelY = floor((y_coord - yllcorner) / cellsize) + 1;

%ij = POI.idx(loc);
ij = sub2ind(size(DTM),pixelY, pixelX); % Location
[j, i] = ind2sub(size(DTM), ij); % Location

Zbas = DTM(j,i); % Altitude
Lat = POI.LAT(loc);
Lon = POI.LON(loc);

%Debugging (for checking)
%line=reshape(DTM,num_cell,1);
%line(ij)


%% FORCING
%==========================================================================
% ERA5
%==========================================================================
forc_file = strcat(forc_path,'Forcing_ERA5_Land_',forcing_point(sel_forc),'_2008_corr_all.mat'); % Put here the path of where you downloaded the repository
load(forc_file); % Load forcing table for the current POI

Date_all=forc.Date; 

%define period and time zone info
x1=find(date_start == Date_all,1);
x2=find(date_end == Date_all,1);


%% Displaying modelling parameters
disp(strcat("Site selected: ", SITE(sel_basin)))
disp(['Forcing selected: ' char(FORCING)])
disp(['Running T&C for pixel: ' id_location])
disp(['Simulation period: ' datestr(date_start) ' to ' datestr(date_end)])

%% Fetch time and do date handling
Date = Date_all(x1:x2);
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

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
clear d1 d2

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
%}
%==========================================================================

% Period of forcing data
forcing = forc(x1:x2,:);
NN= height(forcing);%%% time Step

% Height of virtual station
zatm_hourly = repmat(2.00,height(forcing),1);

%categories    [fir     larch    grass  shrub  BLever    BLdec   ]  
zatm_surface = [18      18       2      2      18        18      ]; %Depend on vegetation
zatm_hourly_on = 0;

%% Precipitation
% Precipitation from ERA5Land
Pr=forcing.Total_Precipitation_HH;
Pr(isnan(Pr))=0;
Pr(Pr<0.01)=0;

%Precipitation from stations


%% Air pressure
% Air Pressure in mbar based on the PDF for variables and parameters of TC
% Pressure comes in Pa from ERA5
Pre=forcing.Pressure/100;    

%% Temperature
% 2m air temperature
Ta=forcing.Temperature;

%% Wind Speed
Ws=forcing.Wind_Speed; Ws(Ws < 0.01) = 0.01;

%% Relative humidity
% Divided by 100 to set the number in the range 0-1
U=forcing.RH/100;

%% Longwave radiation
% N can be cloud cover [-] or longwave incoming radiation [W m-2]
% Latm=forcing.LW_rad_downward_HH; % Latm:Incoming long wave radiation [W m-2]
% N=ones(NN,1); % cloud cover [-]
% N = forcing.LW_rad_downward_HH;
N = forcing.N;

%% Radiation partition
SAD1=forcing.SAD1; SAD2=forcing.SAD2; 
SAB1=forcing.SAB1; SAB2=forcing.SAB2;
PARB=forcing.PARB; PARD=forcing.PARD;

%% Albedo parameters
%Ameas = ones(NN,1);
alpha = 0; % switch for albedo
%Ameas_t=0; % albedo
%Aice_meas_on_hourly = ones(height(forcing),1)/2; % albedo
%Asno_meas_on_hourly = ones(height(forcing),1)/2; % albedo

%% Vapor pressure - Dew Point temperature
%esat/ea/Ds/Tdew
esat=forc.es;   % Vapour pressure at saturation (Pa)
ea=forc.ea;     % Vapour pressure (Pa)
Ds= esat - ea;  % Vapor Pressure Deficit (Pa)
Ds(Ds<0)=0; 
Tdew= forc.Dew_Point_Temp;

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
num_cell=numel(DTM);

%MASK = MASK.*0+1;
MASKn=reshape(MASK,num_cell,1);

if topo == 1
    %load topography data and narrow down to period 
    m = matfile(strcat(folder_path,'3_PointScale_version/2_Forcing/',SITE(sel_basin),'_ShF_',char(FORCING),'_2018.mat')); % ShF matrix created during pre-processing step

    x1_top=find(date_start==m.Top_Date,1);
    x2_top=find(date_end==m.Top_Date,1); 

    ShF = double(squeeze(m.ShF_S(ij,:)));
    ShF = ShF(x1_top:x2_top)';
    
    rho_g = 0.35; %%% Spatial Albedo
    zeta_S = m.zeta_Sts(x1_top:x2_top,1);
    h_S = m.h_Sts(x1_top:x2_top,1);
    clear e1 e2 Top_Date zeta_Sts H_Sts ShF_S
    SvF = m.SvF_S(ij,1); % Sky view factor at pixel ij
    Ct = m.Ct_S(i,j);
    Slo_top = m.Slo_top_S(ij,1); % Slope at pixel ij
    Aspect = m.Aspect_S(ij,1); % Aspect at pixel ij
    clear SvF_S Ct_S Aspect_S

    cos_fst = cos(atan(Slo_top)).*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;

    clear zeta_S 
  
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
switch ksv(ij)
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
        
        Cwat = 0.0; Curb = 0.0; Crock = 0.0; Cbare = 0.0;
        Ccrown = [0.1 0.1 0.8];
        
        %categories   [fir     larch    grass  shrub    BLever    BLdec]
        II =          [0       0        0      1        1         1    ]>0; 

        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types    
 
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

        Cwat = 0.1; Curb = 0.1; Crock = 0.0; Cbare = 0.0;
        Ccrown = [0.7 0.1];  
        cc=length(Ccrown);%% Crown area

        %categories   [fir     larch    grass  shrub  BLever    BLdec]
        II =          [0       0        1      1      0         0    ]>0;  
    
        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types 

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
        Cwat = 0; Curb = 0.0 ; Crock = 0.0; Cbare = 0.0;
        Ccrown = [0.5 0.5];
        
        %categories   [fir     larch    grass  shrub  BLever    BLdec ]
        II =          [0       0        1      1      0         0     ]>0;  
    
        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types 

    case 4 % Evergreen needleaves %   
        %    Case 4 includes from CORINE:
        %      1) Coniferous forest (1%)
        Cwat = 0; Curb = 0.0 ; Crock = 0.0; Cbare = 0.0;
        Ccrown = [1.0];
        
        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]
        II =          [0       0        0      0      1         0      ]>0;  
    
        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types 

    case 5 % Mediterranean shrublands %
        %    Case 5 includes from CORINE:
        %       1)  Transitional woodland-shrub (5%)
        %       2)  Natural grasslands (4.3%)
        %       3)  Sclerophyllous vegetation (0.4%)
        %       4)  Moors and heathland (>0.05%)
        
        Cwat = 0; Curb = 0.0 ; Crock = 0.1; Cbare = 0.1;        
        Ccrown = [0.6 0.1 0.1];        

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]
        II =          [0       0        0      1      1         1      ]>0;  

    
        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types 

    case 6 % Olives %
        %    Case 6 includes from CORINE:
        %       1) Olive groves (4%)

        Cwat = 0; Curb = 0.0 ; Crock = 0.0; Cbare = 0.0;
        Ccrown = [1.0];

        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]
        II =          [0       0        0      1      0         0      ]>0;  
    
        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types 

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

        Cwat = 0; Curb = 0.8 ; Crock = 0.0; Cbare = 0.1;
        Ccrown = [0.1];
                
        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]
        II =          [0       0        1      0      0         0      ]>0;  
        
        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types 

    case 8 % Rock %
        %    Case 8 includes from CORINE:
        %        1) Bare rocks (0.2%)
        %        2) Glaciers and perpetual snow (0%)

        Cwat = 0.0; Curb = 0.0 ; Crock = 1.0; Cbare = 0.0;
        Ccrown = [0.0];
        
        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]  
        II =          [0       0        0      0      0         0      ]>0; 

        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types

    case 9 % Water %
        %    Case 9 includes from CORINE:
        %        1) Water bodies (0.3%)
        %        2) Water courses (0.2%)

        %     It also includes, but not in Tiber basin:
        %       1) Coastal lagoons
        %       2) Estuaries
        %       3) Sea and Ocean


        Cwat = 1.0; Curb = 0.0 ; Crock = 0.0; Cbare = 0.0;
        Ccrown = [0.0];
        
        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]  
        II =          [0       0        0      0      0         0      ]>0;
  
        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types

    case 10 % Bare soils %
        %    Case 10 includes from CORINE:
        %        1) Mineral extraction sites (0.2%)
        %        2) Burnt areas (0.1%)
        %        3) Beaches - dunes - sands (>0.05%)
        %        4) Sparsely vegetated areas (0.7%)
        
        Cwat = 0.0; Curb = 0.0 ; Crock = 0.0; Cbare = 0.9;
        Ccrown = [0.1];
        
        %categories   [fir     larch    grass  shrub  BLever    BLdec  ]  
        II =          [0       0        0      1      0         0      ]>0;
  
        cc = length(Ccrown); % Crown area
        cc_max = length(Ccrown); % one vegetation types

    otherwise
        disp('INDEX FOR SOIL VEGETATION PARAMETER INCONSISTENT')
        return
end

zatm = max(zatm_surface(II)); %choose correct atmospheric reference height
%disp('here')
%% SOIL 
%==========================================================================
% Vector is divided by 100 to put the numbers within 0-1
% Original PSAN, PCLA and PORG come with values between 0-100, g/100g
%==========================================================================
PSAN=reshape(PSAN,num_cell,1)/100; Psan = PSAN(ij); % Soil sand content at pixel ij
PCLA=reshape(PCLA,num_cell,1)/100; Pcla = PCLA(ij); % Soil clay content at pixel ij
PORG=reshape(PORG,num_cell,1)/100; Porg= PORG(ij); % Soil organic content at pixel ij

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


if topo == 1
INIT_COND_v2(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,m.Slo_top_S,ksv,Ca,SNOWD,SNOWALB,out, II, ij);
load(out);
else 
INIT_COND_v2(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Slo_top_S,ksv,Ca,SNOWD,SNOWALB,out, II, ij);
load(out);
end

%end

%% RUN MODEL
%==========================================================================
% PARAM_IC: Define parameter file
% MAIN_FRAME: Contains the model
%==========================================================================
disp(II)
PARAM_IC = strcat(folder_path,'3_PointScale_version/3_Inputs/MOD_PARAM_Multipoint.m');
MAIN_FRAME; % Launch the main frame of T&C. Most of the things happen in this line of code

%% Post-compute calculations
%==========================================================================
% Negative evapotranspiration (ET) is dew
%==========================================================================
% post-compute sublimation from ESN
SSN = ESN.*(Ts<0);
ET(ET < 0) = 0; 

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