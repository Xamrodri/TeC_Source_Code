%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
% Created on Nov 25, 2024
% Author: MAXIMILIANO RODRIGUEZ
% Code originally from: ACHILLE JOUBERTON
% Area of Study: Apennines
% Region: Tiber basin
% Code explanation: This code launches TC model.
% Output variables are managed by the file OUTPUT_MANAGER_DIST_LABEL.
% Some code and names still referes to previous versions of the code.
%==========================================================================

%% CLEAR ALL
clc; clear;

% clear all   % Don't clear all such that it can be run on the HPC cluster without any issues
% delete(gcp('nocreate'))

%% DIRECTORIES
study_name = '2_Pyrenees_distributed';

%% INITIAL CONDITIONS
%==========================================================================

% Names of catchment and clusters
%--------------------------------------------------------------------------
IniCond.SITE = 'Pianello_in_Chiascio';
IniCond.FORCING = "ERA5Land";
IniCond.DeltaGMT= 1; % for Italy

% Name of the folder to save results 
%--------------------------------------------------------------------------
IniCond.run_folder = 'Run_44';

% Restart option
%--------------------------------------------------------------------------
restart.id = 0; % Set to 1 to continue a un-completed T&C run
%Define iter to restart
restart.month = 6;
restart.year = 2001;
restart.run = 'Run_41'; 

% Modelling period
%--------------------------------------------------------------------------
dateRun.start = "01-Jul-2001 00:00:00"; % Starting point of the simulation
dateRun.end = "5-Jul-2001 23:00:00"; % Last timestep of the simulation

% Folder with forcings
%--------------------------------------------------------------------------
%{
20251021: Pianelo in Chiascio
%}
forc_in = '20251021';

% Folder with Terrain inputs
%--------------------------------------------------------------------------
%{
1_Velino
2_Pianello_in_Chiascio
3_Tiber
%}
terrain_in = '2_Pianello_in_Chiascio';

%% DIRECTORIES
%==========================================================================
% All paths are here
%==========================================================================

% Main roots
%--------------------------------------------------------------------------
%Directories.root = '/nfs/scistore18/pelligrp/mrodrigu/'; %HPC
Directories.root = 'C:/Users/mrodrigu/Desktop/19_ISTA/1_Science_MegaWat/'; %Personal computer

% Sub-path for the model
%--------------------------------------------------------------------------
Directories.model = [Directories.root '1_TC/3_Model_Source/2_MegaWat/'];

% Sub-path for outputs
%--------------------------------------------------------------------------
Directories.save = [Directories.root '1_TC/3_Model_Source/2_MegaWat/2_distributed_version/3_Outputs/' IniCond.run_folder '/']; 

% Sub-path for restarting
%--------------------------------------------------------------------------
Directories.restart = [Directories.root '1_TC/3_Model_Source/2_MegaWat/2_distributed_version/3_Outputs/']; 

% Sub-path for forcings
%--------------------------------------------------------------------------
Directories.forcBC = [Directories.root '2_Forcing/3_Downscalling_ERA5/4_Bias_Corrected/' forc_in '/']; 
Directories.forcUN = [Directories.root '2_Forcing/3_Downscalling_ERA5/2_Output_Downscaling/' forc_in '/2_Units/']; 
Directories.forcPD = [Directories.root '2_Forcing/3_Downscalling_ERA5/2_Output_Downscaling/' forc_in '/1_Pure_Downscaling/']; 
Directories.forcRP = [Directories.root '2_Forcing/3_Downscalling_ERA5/5_Radiation_Partition/1_Bias_corrected/' forc_in '/']; 

% Vegetation parameters
%--------------------------------------------------------------------------
Directories.Vegpar = [Directories.root,'1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/3_Inputs/2_Apennine/Parameters_TC.xlsx']; 

% Preprocessing information
%--------------------------------------------------------------------------
Directories.PreProc = [Directories.root,'1_TC/3_Model_Source/2_MegaWat/4_Preparation_files/3_Distributed_GeoInputs_Apennines/8_Output/Pianello_in_Chiascio/20251024_250m/'];
dtm_file = 'dtm_Pianello_in_Chiascio_250m.mat';

% Dependencies
%--------------------------------------------------------------------------
addpath(genpath([Directories.model,'1_Functions'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([Directories.model,'5_Common_inputs'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([Directories.model,'3_Pyrenees_PointScale/2_Forcing'])); % Where is located the meteorological forcing and Shading matrix 
addpath(genpath([Directories.model,'3_Pyrenees_PointScale/3_Inputs'])); % Add path to Ca_Data


%% LOCATION OF OUTPUTS - CREATION OF FOLDERS
%==========================================================================
all_yy = [1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011];

if ~exist(Directories.save, 'dir') 
    disp('Folders do not exist for outputs. Creating new folders')
    mkdir(Directories.save);

    for i=1:length(all_yy)
    mkdir([Directories.save char(num2str(all_yy(i)))]);
    end  

    mkdir([Directories.save 'Initial']);
    mkdir([Directories.save 'Store']);

addpath(genpath(Directories.save)); 
end

folder_path = [Directories.root '1_TC/3_Model_Source/2_MegaWat']; % Put here the path of where you downloaded the repository



%% MODEL INPUTS
%==========================================================================

% Parameters for vegetation
%--------------------------------------------------------------------------
%{
Here are the parameters of the model for vegetation.
opts check the format of the columns. 
Use opts to force all the columns with values to be numeric.
%}
%--------------------------------------------------------------------------

opts = detectImportOptions([Directories.root '1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/3_Inputs/2_Apennine/Parameters_TC.xlsx']);
opts = setvartype(opts, [7:length(opts.VariableTypes)], 'double');

TT_par = readtable([Directories.root '1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/3_Inputs/2_Apennine/Parameters_TC.xlsx'], opts);

%% Some parameters

SITEs = IniCond.SITE; % All of my study catchments
FORCING = 'ERA5';
DeltaGMTs=[1];

% Catchment outlet POI name
outlet_names = ["Pianello_in_Chiascio"];

% Parameters
Tmods = [0]; % temperature modification above clean ice [°C];
Tmaxs = [2]; % Maximum air temperature for precipitation phase scheme (2-dual thresholds)
Tmins = [0]; % Minimum air temperature for precipitation phase scheme (2-dual thresholds)
 

%% Simulation period
x1=datetime(dateRun.start);
x2=datetime(dateRun.end);

%% OUTPUT MANAGER : decide which aggregated results to output
% Veg: Vegetation
% LC: Land Cover class

output_manag = [1,  % Catch_av
                0,  % Catch_std
                1,  % Veg_avg
                0,  % Veg_std
                0,  % LC_avg
                0,  % LC_std
                1,  % Maps
                1]; % POI

%% GLACIER AND AVALANCHE PARAMETERS
% Switch for glacier dynamics (keep to 0 for this case study)
Idyn = [0];

% switch on/off avalanching
Aval = [1];  % 1 to turn on avalanching, 0 to turn off avalanching
a_aval = [0.1]; % avalanche parameters a (Bernhart & Schulz 2010)
C_aval = [145]; % avalanche parameters C 

%% SNOW INITIAL CONDITION
% Initial snow depth
% Set at 2000 masl initially.
diff_IniSND = 1;
fn_IniSnowDepth = 'Cinca_Init_Snow_Depth_virtual.mat';
fn_IniSnowAlbedo = 'Cinca_Init_Snow_Albedo_virtual.mat';

%% Precipitation phase partitioning

%1 = 2-threshold, 2 = Ding 2017, 3 = single-threshold, 4 = Pomeroy 2013, 5
%= Wang 2019, 6 = Jennings 2018

parameterize_phase.OPT_Pr_Part = 2; % Choice of the precipitation phase scheme
parameterize_phase.Tmax = Tmaxs; % Upper air temperature for dual temperature threshold
parameterize_phase.Tmin = Tmins; % Lower air temperature for dual temperature threshold
parameterize_phase.Tconst = 2; % Air temperature for constant thresholds

parameterize_phase_labels = {'2-Ta','Ding','1-Ta','Pomeroy','Wang','Jennings'};
parameterize_phase_label = parameterize_phase_labels(parameterize_phase.OPT_Pr_Part);

% Skin layer thickness for the 2-layer snowpack module
hSTL = 0.003; %m

% Choice of the snow albedo scheme Evars_load
Albsno_method = 5; % 3 doesn't work, 4 is Brock 2000, 5 is Ding 2017

%% Prepare launcher variable based on chosen study site:
SITE = char(SITEs);

DeltaGMT= IniCond.DeltaGMT;

% Solar variables 
t_bef=[1]; % Integration interval for solar variables (hour)
t_aft=[0]; % Integration interval for solar variables (hour)

outlet_name = char(outlet_names);

%% Diplay setting of the incoming T&C model runs:
disp(['Site selected: ' IniCond.SITE])
disp(['Simulation period: ' datestr(x1,'dd-mmm-yyyy HH:MM') ' to ' datestr(x2,'dd-mmm-yyyy HH:MM')])
disp(['Precipitation phase scheme: ' parameterize_phase_label{:}])

if Idyn == 0; disp('Ice dynamics: off'); else; disp('Ice dynamics: on'); end
if Aval == 0; disp('Avalanching: off'); else; disp('Avalanching: on'); end

%% PATHS
%==========================================================================
% Attach folders, launch parallel pool, generate paths
%==========================================================================

outlocation = [Directories.save];

% Create output directory if it doesn't exist already
if ~exist(outlocation, 'dir');
    mkdir(outlocation);
    addpath(genpath(outlocation));
end


%% INPUTS
%==========================================================================
% Load geodata, time handling, carbon data
%==========================================================================
load([Directories.PreProc dtm_file]) %This file is in 2_nputs

% Central lat and lon Lat
%--------------------------------------------------------------------------
UTM_Y = y(floor(length(y)/2));
UTM_X = x(floor(length(x)/2));

[Lat, Lon] = utm2deg(UTM_X, UTM_Y, '33 T');

% Glacier set up
%--------------------------------------------------------------------------
if Idyn == 0
    GLH(GLH>0) = GLH(GLH>0) +400;  % Only use when ice dynamics are off, to avoid glacier disappearance
end

% PATCH, check pre-processing!
GLA_ID(MASK==0) = NaN;
GLA_ID(GLH==0) = NaN;

% Compute bedrock DEM for ice flow
DTM_Bedrock = DTM_orig-GLH; 


%% GENERAL PARAMETER
%==========================================================================

% Time set
%--------------------------------------------------------------------------
Date = x1:hours(1):x2;
N_time_step=length(Date);
Nd_time_step = ceil(N_time_step/24)+1;

% Steps
dt=3600; %% [s]
dth=1; %% [h]
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

% Carbon
%-------------------------------------------------------------------------
load(['Ca_Data.mat']);

d1 = find(abs(Date_CO2-datenum(Date(1)))<1/36);d2 = find(abs(Date_CO2-datenum(Date(end)))<1/36);
Ca=Ca(d1:d2);
clear d1 d2 Date_CO2
Oa= 210000;% Intercellular Partial Pressure Oxygen [umolO2/mol] -

% Solar variables
%--------------------------------------------------------------------------
L_day=zeros(length(Datam),1);
for j=2:24:length(Datam)
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(j,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
end

Lmax_day = max(L_day);
clear('h_S','delta_S','zeta_S','T_sunrise','T_sunset','L_day')

% Variables unknown for now
%--------------------------------------------------------------------------
a_dis=NaN; pow_dis=NaN;
clear a0 gam1 pow0 k2 DTii

% Other parameters
%--------------------------------------------------------------------------
rho_g = 0.35; %%% Spatial Albedo
ms_max = 10; %% Number of soil layers
md_max = 10; %% % Number of debris layers

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


%% TOPOGRAPHIC PARAMETER 
%==========================================================================
%{
Matrix [m_cell x n_cell]
DTM ;T_flow ; cellsize; xllcorner ; yllcorner ; SN ; outlet ; Aacc
%}
%==========================================================================

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
aTop= 1000*ones(m_cell,n_cell)*(cellsize^2)/cellsize; %% [mm] Area/Contour length ratio
Ared=ones(num_cell,1);

% Flow Boundary Condition
%--------------------------------------------------------------------------
Xout = [26]; % Temporal fix to extract discharge (this is the column)
Yout = [14]; % Tempotal fix to extract discharge (this is the row)

Slo_top(Youtlet,Xoutlet)=0.05;
npoint = length(Xout);
Area= (cellsize^2)*sum(sum(MASK)); %% Projected area [m^2]

% Flow potential
%--------------------------------------------------------------------------
T_pot=cell(1,ms_max);
for jk=1:ms_max
    T_pot{jk}= T_flow;
end

% Width channel
%--------------------------------------------------------------------------
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
Kres_Rock =8760; % Bedrock aquifer constant [h]
SPRINGn = SNn; % Spring Location


%% VEGETATION PARAMETERS LOOK-UP TABLE 
%==========================================================================
% For the Apennines there are 10 classes. These 10 classes are represented
% by the vector II. 
%==========================================================================

% Codes from VEG_CODE based on predefined classification of vegetation
%--------------------------------------------------------------------------
ksv=reshape(VEG_CODE,num_cell,1);

POI = table('Size', [10, 7], ...
          'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'double'}, ...
          'VariableNames', {'Class', 'Veg_type', 'Ccrowns','Cwat','Curb','Crock','Cbare'});

POI.Class = (1:10)';

% Representation of each vegetation class
POI.Veg_type = {["BLdec_lowElev_ro2" "grass_Bondone"],                                      % Class 1:  Decidious Broad-leaved forest
                ["grass_Bondone"],                                                          % Class 2:  Grassland/pasture
                ["Crops_WW"  "Crops_WB"  "Crops_S"  "Crops_R" "grass_Bondone"],             % Class 3:  Crops
                ["EvGreen_NeedLeaves_LeBray" "BLdec_highElev_Collelongo"],                  % Class 4:  Evergreen needleaves at high elevation (>700m)
                ["shrub_Balsablanca_A" "shrub_Balsablanca_B" "BLdec_highElev_Collelongo"],  % Class 5:  Mediterranean shrublands
                ["crops_Negrisia" "grass_Bondone"],                                         % Class 6:  Olives
                ["NoVeg"],                                                                  % Class 7:  Urban
                ["NoVeg"],                                                                  % Class 8:  Rock
                ["NoVeg"],                                                                  % Class 9:  Water
                ["grass_Alinya"]                                                            % Class 10: Bare soils
                };

POI.Ccrowns =  {[0.8 0.2],
                [0.9],
                [0.1 0.1 0.2 0.1 0.3],
                [0.5 0.2],
                [0.6 0.05 0.2],
                [0.4 0.2],
                [0.0],
                [0.0],
                [0.0],
                [0.1]};

% From class 1 to 10
%              1     2     3     4     5      6     7     8     9     10
POI.Cwat  = ({[0.0],[0.0],[0.0],[0.0],[0.0], [0.0],[1.0],[0.0],[1.0],[0.0]})';
POI.Curb  = ({[0.0],[0.0],[0.0],[0.0],[0.0], [0.0],[0.0],[0.0],[0.0],[0.0]})';
POI.Crock = ({[0.0],[0.0],[0.0],[0.0],[0.0], [0.0],[0.0],[1.0],[0.0],[0.0]})';
POI.Cbare = ({[0.0],[0.1],[0.2],[0.3],[0.15],[0.4],[0.0],[0.0],[0.0],[0.9]})';


% cc_max calculation
%z = 1;
cc_max = 1;
for z = 1:size(POI,1)
    if cc_max < length(POI.Ccrowns{z})
    cc_max = length(POI.Ccrowns{z});
    end
end

%% SPATIAL INDICES PER LAND COVER CLASS
Kinde            = find(MASK==1);  %%% basin index
idxCode.Veg1     = find(VEG_CODE == 1 & MASK == 1);
idxCode.Veg2     = find(VEG_CODE == 2 & MASK == 1); 
idxCode.Veg3     = find(VEG_CODE == 3 & MASK == 1); 
idxCode.Veg4     = find(VEG_CODE == 4 & MASK == 1); 
idxCode.Veg5     = find(VEG_CODE == 5 & MASK == 1); 
idxCode.Veg6     = find(VEG_CODE == 6 & MASK == 1); 
idxCode.Veg7     = find(VEG_CODE == 7 & MASK == 1); 
idxCode.Veg8     = find(VEG_CODE == 8 & MASK == 1);
idxCode.Veg9     = find(VEG_CODE == 9 & MASK == 1); 
idxCode.Veg10    = find(VEG_CODE == 10 & MASK == 1); 

%% SOIL PARAMETER 
%==========================================================================
%
%==========================================================================

%Pss = [800]; Pwp = [3500]; %% [kPa]
%Kfc = 0.2; %% [mm/h]
Phy = 10000; %% [kPa]
%[Osat,L,Pe,Ks,O33]=Soil_parameters_spatial(PSAN/100,PCLA/100,PORG/1000);
%[Ofc,Oss,Owp,Ohy]=Soil_parametersII_spatial(Osat,L,Pe,Ks,O33,Kfc,Pss,Pwp,Phy);
%plot soil properties
% Eventually we would like to store 'Ohy' and 'Osat'

%%
%clear Pss Pwp Kfc Phy L Pe Ks O33 Ofc Oss Owp
%Osat_OUT = Osat.*MASK; clear Osat
%Osat_OUT = reshape(Osat_OUT,num_cell,1);
%Ohy_OUT  = Ohy.*MASK; clear Ohy
%Ohy_OUT  = reshape(Ohy_OUT,num_cell,1);


%% Soil
PSAN=reshape(PSAN/100,num_cell,1);
PCLA=reshape(PCLA/100,num_cell,1);
PORG=reshape(PORG/100,num_cell,1);

%Zs_OUT=800*ones(num_cell,1);

%% SOLAR PARAMETER 
%==========================================================================
% HZ: Horizon angle array [angular degree]
% Z: Azimuth directions  [angular degree] from N
%==========================================================================

% Computation Horizon Angle
[HZ,Zasp] = Horizon_Angle(DTM_orig,cellsize);

% Sky View Factor and Terrain Configuration Factor
[SvF,Ct] = Sky_View_Factor(DTM_orig,atan(Slo_top)*180/pi,Aspect,HZ,Zasp);
%SvF = reshape(SvF,num_cell,1);

%% INITIAL CONDITION FOR SNOW ALBEDO
% Check if initial snow albedo is given

if ~exist('SNOWALB','var')
    SNOWALB = SNOWD;
    SNOWALB(SNOWD>0) = 0.6;
end 

%% If I am NOT restarting. So, in the normal case.
if restart.id == 0
out = [Directories.save 'Initial/INITIAL_CONDITIONS_' SITE '.mat'];
INIT_COND_v6(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Ca,SNOWD,SNOWALB,out, ...
   POI, TT_par, idxCode, Slo_top);
load(out);
    elseif restart.id == 1  %in case of restart  
    disp(['Parameters restarted from file: ' char(num2str(restart.month)) '_' char(num2str(restart.year)) '-' restart.run])    
    load([Directories.restart restart.run '/Store/Initial_Conditions_' char(num2str(restart.month)) '_' char(num2str(restart.year)) '.mat']);
end


%% NUMERICAL METHODS OPTIONS
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

% Time counter
%--------------------------------------------------------------------------
tic ;
%profile on

%% Workers - Only for personal computer (Windows)
%==========================================================================

%if I am working on my personal computer define this
%{

if ~contains(Directories.root,"nfs") 

    numWorkers = 6;
    
    % Create a parallel pool with the specified number of workers
    poolobj = gcp('nocreate'); % Check if a pool already exists
    if isempty(poolobj)
        parpool(numWorkers); % Create a new pool if one doesn't exist
    elseif poolobj.NumWorkers ~= numWorkers
        % If a pool exists with a different size, delete it and create a new one
        delete(poolobj);
        parpool(numWorkers);
    end

end
%}

%save([Directories.save 'Store/Before_t.mat'])


% Label for the creation of outputs
%--------------------------------------------------------------------------
output_creation = 0; 

%% Iterating on time
%==========================================================================
fts = 2;
for t=fts:N_time_step   
    
%waitbar(t/N_time_step,bau)
disp(['Iter: ' char(num2str(t))]);

% CHECK WELL WHY THE MODEL USES t-1 instead of t
Datam_S=Datam(t-1,:);

% Year and month to load the forcing
%----------------------------------------------------------------------
yy = char(num2str(Datam_S(1,1)));
mth = char(num2str(Datam_S(1,2)));

%% Matrix for storing - Initializing outputs
%==========================================================================
% Made by month
% Only for the first day of the month
% Because of this, the modelling must start on day 1 for now
%==========================================================================

if output_creation == 0
disp("Creating matrices for storing")

% function  Date_generator requires yy and mth as numbers
date_forMonth = Date_generator(str2num(yy),str2num(mth));
t_store = date_forMonth(end); % Date to store data

xdim = size(MASK, 1);
ydim = size(MASK, 2);
zdim = length(date_forMonth); 

ANPP_H_spatial         = single(zeros(xdim , ydim, zdim)); % Above ground net primary production
ANPP_L_spatial         = single(zeros(xdim , ydim, zdim));
An_H_spatial           = single(zeros(xdim , ydim, zdim)); % Net assimilation Vegetation
An_L_spatial           = single(zeros(xdim , ydim, zdim));
%B_H_spatial            = single(zeros(xdim , ydim, zdim)); % Carbon pool biomass
%B_L_spatial            = single(zeros(xdim , ydim, zdim));
%Ca_spatial             = single(zeros(xdim , ydim, zdim)); % CO2 atmospheric concentration
%Cicew_spatial          = single(zeros(xdim , ydim, zdim)); % Boolean operator for presence or absence of frozen water
%Cice_spatial           = single(zeros(xdim , ydim, zdim)); % Boolean operator for presence or absence of ice water
%Ci_shdH_spatial        = single(zeros(xdim , ydim, zdim)); % CO2 sunlit leaf internal concentration
%Ci_shdL_spatial        = single(zeros(xdim , ydim, zdim));
%Ci_sunH_spatial        = single(zeros(xdim , ydim, zdim));
%Ci_sunL_spatial        = single(zeros(xdim , ydim, zdim));
CK1_spatial            = single(zeros(xdim , ydim, zdim)); % Check on Mass Balance
Csnow_spatial          = single(zeros(xdim , ydim, zdim)); % Boolean operator for presence or absence of snow above frozen water 
Csno_spatial           = single(zeros(xdim , ydim, zdim)); % Boolean operator for presence or absence of snow 
dQ_S_spatial           = single(zeros(xdim , ydim, zdim)); % Residual from energy budget
DQ_S_spatial           = single(zeros(xdim , ydim, zdim)); % Residual of the energy budget
Dr_H_spatial           = single(zeros(xdim , ydim, zdim)); % Total Drainage from intercepted water
Dr_L_spatial           = single(zeros(xdim , ydim, zdim)); %
Ds_spatial             = single(zeros(xdim , ydim, zdim)); % Vapor Pressure Deficit
DT_S_spatial           = single(zeros(xdim , ydim, zdim)); % Residual temperature difference in the energy budget
dw_SNO_spatial         = single(zeros(xdim , ydim, zdim)); % Fraction of leaf covered by snow
ea_spatial             = single(zeros(xdim , ydim, zdim)); % Vapor Pressure
EG_spatial             = single(zeros(xdim , ydim, zdim)); % Evaporation from Bare soil
EICE_spatial           = single(zeros(xdim , ydim, zdim)); % Evaporation/sublimation from Ice
EIn_H_spatial          = single(zeros(xdim , ydim, zdim)); % Evaporation from intercepted water
EIn_L_spatial          = single(zeros(xdim , ydim, zdim));
EIn_rock_spatial       = single(zeros(xdim , ydim, zdim));
er_spatial             = single(zeros(xdim , ydim, zdim)); % Splash erosion
ESN_spatial            = single(zeros(xdim , ydim, zdim)); % Evaporation from the snowpack at the ground
SSN_spatial            = single(zeros(xdim , ydim, zdim));
EWAT_spatial           = single(zeros(xdim , ydim, zdim)); % Evaporation from water and ponds
FROCK_spatial          = single(zeros(xdim , ydim, zdim)); % Storage in fractured rocks
f_spatial              = single(zeros(xdim , ydim, zdim)); % Infiltration
Gfin_spatial           = single(zeros(xdim , ydim, zdim)); % Ground Heat Flux heat diffusion
GPP_H_spatial          = single(zeros(xdim , ydim, zdim));
GPP_L_spatial          = single(zeros(xdim , ydim, zdim));
G_spatial              = single(zeros(xdim , ydim, zdim)); % Ground Heat Flux force restore method
hc_H_spatial           = single(zeros(xdim , ydim, zdim)); % Vegetation Height
hc_L_spatial           = single(zeros(xdim , ydim, zdim));
H_spatial              = single(zeros(xdim , ydim, zdim));
ICE_D_spatial          = single(zeros(xdim , ydim, zdim)); % Ice thickness
ICE_spatial            = single(zeros(xdim , ydim, zdim)); % Ice water equivalent
Imelt_spatial          = single(zeros(xdim , ydim, zdim));
Inveg_spatial          = single(zeros(xdim , ydim, zdim));
In_rock_spatial        = single(zeros(xdim , ydim, zdim));
In_spatial             = single(zeros(xdim , ydim, zdim)); % Intercepted water (storage) 
In_SWE_spatial         = single(zeros(xdim , ydim, zdim)); % Intercepted snow water equivalent (storage)
IP_wc_spatial          = single(zeros(xdim , ydim, zdim)); % Ice pack water content
LAIdead_H_spatial      = single(zeros(xdim , ydim, zdim)); % Dead Leaf Area Index
LAIdead_L_spatial      = single(zeros(xdim , ydim, zdim));
LAI_H_spatial          = single(zeros(xdim , ydim, zdim)); % Leaf Area Index
LAI_L_spatial          = single(zeros(xdim , ydim, zdim));
Lk_rock_spatial        = single(zeros(xdim , ydim, zdim)); % Leakage rock surface to bedrock (recharge)
Lk_spatial             = single(zeros(xdim , ydim, zdim)); % Bottom Leakage soil to bedrock (recharge)
Lk_wat_spatial         = single(zeros(xdim , ydim, zdim)); % Leakage water pond to bedrock (recharge)
NICe_spatial           = single(zeros(xdim , ydim, zdim));
NPP_H_spatial          = single(zeros(xdim , ydim, zdim)); % Net Primary Production
NPP_L_spatial          = single(zeros(xdim , ydim, zdim));
NDVI_spatial           = single(zeros(xdim , ydim, zdim));
OF_spatial             = single(zeros(xdim , ydim, zdim)); % Soil Moisture First Soil Layer
OH_spatial             = single(zeros(xdim , ydim, zdim)); % Soil Moisture available to roots
OL_spatial             = single(zeros(xdim , ydim, zdim));
OS_spatial             = single(zeros(xdim , ydim, zdim)); % Soil Moisture for Bare Evaporation Layers
O_spatial              = single(zeros(xdim , ydim, zdim)); % Soil Moisture – Soil Water Content
PAR_spatial            = single(zeros(xdim , ydim, zdim));
%Pre_spatial            = single(zeros(xdim , ydim, zdim)); % Atmospheric Pressure
Pr_liq_spatial         = single(zeros(xdim , ydim, zdim)); % Liquid Precipitation
Pr_sno_spatial         = single(zeros(xdim , ydim, zdim)); % Solid (snow)
Pr_spatial             = single(zeros(xdim , ydim, zdim)); % Precipitation
QE_spatial             = single(zeros(xdim , ydim, zdim)); % Latent Heat
Qfm_spatial            = single(zeros(xdim , ydim, zdim)); % Heat for freezing or melting
Qlat_in_spatial        = single(zeros(xdim , ydim, zdim));
Qlat_out_spatial       = single(zeros(xdim , ydim, zdim));
Qv_spatial             = single(zeros(xdim , ydim, zdim)); % Heat advected by Precipitation
%Q_channel_spatial      = single(zeros(xdim , ydim, zdim));
%q_runon_spatial        = single(zeros(xdim , ydim, zdim)); % Runon
RA_H_spatial           = single(zeros(xdim , ydim, zdim)); % Autotrophic Respiration
RA_L_spatial           = single(zeros(xdim , ydim, zdim));
Rdark_H_spatial        = single(zeros(xdim , ydim, zdim)); % Leaf Dark Respiration
Rdark_L_spatial        = single(zeros(xdim , ydim, zdim));
Rd_spatial             = single(zeros(xdim , ydim, zdim)); % Saturation excess runoff
Rg_H_spatial           = single(zeros(xdim , ydim, zdim)); % Growth Respiration
Rg_L_spatial           = single(zeros(xdim , ydim, zdim));
Rh_spatial             = single(zeros(xdim , ydim, zdim)); % Infiltration excess runoff
Rmc_H_spatial          = single(zeros(xdim , ydim, zdim)); % Maintenance Respiration Carbohydrate reserve
Rmc_L_spatial          = single(zeros(xdim , ydim, zdim));
Rmr_H_spatial          = single(zeros(xdim , ydim, zdim)); % Maintenance Respiration roots
Rmr_L_spatial          = single(zeros(xdim , ydim, zdim));
Rms_H_spatial          = single(zeros(xdim , ydim, zdim)); % Maintenance Respiration sapwood
Rms_L_spatial          = single(zeros(xdim , ydim, zdim));
Rn_spatial             = single(zeros(xdim , ydim, zdim)); % Net radiation
ros_spatial            = single(zeros(xdim , ydim, zdim)); % Snow density
Rsw_spatial            = single(zeros(xdim , ydim, zdim));
SAI_H_spatial          = single(zeros(xdim , ydim, zdim)); % Stem Area Index
SAI_L_spatial          = single(zeros(xdim , ydim, zdim));
SAT_spatial            = single(zeros(xdim , ydim, zdim));
Smelt_spatial          = single(zeros(xdim , ydim, zdim));
SND_spatial            = single(zeros(xdim , ydim, zdim)); % Snow Depth
snow_albedo_spatial    = single(zeros(xdim , ydim, zdim));
SP_wc_spatial          = single(zeros(xdim , ydim, zdim)); % Snowpack water content
SWE_avalanched_spatial = single(zeros(xdim , ydim, zdim));
SWE_spatial            = single(zeros(xdim , ydim, zdim));
Ta_spatial             = single(zeros(xdim , ydim, zdim));
Tdamp_spatial          = single(zeros(xdim , ydim, zdim)); % Soil/snow Temperature at Dampening depth
Tdew_spatial           = single(zeros(xdim , ydim, zdim));
Tdp_spatial            = single(zeros(xdim , ydim, zdim)); % Soil Temperature of the layer
Tice_spatial           = single(zeros(xdim , ydim, zdim));
TsVEG_spatial          = single(zeros(xdim , ydim, zdim));
Ts_spatial             = single(zeros(xdim , ydim, zdim)); % Soil/snow Prognostic Temperature for the energy balance
T_H_spatial            = single(zeros(xdim , ydim, zdim)); % Transpiration
T_L_spatial            = single(zeros(xdim , ydim, zdim));
U_SWE_spatial          = single(zeros(xdim , ydim, zdim)); % Unloaded snow water equivalent from intercepted snow
Vice_spatial           = single(zeros(xdim , ydim, zdim));
V_spatial              = single(zeros(xdim , ydim, zdim)); % Volume of water stored in the soil layer
WAT_spatial            = single(zeros(xdim , ydim, zdim)); % Volume of water in the lakes/ponds
WIS_spatial            = single(zeros(xdim , ydim, zdim)); % Water flux incoming to the soil
WR_IP_spatial          = single(zeros(xdim , ydim, zdim)); % Water released from the ice pack
WR_SP_spatial          = single(zeros(xdim , ydim, zdim)); % Water released from the snow pack
%Ws_spatial             = single(zeros(xdim , ydim, zdim)); % Wind speed
ELitter_spatial        = single(zeros(xdim , ydim, zdim)); % ET litter
ESN_In_spatial         = single(zeros(xdim , ydim, zdim)); % ET litter
EIn_urb_spatial        = single(zeros(xdim , ydim, zdim)); % ET litter

% Series
QpointC_series         = single(zeros(zdim, length(Xout))); % Autotrophic Respiration

% Changing the label to not create it again at the next hour
%--------------------------------------------------------------------------
output_creation = 1;
end


    %% DISTRIBUTED FORCINGS
    %======================================================================
    % Forcing from stations and distributed across the catchment.
    % load spatial input data on first day and hour of month 
    %======================================================================      

    %% Stations

    % Load new forcings if the forcings do not exists or the year and
    % months are different from the ones loaded
    if ~exist('t2m', 'var') || str2num(yy) ~= year_loaded || str2num(mth) ~= month_loaded
    
    disp(['New forcing loaded for period: ' char(num2str(yy)) '-' char(num2str(mth)) ])    

    % Forcings from bias correction
    %----------------------------------------------------------------------
    t2m_file= [Directories.forcBC yy '/t2m_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    d2m_file= [Directories.forcBC yy '/d2m_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    tp_file= [Directories.forcBC yy '/tp_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    
    % Forcings from units
    %----------------------------------------------------------------------
    RH_file= [Directories.forcUN yy '/RH_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    strd_file= [Directories.forcUN yy '/strd_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository

    % Forcings from raw downscaling
    %----------------------------------------------------------------------
    sp_file= [Directories.forcPD yy '/sp_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    ws10_file= [Directories.forcPD yy '/ws10_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository

    % Forcings from Radiation Partition
    %----------------------------------------------------------------------
    N_file= [Directories.forcRP yy '/N_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    PARB_file= [Directories.forcRP yy '/PARB_3D_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    PARD_file= [Directories.forcRP yy '/PARD_3D_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    SAB1_file= [Directories.forcRP yy '/SAB1_3D_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    SAB2_file= [Directories.forcRP yy '/SAB2_3D_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    SAD1_file= [Directories.forcRP yy '/SAD1_3D_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository
    SAD2_file= [Directories.forcRP yy '/SAD2_3D_' yy '_' mth '.mat']; % Put here the path of where you downloaded the repository

    % Opening files
    %----------------------------------------------------------------------
    load(t2m_file);    load(d2m_file);    load(tp_file);     load(RH_file);
    load(sp_file);     load(ws10_file);   load(N_file);      load(PARB_file);
    load(PARD_file);   load(SAB1_file);   load(SAB2_file);   load(SAD1_file);
    load(SAD2_file);   load(strd_file);

    %forc(forc.t2m<forc.d2m,:)
    
    % Decompressing - Only for big data sets
    %--------------------------------------------------------------------------
    t2m        = Data_compressor_grid(t2m_3D_bc_c,"t2m","decompress");
    d2m        = Data_compressor_grid(d2m_3D_bc_c,"d2m","decompress");
    tp         = Data_compressor_grid(tp_3D_bc_c,"t2m","decompress");
    RH         = Data_compressor_grid(RH_3D_c,"RH","decompress");
    sp         = Data_compressor_grid(sp_3D_c,"sp","decompress");
    ws10       = Data_compressor_grid(ws10_3D_c,"ws10","decompress");
    strd       = Data_compressor_grid(strd_3D_c,"strd","decompress");
    PARB_model = Data_compressor_grid(PARB,"PARB","decompress");
    PARD_model = Data_compressor_grid(PARD,"PARD","decompress");
    SAB1_model = Data_compressor_grid(SAB1,"SAB1","decompress");
    SAB2_model = Data_compressor_grid(SAB2,"SAB2","decompress");
    SAD1_model = Data_compressor_grid(SAD1,"SAD1","decompress");
    SAD2_model = Data_compressor_grid(SAD2,"SAD2","decompress");
   
    % Store year and month loaded
    %----------------------------------------------------------------------
    year_loaded = str2num(yy);
    month_loaded = str2num(mth);

    end
    
    %% Finding the row in the forcing of the modeling date      
    t_forc = find(Date(t-1) == date_forMonth);

    % Other parameters
    %----------------------------------------------------------------------
    t_bef=1; t_aft=0; % otherwise problems when loading SWPART        
    
    %% Temperature      
    % Using stations to obtain distributed temperature along the catchment
    Ta_S = t2m(:,:,t_forc);
    Ta_S = double(reshape(Ta_S,num_cell,1));
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
    N_S = strd(:,:,t_forc);
    N_S = double(reshape(N_S,num_cell,1));
    N_S(MASKn==0) = NaN;

    %% Relative humidity (ERA5)
    U_S = RH(:,:,t_forc)/100;
    U_S = double(reshape(U_S,num_cell,1));
    U_S(MASKn==0) = NaN;

    %% Wind speed
    Ws_S = ws10(:,:,t_forc);
    Ws_S = double(reshape(Ws_S,num_cell,1));
    Ws_S(Ws_S < 0.01) = 0.01;
    Ws_S(MASKn==0) = NaN;

    %% Precipitation
    Pr_S = tp(:,:,t_forc);
    Pr_S = double(reshape(Pr_S,num_cell,1));
    Pr_S(isnan(Pr_S)) = 0;

    %% Air pressure
    Pre_S = sp(:,:,t_forc)/100;
    Pre_S = double(reshape(Pre_S,num_cell,1));
    

    % Reshape - This because of the problem in the Hydrological module
    %----------------------------------------------------------------------
    Slo_top2 = reshape(Slo_top,num_cell,1); % Creating aux variable for the hydrological module
    aTop = reshape(aTop,num_cell,1);
    %Slo_top = reshape( Slo_top,num_cell,1);

    %% Radiation Partition
    
    [jDay]= JULIAN_DAY(Datam_S);
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day]= SetSunVariables(Datam_S,DeltaGMT,Lon,Lat,t_bef,t_aft);
    [ShF] = Shadow_Effect(DTM,h_S,zeta_S,HZ,Zasp);

    %needed, if terrain effects have not been considered during pre-processing
    cos_fst = cos(atan(Slo_top))*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;        

    SAD1_S = SAD1_model(:,:,t_forc);
    SAD2_S = SAD2_model(:,:,t_forc);
    SAB1_S = SAB1_model(:,:,t_forc);
    SAB2_S = SAB2_model(:,:,t_forc);
    PARD_S = PARB_model(:,:,t_forc);
    PARB_S = PARD_model(:,:,t_forc);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAD1_S = SAD1_S.*SvF + Ct.*rho_g.*(SAD1_S./sin(h_S).*cos_fst + (1-SvF).*SAD1_S);
    SAD2_S = SAD2_S.*SvF + Ct.*rho_g.*(SAD2_S./sin(h_S).*cos_fst + (1-SvF).*SAD2_S);
    PARD_S = PARD_S.*SvF + Ct.*rho_g.*(PARD_S./sin(h_S).*cos_fst + (1-SvF).*PARD_S);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAB1_S =(SAB1_S./sin(h_S)).*cos_fst.*ShF;
    SAB2_S =(SAB2_S./sin(h_S)).*cos_fst.*ShF;
    PARB_S = (PARB_S./sin(h_S)).*cos_fst.*ShF;        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

    SAB1_S  =  double(reshape(SAB1_S,num_cell,1));
    SAB2_S  =  double(reshape(SAB2_S,num_cell,1));
    SAD1_S  =  double(reshape(SAD1_S,num_cell,1));
    SAD2_S  =  double(reshape(SAD2_S,num_cell,1));
    PARB_S  =  double(reshape(PARB_S,num_cell,1));
    PARD_S  =  double(reshape(PARD_S,num_cell,1));
  
       
    %% Carbon
    Ca_S = Ca(t)*MASKn;
    IrD_S =  MASK*0; IrD_S=reshape(IrD_S,num_cell,1);
    Salt_S =  MASK*0; Salt_S=reshape(Salt_S,num_cell,1);
    %N_S = N(t)*MASKn;
    
    %% Vapor pressure
    % esat/ea/Ds/Tdew maps
    a=17.27; b=237.3;
    esat_S=611*exp(a*Ta_S./(b+Ta_S)); % Vapour pressure at saturation (Pa)
    ea_S=U_S.*esat_S;                 % Vapour pressure (Pa)
    Ds_S= esat_S - ea_S;              % Vapor Pressure Deficit (Pa)
    Ds_S(Ds_S<0)=0; 
    xr=a*Ta_S./(b+Ta_S)+log10(U_S);
    Tdew_S=b*xr./(a-xr);              % Presumed dewpoint temperature (C)
    clear a b xr;
    
    esat_S =reshape(esat_S,num_cell,1);
    ea_S   =reshape(ea_S,num_cell,1);
    Ds_S   =reshape(Ds_S,num_cell,1);
    Tdew_S =reshape(Tdew_S,num_cell,1);
      
    %{
    forcing.t2m = t2m;
    forcing.d2m = d2m;
    forcing.tp = tp;
    forcing.ws10 = ws10;
    forcing.sp = sp; 

    check_var(forcing, ...
    ["t2m" ... % Temperature
    "d2m" ...  % Dew Point temperature
    "tp" ...   % Precipitation
    "ssrd" ... % Downward short wave radiation
    "strd" ... % Downdward Long wave radiation
    "ws10" ... % Wind speed
    "sp" ...   % Air pressure
    "es" ...   % saturation vapor pressure
    "ea" ...   % actual vapor pressure
    "RH" ...   % Relative humidity
    "SAD1" ... % SAD1
    "SAD2" ... % SAD2
    "SAB1" ... % SAB1
    "SAB2" ... % SAB2
    "PARB" ... % PARB
    "PARD" ... % PARD
    "N" ...    % Cloudiness
    ],Point)
    %}
    
    %% SPATIAL INITIALIZATION VECTOR PREDEFINING
    %======================================================================
   
    if t == 2
        % General vegetation/hydrology
        %------------------------------------------------------------------
        alp_soil=	     alp_soiltm1;
        b_soil=          b_soiltm1;
        Bam=             Bamtm1;
        Bem=             Bemtm1;
        BLit=            BLittm1;
        Ccrown_t=	     Ccrown_t_tm1;
        Cice=            Cicetm1;
        Cicew=           Cicewtm1;
        CK1=             CK1tm1;
        Csno=            Csnotm1;
        Csnow=           Csnowtm1;
        dQ_S=            dQ_Stm1;
        DQ_S=            DQ_Stm1;
        dQVEG=           dQVEGtm1;
        DT_S=            DT_Stm1;
        dw_SNO=          dw_SNOtm1;
        e_sno=           e_snotm1;
        EG=              EGtm1;
        EICE=            EICEtm1;
        EIn_rock=	     EIn_rocktm1;
        EIn_urb=	     EIn_urbtm1;
        EK=              EKtm1;
        ELitter=	     ELittertm1;
        er=              ertm1;
        ESN_In=          ESN_Intm1;
        SSN_In =         SSN_Intm1;
        ESN=             ESNtm1; 
        SSN=             SSNtm1;  
        EWAT=            EWATtm1;
        FROCK=           FROCKtm1;
        f=               ftm1;
        Gfin=            Gfintm1;
        G=               Gtm1;
        H=               Htm1;
        HV=              HVtm1;
        ICE_D=           ICE_Dtm1;
        ICE=             ICEtm1;
        Imelt=           Imelttm1;
        In_H=            In_Htm1;
        In_Litter=	     In_Littertm1;
        In_L=            In_Ltm1;
        In_rock=	     In_rocktm1;
        In_SWE=          In_SWEtm1;
        In_urb=          In_urbtm1;
        IP_wc=           IP_wctm1;
        Lk_rock=	     Lk_rocktm1;
        Lk_wat=          Lk_wattm1;
        Lk=              Lktm1;
        Lpho=            Lphotm1;
        NavlI=           NavlItm1;
        NDVI=            NDVItm1;
        NIce=            NIcetm1;
        NIn_SWE=	     NIn_SWEtm1;
        OF=              OFtm1;
        Oice=            Oicetm1;
        OS=              OStm1;
        O=               Otm1;
        POT=             POTtm1;
        Pr_liq=          Pr_liqtm1;
        Pr_sno=          Pr_snotm1;
        Q_channel=	     Q_channel;
        Q_exit=          Q_exit;
        q_runon=         q_runon;
        QE=              QEtm1;
        QEV=             QEVtm1;
        Qfm=             Qfmtm1;
        Qi_in=           Qi_in;
        %Qi_in_Ro=        Qi_out;
        Qi_out_Rout=     Qi_out_Rout;
        Qi_out=          Qi_outtm1;
        Qsub_exit=	     Qsub_exit;
        Qv=              Qvtm1;
        r_litter=	     r_littertm1;
        r_soil=          r_soiltm1;
        ra=              ratm1;
        Rd=              Rdtm1;
        Rh=              Rhtm1;
        Rn=              Rntm1;
        ros=             rostm1;
        SE_rock=	     SE_rocktm1;
        SE_urb=          SE_urbtm1;
        Slo_head=	     Slo_head;
        Smelt=           Smelttm1;
        SND=             SNDtm1;
        snow_albedo=     snow_albedotm1;
        soil_albedo=     soil_albedotm1;
        SP_wc=           SP_wctm1;
        surface_albedo=  surface_albedotm1;
        SWE=             SWEtm1;
        SWE_avalanched=  SWE_avalanchedtm1;
        t_sls=           t_slstm1;
        tau_sno=	     tau_snotm1;
        Tdamp=           Tdamptm1;
        Tdeb=            Tdebtm1;
        Tdp=             Tdptm1;
        Tdp_snow =       Tdp_snowtm1;
        Tice=            Ticetm1;
        Tstm0=           Tstm0;
        Ts=              Tstm1;
        Ts_under =       Ts_undertm1;
        TsVEG=           TsVEGtm1;
        U_SWE=           U_SWEtm1;
        Vice=            Vicetm1;
        V=               Vtm1;
        WAT=             WATtm1;
        WIS=             WIStm1;
        WR_IP=           WR_IPtm1;
        WR_SP=           WR_SPtm1;
        Ws_under=	     Ws_undertm1;
        ZWT=             ZWTtm1;
        
        % Specifications for high & low vegetation
        %------------------------------------------------------------------
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
    %P = (Xout - 1) * 133 + Yout;

    %Follow=MASK; %Debugging
    parfor ij= 1:num_cell % this is a parfor
        % Good practice to use a simple for loop for debugging/testing
        % ij is the index to go pixel by pixel through the mask
        % ij=1:num_cell
        %disp(strcat('in the loop', ij))
        
        %% Debugging for a for loop
        %Follow(ij) = 2222; %22307 - 22308 (row 107 - col 151)
        %disp('bye')
        %{
        if ismember(ij, [num_cell/8, num_cell/4, num_cell/2, 3*num_cell/4])   
        disp(ij)
        end        
        %Ta,Ts,Pre,zatm,disp_h,zom,zoh,Ws,ea
        %} 
        % =================================================================
        
        if MASKn(ij)== 1
            
            %disp(['Cell: ' char(num2str(ij)) ', Veg Type: ' char(num2str(ksv(ij)))])
            Elev=DTM(ij);
            %[i,j] = ind2sub([m_cell,n_cell],ij);

            % BOUNDARY CONDITION  
            % INTRODUCED SOIL AND VEG. for ij
            %--------------------------------------------------------------
            [aR,             Zs,             EvL_Zs,       Inf_Zs,     Bio_Zs,      Zinf, ...
             RfH_Zs,         RfL_Zs,         dz,           Ks_Zs,      Dz,          ms, ...
             Kbot,           Krock,          zatm,         Ccrown,     Cbare,       Crock, ...
             Curb,           Cwat,           Color_Class,  OM_H,       OM_L,        PFT_opt_H, ...
             PFT_opt_L,      d_leaf_H,       d_leaf_L,     SPAR,       Phy,         Soil_Param, ...
             Interc_Param,   SnowIce_Param,  VegH_Param,   VegL_Param, fpr,         VegH_Param_Dyn, ...
             VegL_Param_Dyn, Stoich_H,       aSE_H,        Stoich_L,   aSE_L,       fab_H,...
             fbe_H,          fab_L,          fbe_L,        ZR95_H,     ZR95_L,      In_max_urb,...
             In_max_rock,    K_usle,         Urb_Par,      Deb_Par,    Zs_deb,      Sllit,...
             Kct,            ExEM,           ParEx_H,      Mpar_H,     ParEx_L,     Mpar_L,] = ....
                                    PARAMETERS_SOIL2( ...
             ksv(ij),        PSAN(ij),       PCLA(ij),     PORG(ij),   DEB_MAP(ij),  md_max,...
             Afirn(ij),      SOIL_TH(ij),    POI,        TT_par);
             
            
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
              
               
        %% VEGETATION MODULE
        %==================================================================
        % FUNCTION: VEGETATION_MODULE_PAR
        %==================================================================

       [LAI_H(ij,:),          B_H(ij,:,:),        NPP_H(ij,:),          ANPP_H(ij,:),        Rg_H(ij,:), ...
        RA_H(ij,:),           Rms_H(ij,:),        Rmr_H(ij,:),          Rmc_H(ij,:),         PHE_S_H(ij,:),...
        dflo_H(ij,:),         AgeL_H(ij,:),       e_rel_H(ij,:),        e_relN_H(ij,:),      LAI_L(ij,:),...
        B_L(ij,:,:),          NPP_L(ij,:),        ANPP_L(ij,:),         Rg_L(ij,:),          RA_L(ij,:),...
        Rms_L(ij,:),          Rmr_L(ij,:),        Rmc_L(ij,:),          PHE_S_L(ij,:),       dflo_L(ij,:), ...
        AgeL_L(ij,:),         e_rel_L(ij,:),      e_relN_L(ij,:),       SAI_H(ij,:),         hc_H(ij,:), ...
        SAI_L(ij,:),          hc_L(ij,:),         LAIdead_H(ij,:),      NBLeaf_H(ij,:),      Sr_H(ij,:), ...
        Slf_H(ij,:),          Sfr_H(ij,:),        Sll_H(ij,:),          Swm_H(ij,:),         Rexmy_H(ij,:,:), ...
        NupI_H(ij,:,:),       NuLit_H(ij,:,:),    LAIdead_L(ij,:),      NBLeaf_L(ij,:),      Sr_L(ij,:), ...
        Slf_L(ij,:),          Sfr_L(ij,:),        Sll_L(ij,:),          Swm_L(ij,:),         Rexmy_L(ij,:,:),... 
        NupI_L(ij,:,:),       NuLit_L(ij,:,:),    Rrootl_H(ij,:),       AgeDL_H(ij,:),       Bfac_dayH(ij,:), ...
        Bfac_weekH(ij,:),     NPPI_H(ij,:),       TdpI_H(ij,:),         PARI_H(ij,:,:),      NBLI_H(ij,:), ...
        RB_H(ij,:,:),         FNC_H(ij,:),        Nreserve_H(ij,:),     Preserve_H(ij,:),    Kreserve_H(ij,:), ...
        rNc_H(ij,:),          rPc_H(ij,:),        rKc_H(ij,:),          ManIH(ij,:),         Rrootl_L(ij,:), ...
        AgeDL_L(ij,:),        Bfac_dayL(ij,:),    Bfac_weekL(ij,:),     NPPI_L(ij,:),        TdpI_L(ij,:), ...
        PARI_L(ij,:,:),       NBLI_L(ij,:),       RB_L(ij,:,:),         FNC_L(ij,:),         Nreserve_L(ij,:), ...
        Preserve_L(ij,:),     Kreserve_L(ij,:),   rNc_L(ij,:),          rPc_L(ij,:),         rKc_L(ij,:), ...
        ManIL(ij,:),          TexC_H(ij,:),       TexN_H(ij,:),         TexP_H(ij,:),        TexK_H(ij,:), ...
        TNIT_H(ij,:),         TPHO_H(ij,:),       TPOT_H(ij,:),         SupN_H(ij,:),        SupP_H(ij,:), ...
        SupK_H(ij,:),         ISOIL_H(ij,:,:),    TexC_L(ij,:),         TexN_L(ij,:),        TexP_L(ij,:), ...
        TexK_L(ij,:),         TNIT_L(ij,:),       TPHO_L(ij,:),         TPOT_L(ij,:),        SupN_L(ij,:), ...
        SupP_L(ij,:),         SupK_L(ij,:),       ISOIL_L(ij,:,:),      BA_H(ij,:),          Tden_H(ij,:), ...
        AgePl_H(ij,:),        BA_L(ij,:),         Tden_L(ij,:),         AgePl_L(ij,:),       Ccrown_t(ij,:)]= ...
                          VEGETATION_MODULE_PAR( ...
        cc_max,               Ccrown,              ZR95_H,              ZR95_L,              B_Htm1(ij,:,:),...
        PHE_S_Htm1(ij,:),     dflo_Htm1(ij,:),     AgeL_Htm1(ij,:),     AgeDL_Htm1(ij,:),    Ta_t(ij,:), ...
        PAR_t(ij,:),          Tdp_H_t(ij,:,:),     Psi_x_H_t(ij,:,:),   Psi_l_H_t(ij,:,:),   An_H_t(ij,:,:), ...
        Rdark_H_t(ij,:,:),    NPP_Htm1(ij,:),      jDay,                Datam_S,             NPPI_Htm1(ij,:), ...
        TdpI_Htm1(ij,:),      Bfac_weekHtm1(ij,:), Stoich_H,            aSE_H,               VegH_Param_Dyn,...
        Nreserve_Htm1(ij,:),  Preserve_Htm1(ij,:), Kreserve_Htm1(ij,:), Nuptake_H(ij,:),     Puptake_H(ij,:), ...
        Kuptake_H(ij,:),      FNC_Htm1(ij,:),      Tden_Htm1(ij,:),     AgePl_Htm1(ij,:),    fab_H, ...
        fbe_H,                ParEx_H,             Mpar_H,              TBio_H(ij,:),        SAI_Htm1(ij,:), ...
        hc_Htm1(ij,:),        B_Ltm1(ij,:,:),      PHE_S_Ltm1(ij,:),    dflo_Ltm1(ij,:),     AgeL_Ltm1(ij,:), ...
        AgeDL_Ltm1(ij,:),     Tdp_L_t(ij,:,:),     Psi_x_L_t(ij,:,:),   Psi_l_L_t(ij,:,:),   An_L_t(ij,:,:), ...
        Rdark_L_t(ij,:,:),    NPP_Ltm1(ij,:),      NPPI_Ltm1(ij,:),     TdpI_Ltm1(ij,:),     Bfac_weekLtm1(ij,:),...
        NupI_Htm1(ij,:,:),    NupI_Ltm1(ij,:,:),   NuLit_Htm1(ij,:,:),  NuLit_Ltm1(ij,:,:),  NBLeaf_Htm1(ij,:), ...
        NBLeaf_Ltm1(ij,:),    PARI_Htm1(ij,:,:),   NBLI_Htm1(ij,:),     PARI_Ltm1(ij,:,:),   NBLI_Ltm1(ij,:),...
        Stoich_L,             aSE_L,               VegL_Param_Dyn,      NavlI(ij,:),         Bam(ij), ...
        Bem(ij),              Ccrown_t_tm1(ij,:),  Nreserve_Ltm1(ij,:), Preserve_Ltm1(ij,:), Kreserve_Ltm1(ij,:), ...
        Nuptake_L(ij,:),      Puptake_L(ij,:),     Kuptake_L(ij,:),     FNC_Ltm1(ij,:),      Tden_Ltm1(ij,:), ...
        AgePl_Ltm1(ij,:),     fab_L,               fbe_L,               ParEx_L,             Mpar_L, ...
        TBio_L(ij,:),         SAI_Ltm1(ij,:),      hc_Ltm1(ij,:),       ExEM,                Lmax_day, ...
        L_day,                Se_bio,              Tdp_bio,             OPT_EnvLimitGrowth,  OPT_VD, ...
        OPT_VCA,              OPT_ALLOME,          OPT_SoilBiogeochemistry);
    
                BLit(ij,:)= 0.0 ; % %% %%[kg DM / m2]
            end        

            %% HYDROLOGY MODULE
            %==============================================================
            % FUNTION: HYDROLOGY_MODULE_PAR
            %==============================================================                                                 
    %try
 %{   
mm = {Vtm1,        Oicetm1,       aR,           Zs,                  EvL_Zs,        ...
     Inf_Zs, ...
     Zinf,         RfH_Zs,        RfL_Zs,       dz,                  Dz, ...
     ms,           Kbot,          Pr_S,         Ta_S,                Ds_S, ...
     Ws_S,         zatm,          Tstm1,        dt,                  dth, ...
     ea_S,         N_S,           Pre_S,        Tstm0,               LAI_H, ...
     SAI_H,        LAI_L,         SAI_L,        LAIdead_H,           LAIdead_L,...
     Rrootl_H,     Rrootl_L,      BLit,         Sllit,               Kct,...
     Datam_S,      DeltaGMT,      Lon,          Lat,                 t_bef, ...
     t_aft,        Ccrown,        Cbare,        Crock,               Curb, ...
     Cwat,         SAB1_S,        SAB2_S,       SAD1_S,              SAD2_S, ...
     PARB_S,       PARD_S,        SvF,          SNDtm1,              snow_albedotm1, ...
     Color_Class,  OM_H,          OM_L,         PFT_opt_H,           PFT_opt_L, ...
     hc_H,         hc_L,          d_leaf_H,     d_leaf_L,            Soil_Param, ...
     Interc_Param, SnowIce_Param, VegH_Param,   VegL_Param,          Ca_S, ...
     Oa,           Citm1_sunH,    Citm1_shdH,   Citm1_sunL,          Citm1_shdL,...
     e_rel_H,      e_relN_H,      e_rel_L,      e_relN_L,            e_snotm1, ...
     In_Htm1,      In_Ltm1,       In_Littertm1, In_urbtm1,           In_rocktm1(ij), ...
     SWEtm1,       In_SWEtm1,     Tdebtm1,      Ticetm1,             Tdptm1(ij,:), ...
     Tdp_snowtm1,  Tdamptm1,      Ts_undertm1,  WATtm1,              ICEtm1(ij), ...
     IP_wctm1,     ICE_Dtm1,      Cicewtm1,     Vx_Htm1,             Vl_Htm1(ij,:), ...
     Vx_Ltm1,      Vl_Ltm1,       Psi_x_Htm1,   Psi_l_Htm1,          Psi_x_Ltm1(ij,:), ...
     Psi_l_Ltm1,   ZR95_H,        ZR95_L,       FROCKtm1,            Krock,...
     Urb_Par, ...
     Deb_Par,      Zs_deb,        Tdew_S,       t_slstm1,            rostm1, ...
     SP_wctm1,     fpr,           IrD_S,        In_max_urb,          In_max_rock, ...
     K_usle,       tau_snotm1,    Ta_day,       Slo_top2,            Slo_head, ...
     Asur,         Ared,          aTop,         EKtm1,               q_runon, ...
     Qi_in,        Ws_undertm1,   Pr_sno_t,     pow_dis,             a_dis, ...
     Salt_S,       SPAR,SNn,      OPT_min_SPD,  OPT_VegSnow,         OPT_SoilTemp, ...
     OPT_PlantHydr,Opt_CR,        Opt_ST,       Opt_ST2,             OPT_SM, ...
     OPT_STh,      OPT_FR_SOIL,   OPT_PH,       parameterize_phase,  hSTL, ...
     Albsno_method};

sizes = cellfun(@size, mm, 'UniformOutput', false);
stringCell = cellfun(@(x) sprintf('[%.0f,%.0f]', x(1), x(2)), sizes, 'UniformOutput', false);

T = table('Size', [31,5], ...
    'VariableTypes', {'string', 'string', 'string', 'string', 'string'});

%i = 5;
z = 1;
y = 1;
        for i = 1:length(sizes)
        T(z,y) = stringCell(i);
        y = y+1;
            if mod(i, 5) == 0ij
            z = z+1;
            y=1;
            end
        
        end
 
    if t == 2    
    writetable(T, [Directories.save ' output_data.csv']);
    end
    
 %}

    [V(ij,:),           O(ij,:),          Vice(ij,:),       Oice(ij,:),          ZWT(ij), ...         % 1
     OF(ij),            OS(ij),           OH(ij,:),         OL(ij,:),            Psi_s_H(ij,:),...    % 2
     Psi_s_L(ij,:),     Rd(ij),           Qi_out(ij,:),     Rh(ij),              Lk(ij), ...          % 3
     f(ij),             WIS(ij),          Ts(ij),           Csno(ij),            Cice(ij), ...        % 4 
     NDVI(ij),          Pr_sno(ij),       Pr_liq(ij),       rb_H(ij,:),          rb_L(ij,:), ...      % 5 
     rs_sunH(ij,:),     rs_sunL(ij,:),    rs_shdH(ij,:),    rs_shdL(ij,:),       rap_H(ij,:), ...     % 6
     rap_L(ij,:),       r_soil(ij),       b_soil(ij),       alp_soil(ij),        ra(ij), ...          % 7
     r_litter(ij,:),    WR_SP(ij),        U_SWE(ij),        NIn_SWE(ij),         dQ_S(ij), ...        % 8 
     DQ_S(ij),          DT_S(ij),         WAT(ij),          ICE(ij),             ICE_D(ij), ...       % 9
     IP_wc(ij),         WR_IP(ij),        NIce(ij),         Cicew(ij),           Csnow(ij), ...       % 10
     FROCK(ij),         Dr_H(ij,:),       Dr_L(ij,:),       SE_rock(ij),         SE_urb(ij), ...      % 11
     Lk_wat(ij),        Lk_rock(ij),      An_L(ij,:),       An_H(ij,:),          Rdark_L(ij,:), ...   % 12
     Rdark_H(ij,:),     Ci_sunH(ij,:),    Ci_sunL(ij,:),    Ci_shdH(ij,:),       Ci_shdL(ij,:), ...   % 13
     Rn(ij),            H(ij),            QE(ij),           Qv(ij),              Lpho(ij), ...        % 14
     T_H(ij,:),         T_L(ij,:),        EIn_H(ij,:),      EIn_L(ij,:),         EG(ij), ...          % 15
     ELitter(ij),       ESN(ij),          ESN_In(ij),       EWAT(ij),            EICE(ij), ...        % 16
     EIn_urb(ij),       EIn_rock(ij),     dw_SNO(ij),       Imelt(ij),           Smelt(ij), ...       % 17
     G(ij),             Gfin(ij),         Tdp(ij,:),        Tdp_snow(ij,:),      Tdeb(ij,:), ...      % 18
     Tice(ij),          Tdamp(ij),        Tdp_H(ij,:),      Tdp_L(ij,:),         SWE(ij), ...         % 19
     SND(ij),           ros(ij),          In_SWE(ij),       SP_wc(ij),           Qfm(ij), ...         % 20 
     t_sls(ij),         In_H(ij,:),       In_L(ij,:),       In_Litter(ij),       In_urb(ij), ...      % 21
     In_rock(ij),       gsr_H(ij,:),      Psi_x_H(ij,:),    Psi_l_H(ij,:),       Jsx_H(ij,:), ...     % 22
     Jxl_H(ij,:),       Kleaf_H(ij,:),    Kx_H(ij,:),       Vx_H(ij,:),          Vl_H(ij,:),...       % 23
     gsr_L(ij,:),       Psi_x_L(ij,:),    Psi_l_L(ij,:),    Jsx_L(ij,:),         Jxl_L(ij,:), ...     % 24
     Kleaf_L(ij,:),     Kx_L(ij,:),       Vx_L(ij,:),       Vl_L(ij,:),          fapar_H(ij,:), ...   % 25
     fapar_L(ij,:),     SIF_H(ij,:),      SIF_L(ij,:),      Ws_under(ij),        er(ij), ...          % 26
     snow_albedo(ij,:), tau_sno(ij),      e_sno(ij),        HV(ij),              QEV(ij), ...         % 27
     dQVEG(ij),         TsVEG(ij),        Ts_under(ij),     EK(ij),              POT(ij,:), ...       % 28
     CK1(ij)] = ...
                HYDROLOGY_MODULE_PAR(...
     Vtm1(ij,:),            Oicetm1(ij,:),       aR,               Zs,                  EvL_Zs, ...
     Inf_Zs,                Zinf,                RfH_Zs,           RfL_Zs,              dz,  ...
     Dz,                    ms,                  Kbot,             Pr_S(ij),            Ta_S(ij), ...
     Ds_S(ij),              Ws_S(ij),            zatm,             Tstm1(ij),           dt, ...
     dth,                   ea_S(ij),            N_S(ij),          Pre_S(ij),           Tstm0(ij), ...
     LAI_H(ij,:),           SAI_H(ij,:),         LAI_L(ij,:),      SAI_L(ij,:),         LAIdead_H(ij,:), ...
     LAIdead_L(ij,:),       Rrootl_H(ij,:),      Rrootl_L(ij,:),   BLit(ij,:),          Sllit,   ...
     Kct,                   Datam_S,             DeltaGMT,         Lon,                 Lat,     ...
     t_bef,                 t_aft,               Ccrown,           Cbare,               Crock,  ...
     Curb,                  Cwat,                SAB1_S(ij),       SAB2_S(ij),          SAD1_S(ij), ...
     SAD2_S(ij),            PARB_S(ij),          PARD_S(ij),       SvF(ij),             SNDtm1(ij),  ...
     snow_albedotm1(ij,:),  Color_Class,         OM_H,             OM_L,                PFT_opt_H,     ...
     PFT_opt_L,             hc_H(ij,:),          hc_L(ij,:),       d_leaf_H,            d_leaf_L,     ...
     Soil_Param,            Interc_Param,        SnowIce_Param,    VegH_Param,          VegL_Param,  ...
     Ca_S(ij),              Oa,                  Citm1_sunH(ij,:), Citm1_shdH(ij,:),    Citm1_sunL(ij,:), ...
     Citm1_shdL(ij,:),      e_rel_H(ij,:),       e_relN_H(ij,:),   e_rel_L(ij,:),       e_relN_L(ij,:), ...
     e_snotm1(ij),          In_Htm1(ij,:),       In_Ltm1(ij,:),    In_Littertm1(ij),    In_urbtm1(ij), ...
     In_rocktm1(ij),        SWEtm1(ij),          In_SWEtm1(ij),    Tdebtm1(ij,:),       Ticetm1(ij), ...
     Tdptm1(ij,:),          Tdp_snowtm1(ij,:),   Tdamptm1(ij),     Ts_undertm1(ij),     WATtm1(ij), ...
     ICEtm1(ij),            IP_wctm1(ij),        ICE_Dtm1(ij),     Cicewtm1(ij),        Vx_Htm1(ij,:), ...
     Vl_Htm1(ij,:),         Vx_Ltm1(ij,:),       Vl_Ltm1(ij,:),    Psi_x_Htm1(ij,:),    Psi_l_Htm1(ij,:), ...
     Psi_x_Ltm1(ij,:),      Psi_l_Ltm1(ij,:),    ZR95_H,           ZR95_L,              FROCKtm1(ij),  ...
     Krock,                 Urb_Par,             Deb_Par,          Zs_deb,              Tdew_S(ij), ...
     t_slstm1(ij),          rostm1(ij),          SP_wctm1(ij),     fpr,                 IrD_S(ij), ...
     In_max_urb,            In_max_rock,         K_usle,           tau_snotm1(ij),      Ta_day(ij,:), ...
     Slo_top2(ij),          Slo_head(ij,:),      Asur(ij),         Ared(ij),            aTop(ij), ...
     EKtm1(ij),             q_runon(ij),         Qi_in(ij,:),      Ws_undertm1(ij),     Pr_sno_t(ij,:), ...
     pow_dis,               a_dis,               Salt_S(ij),       SPAR,                SNn(ij), ...
     OPT_min_SPD,           OPT_VegSnow,         OPT_SoilTemp,     OPT_PlantHydr,       Opt_CR, ...
     Opt_ST,                Opt_ST2,             OPT_SM,           OPT_STh,             OPT_FR_SOIL, ...
     OPT_PH,                parameterize_phase,  hSTL,             Albsno_method);
    %catch ME
    %disp(['Error in HYDROLOGY_MODULE_PAR in ij:' char(num2str(ij))])
    
    % Get the full error identifier (type/ID)
    %errorID = ME.identifier;
    %disp(['Error Type (ID): ', errorID]);

    %end
    
        end
   end %% END OF PARFOR
    
%% LATER STEPS AFTER PARFOR
%==========================================================================
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
    
    [q_runon,    Q_channel,    Qi_in_Rout,   Slo_pot,     Q_exit, ...
     Qsub_exit,  T_pot,        QpointH,      QpointC,     UpointH, ...
     UpointC]= ...
                     ROUTING_MODULE( ...
     dt,         dth,          Rd,           Rh,          Qi_out_Rout, ...
     Q_channel,  cellsize,     Area,         DTM,         NMAN_H, ...
     NMAN_C,     MRough,       WC,           SN,          T_flow, ...
     T_pot,      Slo_top,      ms_max,       POT,         ZWT, ...
     OPT_HEAD,   Xout,         Yout); 
    
 %Q_channel(14, 26)
 %Qi_out_Rout(14,26,8)
 %Rd(14,26)
 %Rh(14,26)

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
%==========================================================================
    Cicewtm1=          Cicew;
    Citm1_shdH=        Ci_shdH;
    Citm1_shdL=        Ci_shdL;
    Citm1_sunH=        Ci_sunH;
    Citm1_sunL=        Ci_sunL;
    e_snotm1=          e_sno;
    EKtm1=             EK;
    FROCKtm1=          FROCK;
    ICE_Dtm1=          ICE_D;
    ICEtm1=            ICE;
    In_Htm1=           In_H ;
    In_Littertm1=      In_Litter;
    In_Ltm1=           In_L;
    In_rocktm1=        In_rock;
    In_SWEtm1=         In_SWE;
    In_urbtm1=         In_urb;
    IP_wctm1=          IP_wc;
    Oicetm1=           Oice;
    Psi_l_Htm1=	       Psi_l_H;
    Psi_l_Ltm1=	       Psi_l_L;
    Psi_x_Htm1=	       Psi_x_H;
    Psi_x_Ltm1=	       Psi_x_L;
    rostm1=            ros;
    SNDtm1=            SND ;
    snow_albedotm1=    snow_albedo ;
    soil_albedotm1=    soil_albedo ;
    SP_wctm1=          SP_wc;
    surface_albedotm1= surface_albedo ;
    SWEtm1=            SWE;
    SWE_avalanchedtm1= SWE_avalanched ;
    t_slstm1=          t_sls;
    tau_snotm1=        tau_sno;
    Tdamptm1=          Tdamp;
    Tdebtm1=           Tdeb;
    Tdp_snowtm1 =      Tdp_snow;
    Tdptm1=            Tdp;
    Ticetm1=           Tice;
    %%% order!
    Tstm1=             Ts;
    Tstm0=             2*Ts-Tstm1;
    Ts_undertm1=       Ts_under;
    %%%
    Vl_Htm1=           Vl_H;
    Vl_Ltm1=           Vl_L;
    Vtm1=              V;
    Vx_Htm1=           Vx_H;
    Vx_Ltm1=           Vx_L;
    WATtm1=            WAT;
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
    
    
    % MEMORY CONDITION FOR VEGETATION  MODEL
    %----------------------------------------------------------------------
    % 1 day
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
       
    % Ccrown_OUT and EVcode are used in the OUTPUT_MANAGER_DIST_LABEL.m
    %--------------------------------------------------------------------------
    % This only for the OUTPUT_MANAGER_DIST_LABEL
    %Ccrown_OUT =[ 1 ; 1; 1 ; 0.9 ; 1; 1; 0];  %% Ccrown fraction for Plant Functional Types (PFT)
    %EVcode = [ 1 2 3 4 5 6 7 8 9 10];             %% code of each PFTs
    %Ccrown_OUT2 = POI.Ccrowns;
    %{
    %Mike
    Ccrown_OUT = [ 1   0 ; 0.5  0.5 ; 1   0 ; 1   0];  %% Ccrown fraction for PFT
    EVcode     = [ 1      2           3       4    ];  %% Code of each PFT
    cc_max = size(Ccrown_OUT,2);                       %% Number of vegetation types per cell (max.)
    %}

    % Generating the vector for the calculation of the outputs by crowns
    %----------------------------------------------------------------------
    max_len = max(cellfun(@numel, POI.Ccrowns));
    padded_cells = cellfun(@(v) [v, zeros(1, max_len - numel(v))], POI.Ccrowns, 'UniformOutput', false);
    Ccrown_OUT = cell2mat(padded_cells);
    
    Ccrown_OUTv = Ccrown_OUT(ksv, :);
    %EVcode     = POI.Class;
  
    % Assignation to temporal matrix
    %----------------------------------------------------------------------
    
    t_date = datetime(Datam_S(1), Datam_S(2), Datam_S(3), Datam_S(4), 0, 0);
    t_assign =  find(t_date  ==  date_forMonth);

    % Outputs for vegetation
    ANPP_H_spatial(:,:,t_assign)   = reshape(sum(ANPP_H.*Ccrown_OUTv,2),xdim,ydim); 
    ANPP_L_spatial(:,:,t_assign)   = reshape(sum(ANPP_L.*Ccrown_OUTv,2),xdim,ydim);
    LAI_H_spatial(:,:,t_assign)    = reshape(sum(LAI_H.*Ccrown_OUTv,2),xdim,ydim);
    LAI_L_spatial(:,:,t_assign)    = reshape(sum(LAI_L.*Ccrown_OUTv,2),xdim,ydim);
    NDVI_spatial(:,:,t_assign)     = reshape(NDVI,xdim,ydim);

    % Output for ET
    EG_spatial(:,:,t_assign)       = reshape(EG,xdim,ydim);
    EICE_spatial(:,:,t_assign)     = reshape(EICE,xdim,ydim);    
    EIn_H_spatial(:,:,t_assign)    = reshape(sum(EIn_H.*Ccrown_OUTv,2),xdim,ydim);
    EIn_L_spatial(:,:,t_assign)    = reshape(sum(EIn_L.*Ccrown_OUTv,2),xdim,ydim);
    EIn_rock_spatial(:,:,t_assign) = reshape(EIn_rock,xdim,ydim); 
    EIn_urb_spatial(:,:,t_assign)  = reshape(EIn_urb,xdim,ydim); 
    EWAT_spatial(:,:,t_assign)     = reshape(EWAT,xdim,ydim); 
    T_H_spatial(:,:,t_assign)      = reshape(sum(T_H.*Ccrown_OUTv,2),xdim,ydim); 
    T_L_spatial(:,:,t_assign)      = reshape(sum(T_L.*Ccrown_OUTv,2),xdim,ydim); 
    SSN_spatial(:,:,t_assign)      = reshape(SSN,xdim,ydim);  
    ESN_spatial(:,:,t_assign)      = reshape(ESN,xdim,ydim);
    ELitter_spatial(:,:,t_assign)  = reshape(ELitter,xdim,ydim);
    ESN_In_spatial(:,:,t_assign)   = reshape(ESN_In,xdim,ydim);

    % Outputs for snow
    SND_spatial(:,:,t_assign)      = reshape(SND,xdim,ydim);
    SWE_spatial(:,:,t_assign)      = reshape(SWE,xdim,ydim);    
        
    % Outputs for groundwater
    FROCK_spatial(:,:,t_assign)    = reshape(FROCK,xdim,ydim); 
    f_spatial(:,:,t_assign)        = reshape(f,xdim,ydim); 
    Gfin_spatial(:,:,t_assign)     = reshape(Gfin,xdim,ydim); 
    %GPP_H_spatial(:,:,t_assign) 
    %GPP_L_spatial(:,:,t_assign) 
    G_spatial(:,:,t_assign)        = reshape(G,xdim,ydim);
    hc_H_spatial(:,:,t_assign)     = reshape(sum(hc_H.*Ccrown_OUTv,2),xdim,ydim);
    hc_L_spatial(:,:,t_assign)     = reshape(sum(hc_H.*Ccrown_OUTv,2),xdim,ydim);

    % Outputs runoff
    %Q_channel_spatial(:,:,t_assign)       = reshape(Q_channel,xdim,ydim);
    %q_runon_spatial(:,:,t_assign)         = reshape(q_runon,xdim,ydim);
    QpointC_series(t_assign,:)     = QpointC;
    
    q_runon_pp   = reshape(q_runon,xdim,ydim);
    Q_channel_pp = reshape(Q_channel,xdim,ydim);

    %disp(q_runon_pp(Yout,Xout))
    %disp(Q_channel_pp(Yout,Xout))
    %disp(QpointC)

    Rd_spatial(:,:,t_assign)       = reshape(Rd,xdim,ydim);

    %Rd_spatial(Yout,Xout, 50);
    
    %% SAVING
    %======================================================================
    % Only save when it is the last day of the month or the last iteration
    %======================================================================
    if t_date  == t_store | t == N_time_step

    % Saving
    %----------------------------------------------------------------------
    disp('Storing results in .mat file')

    save([Directories.save yy '/3D_Outputs_' yy '_' mth '.mat'], ...
         'ANPP_H_spatial',   'ANPP_L_spatial',    'EG_spatial',      'EIn_H_spatial',  'EIn_L_spatial', ...
         'EIn_rock_spatial', 'EIn_urb_spatial',   'EWAT_spatial',    'T_H_spatial',    'T_L_spatial', ...
         'SSN_spatial',      'ESN_spatial',       'ELitter_spatial', 'ESN_In_spatial', 'FROCK_spatial', ...
         'f_spatial',        'Gfin_spatial',      'G_spatial',       'hc_H_spatial',   'hc_L_spatial', ...
         'EICE_spatial',     'LAI_H_spatial',     'LAI_L_spatial', ...
         'NDVI_spatial',     'SND_spatial' ,      'SWE_spatial', ...
         'Rd_spatial',       'QpointC_series' ...
          );

    % Changing the label to create again the matrices for storing
    %--------------------------------------------------------------------------
    output_creation = 0;
    
    %% Memory 
    %----------------------------------------------------------------------
    % Information about the memory used by variables when saving occurs
    %----------------------------------------------------------------------
    
    % Get a structure array of all variables in the current workspace
    %--------------------------------------------------------------------------
    vars = whos;
    
    % Iterate through the variables and sum their sizes
    %--------------------------------------------------------------------------
    totalSize = 0;
    for i = 1:length(vars)
        totalSize = totalSize + vars(i).bytes;
    end
    
    % Display variables size
    %--------------------------------------------------------------------------
    disp(['Memory used by variables at the end of MAIN_FRAME: ', num2str(totalSize/1e6), ' MB']);
    
    % Saving actual conditions in the model to restart if needed
    %----------------------------------------------------------------------
    % run(['INIT_COND_MID.m'])
       
    %save([Directories.save 'Store/Initial_Conditions_' ...
    %  char(num2str(day(t_store))) '_' char(num2str(month(t_store))) '_' char(num2str(year(t_store))) '.mat'])

    INIT_COND_MID(t_store, Directories, ...
     alp_soil,      b_soil,          Bam,              Bem,            BLit, ...
     Ccrown_t,      Cice,            Cicew,            CK1,            Csno,...
     Csnow,         dQ_S,            DQ_S,             dQVEG,          DT_S,...
     dw_SNO,        e_sno,           EG,               EICE,           EIn_rock,...
     EIn_urb,       EK,              ELitter,          er,             ESN_In,...
     SSN_In,        ESN,             SSN,              EWAT,           FROCK,...
     f,             Gfin,            G,                H,              HV,...
     ICE_D,         ICE,             Imelt,            In_H,           In_Litter,...
     In_L,          In_rock,         In_SWE,           In_urb,         IP_wc,...
     Lk_rock,       Lk_wat,          Lk,               Lpho,           NavlI,...
     NDVI,          NIce,            NIn_SWE,          OF,             Oice,...
     OS,            O,               POT,              Pr_liq,         Pr_sno,...
     Q_channel,     Q_exit,          q_runon,          QE,             QEV,...
     Qfm,           Qi_in,                             Qi_out_Rout,    Qi_out,...
     Qsub_exit,     Qv,              r_litter,         r_soil,         ra,...
     Rd,            Rh,              Rn,               ros,            SE_rock,...
     SE_urb,        Slo_head,        Smelt,            SND,            snow_albedo, ...
     soil_albedo,   SP_wc,           surface_albedo,   SWE,            SWE_avalanched, ...
     t_sls,         tau_sno,         Tdamp,            Tdeb,           Tdp, ...
     Tdp_snow,      Tice,            Tstm0,            Ts,             Ts_under, ...
     TsVEG,         U_SWE,           Vice,             V,              WAT, ...
     WIS,           WR_IP,           WR_SP,            Ws_under,       ZWT, ...
     AgeDL_H,       AgeDL_L,         AgeL_H,           AgeL_L,         AgePl_H,...
     AgePl_L,       An_H,            An_L,             ANPP_H,         ANPP_L,...
     B_H,           B_L,             BA_H,             BA_L,           Bfac_dayH,...
     Bfac_dayL,     Bfac_weekH,      Bfac_weekL,       Ci_shdH,        Ci_shdL,...
     Ci_sunH,       Ci_sunL,         dflo_H,           dflo_L,         Dr_H,...
     Dr_L,          e_rel_H,         e_rel_L,          e_relN_H,       e_relN_L,...
     EIn_H,         EIn_L,           fapar_H,          fapar_L,        FNC_H,...
     FNC_L,         gsr_H,           gsr_L,            hc_H,           hc_L,...
     ISOIL_H,       ISOIL_L,         Jsx_H,            Jsx_L,          Jxl_H, ...
     Jxl_L,         Kleaf_H,         Kleaf_L,          Kreserve_H,     Kreserve_L, ...
     Kuptake_H,     Kuptake_L,       Kx_H,             Kx_L,           LAI_H,  ...
     LAI_L,         LAIdead_H,       LAIdead_L,        ManIH,          ManIL, ...
     NBLeaf_H,      NBLeaf_L,        NBLI_H,           NBLI_L,         NPP_H, ...
     NPP_L,         NPPI_H,          NPPI_L,           Nreserve_H,     Nreserve_L, ...
     NuLit_H,       NuLit_L,         NupI_H,           NupI_L,         Nuptake_H, ...
     Nuptake_L,     OH,              OL,               PARI_H,         PARI_L, ...
     PHE_S_H,       PHE_S_L,         Preserve_H,       Preserve_L,     Psi_l_H, ...
     Psi_l_L,       Psi_s_H,         Psi_s_L,          Psi_x_H,        Psi_x_L, ...
     Puptake_H,     Puptake_L,       RA_H,             RA_L,           rap_H, ...
     rap_L,         RB_H,            rb_H,             RB_L,           rb_L, ...
     Rdark_H,       Rdark_L,         Rexmy_H,          Rexmy_L,        Rg_H, ...
     Rg_L,          rKc_H,           rKc_L,            Rmc_H,          Rmc_L, ...
     Rmr_H,         Rmr_L,           Rms_H,            Rms_L,          rNc_H, ...
     rNc_L,         rPc_H,           rPc_L,            Rrootl_H,       Rrootl_L, ...
     rs_shdH,       rs_shdL,         rs_sunH,          rs_sunL,        SAI_H, ...
     SAI_L,         Sfr_H,           Sfr_L,            SIF_H,          SIF_L, ...
     Slf_H,         Slf_L,           Sll_H,            Sll_L,          Sr_H, ...
     Sr_L,          SupK_H,          SupK_L,           SupN_H,         SupN_L, ...
     SupP_H,        SupP_L,          Swm_H,            Swm_L,          T_H, ...
     T_L,           TBio_H,          TBio_L,           Tden_H,         Tden_L, ...
     Tdp_H,         Tdp_L,           TdpI_H,           TdpI_L,         TexC_H, ...
     TexC_L,        TexK_H,          TexK_L,           TexN_H,         TexN_L, ...
     TexP_H,        TexP_L,          TNIT_H,           TNIT_L,         TPHO_H, ...
     TPHO_L,        TPOT_H,          TPOT_L,           Vl_H,           Vl_L, ...
     Vx_H,          Vx_L,            An_H_t,           An_L_t,         O_t, ...
     PAR_t,         Pr_sno_t,        Psi_l_H_t,        Psi_l_L_t,      Psi_x_H_t, ...
     Psi_x_L_t,     Rdark_H_t,       Rdark_L_t,        Ta_t,           Tdp_H_t, ...
     Tdp_L_t,       Tdp_t,           V_t,              Ared,           ICEym1, ...
     SNOWALB,       SWEym1 ...         
        );


    end

    %run(['OUTPUT_MANAGER_DIST_LABEL.m']);

    % Save the workspace at frequent interval. Very useful in case it crashes 
    %if  mod(t,25)==0    
    %    save([outlocation, Fstep], '-regexp', '^(?!(FF_BC|LWIN_BC|PARB|PARD|PP_BC|PRESS_BC|RH_BC|SAB1|SAB2|SAD1|SAD2|TA_BC|WS)$).');
    %end   

    %if  mod(t,8760)==0  ||  t==N_time_step
    %    Fstep2= strcat(Fstep,'_',num2str(t));
    %    save([outlocation, Fstep2], '-regexp', '^(?!(FF_BC|LWIN_BC|PARB|PARD|PP_BC|PRESS_BC|RH_BC|SAB1|SAB2|SAD1|SAD2|TA_BC|WS)$).');
    %end


end
%close(bau)

%% Display
Computational_Time =toc;
disp('COMPUTATIONAL TIME [h] ')
disp(Computational_Time/3600)

%Q_channel