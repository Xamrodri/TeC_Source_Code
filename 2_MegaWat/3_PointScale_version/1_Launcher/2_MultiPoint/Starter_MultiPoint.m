%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%         MULTIPOINT SCALE MODEL                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
%{
Created on April 16, 2024
Author: MAXIMILIANO RODRIGUEZ
Region: Apennines

Code explanation: 
    This code launches the Point Scale version of TC model.
    
%}
%==========================================================================


%% Clear all
clc; clear all;


%% Names of cacthment
SITE = 'Velino';

%Name of the folder to save results and not overwrite previous runs
run_folder = 'Run_8';

%% DIRECTORY
%==========================================================================
%
%==========================================================================
%% Main roots
%root = '/nfs/scistore18/pelligrp/mrodrigu/' %HPC
root = 'C:/Users/mrodrigu/Desktop/19_ISTA/'; %Personal computer
folder_path = [root '1_TC/3_Model_Source/2_MegaWat/']; % Put here the path of where you downloaded the repository

%% Dependencies
addpath(genpath([folder_path,'1_Functions'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path,'5_Common_inputs'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path,'3_Pyrenees_PointScale/2_Forcing'])); % Where is located the meteorological forcing and Shading matrix 
addpath(genpath([folder_path,'3_Pyrenees_PointScale/3_Inputs'])); % Add path to Ca_Data

%% Subfolders
path_model = [root '1_TC/3_Model_Source/2_MegaWat/'];

% Location of the launcher
folder_launcher = [path_model '3_PointScale_version/1_Launcher/Experimental/']; % Put here the path of where you downloaded the repository
addpath(genpath(folder_launcher)); % Add path to Ca_Data

%% Location for outputs - Creation
% Create the directory where model outputs will be stored
outlocation = [folder_path,'3_PointScale_version/4_Outputs/' run_folder '/'];

if ~exist(outlocation, 'dir'); 
mkdir(outlocation); 
    mkdir([outlocation '1_Hourly/']); 
    mkdir([outlocation '2_Daily/']); 
    mkdir([outlocation '3_Params/']); 
    mkdir([outlocation '4_INIT_COND/']); 
    mkdir([outlocation '5_Env/']); 
addpath(genpath(outlocation)); 
end

%% Modelling period
%==========================================================================
% Set
%==========================================================================

date_start = ["01-Jan-2000 00:00:00"]; % Starting point of the simulation
date_end = ["31-Dec-2000 23:00:00"]; % Last timestep of the simulation


%% Main file with points 
%==========================================================================
% In this version, POI keeps all the data needed to run TC.
%==========================================================================
POI = readtable([path_model '3_PointScale_version/3_Inputs/2_Apennine/Velino_MultiPoints.txt']); %import table with points info
UTM_zone = 33; % for Italy
[POI.LAT, POI.LON] = utm2ll(POI.UTM_X, POI.UTM_Y, UTM_zone);

% names of points
names = string(POI.Name);

%% Load DEM and extraction of matrices - Only DTM and VEG_CODE here
%==========================================================================
%{
There are different dtm
        1) dtm_Cinca_Mid_250m.mat
        2) dtm_Tiber_250m.mat
        3) dtm_Velino_250m.mat
DTM is not passed in the function below, but it is reopened inside
run_Point_Pro.m
%}
%==========================================================================

% Load preprocessed data
dtm_file = 'dtm_Velino_250m.mat'; 
mm = load([path_model,'5_Common_inputs/',SITE,'/',dtm_file]); % Distributed maps pre-processing. Useful here to get the DTM and initial snow depth

% DEM and features
DTM = mm.DTM_orig; % Use the full DEM in case running POI outside of mask
DTM(isnan(DTM)) = 0; %%??

[m_cell,n_cell]=size(DTM);
num_cell=numel(DTM);

xllcorner = mm.xllcorner;
yllcorner = mm.yllcorner;
cellsize = mm.cellsize;
VEG_CODE = mm.VEG_CODE; 
%% PIXELS 
% k = 10;
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

%% LAND COVER PARTITION
%==========================================================================
%
%==========================================================================

%% Matrix
ksv=reshape(VEG_CODE,num_cell,1);

%% Surface elevation for vegetation
% categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec   NoVeg]  
zatm_surface = [18      18       2        2        2      18        18      2    ]; %Depend on vegetation
zatm_hourly_on = 0;

%% Loop for vegetation classification
%k=98
for k = 1:height(POI)
%disp(k)
    if ~strcmp(POI.Feature{k}, 'Snow_station') %If the POI is a snow depth station, then it is just bare soil

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
        POI.Ccrown(k) = {[0.2 0.8]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 1; 

        %categories    [fir     larch    grass_A  grass_B  shrub    BLever    BLdec  NoVeg]
        POI.II(k) =    {[0       0       0        0        1        0         1      0    ]>0}; 
        
        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));

 
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
        POI.VEG_CLASS(k) = 2; 
        
        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec   NoVeg]
        POI.II(k) =  {[0       0        1        0        1      0         0       0    ]>0};  

        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 


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
        POI.VEG_CLASS(k) = 3; 

        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg]
        POI.II(k) =  {[0       0        1        0        1      0         0      0    ]>0};
        
        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 
    

        case 4 % Evergreen needleaves at high elevation (>700m) %   
        %    Case 4 includes from CORINE:
        %      1) Coniferous forest (1%)
        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {[0.7 0.3]};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 4; 

        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg]
        POI.II(k) =  {[1       0        0        0        1      0         0      0    ]>0};  
    
        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 


    case 5 % Mediterranean shrublands %
        %    Case 5 includes from CORINE:
        %       1)  Transitional woodland-shrub (5%)
        %       2)  Natural grasslands (4.3%)
        %       3)  Sclerophyllous vegetation (0.4%)
        %       4)  Moors and heathland (>0.05%)
        
        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.2;
        POI.Ccrown(k) = {[0.8]};
   
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 5; 
        
        %% DEBUGGER
        %THERE IS A PROBLEM WITH SHRUBS

        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg]
        POI.II(k) =  {[0       0        0        1        0      0         0      0    ]>0};  
          
        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 
    

    case 6 % Olives %
        %    Case 6 includes from CORINE:
        %       1) Olive groves (4%)

        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {1.0};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 6; 

        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg]
        POI.II(k) =  {[0       0        0        0        1      0         0      0    ]>0};  
    
        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 

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

        POI.Cwat(k) = 0.0; POI.Curb(k) = 1.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {0.0};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 7; 
        
        %%CHECK  WHY IS NEEDED TO SET THE VEGETATION PARAMTERS
        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg]
        POI.II(k) =  {[1       0        0        0        0      0         0      0    ]>0};  
        
        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 

    
    %% CHECK THE NEED OF SETTING II AS 1 WHEN CCROCK = 1 
    %======================================================================
    %{
    Problem is with these parameters defined in the MAIN_FRAME_Pro
    % Vegetation Optical Parameter
    [PFT_opt_H(i)]=Veg_Optical_Parameter(OPT_PROP_H(i));
    [PFT_opt_L(i)]=Veg_Optical_Parameter(OPT_PROP_L(i));

    [Stoich_H(i)]=Veg_Stoichiometric_Parameter(Nl_H(i));
    [ParEx_H(i)]=Exudation_Parameter(0);
    [Mpar_H(i)]=Vegetation_Management_Parameter;

    [Stoich_L(i)]=Veg_Stoichiometric_Parameter(Nl_L(i));
    [ParEx_L(i)]=Exudation_Parameter(0);
    [Mpar_L(i)]=Vegetation_Management_Parameter;

    %}
    %======================================================================


    case 8 % Rock %
        %    Case 8 includes from CORINE:
        %        1) Bare rocks (0.2%)
        %        2) Glaciers and perpetual snow (0%)

        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 1.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {0.0};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 8; 

        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg ]  
        POI.II(k) =  {[0       0        1        0        0      0         0      0     ]>0}; 

        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 

    case 9 % Water %
        %    Case 9 includes from CORINE:
        %        1) Water bodies (0.3%)
        %        2) Water courses (0.2%)

        %     It also includes, but not in Tiber basin:
        %       1) Coastal lagoons
        %       2) Estuaries
        %       3) Sea and Ocean


        POI.Cwat(k) = 1.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        POI.Ccrown(k) = {0.0};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 9; 

        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg]  
        POI.II(k) =  {[0       0        1        0        0      0         0      0    ]>0};
  
        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 

    case 10 % Bare soils %
        %    Case 10 includes from CORINE:
        %        1) Mineral extraction sites (0.2%)
        %        2) Burnt areas (0.1%)
        %        3) Beaches - dunes - sands (>0.05%)
        %        4) Sparsely vegetated areas (0.7%)
        
        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.9;
        POI.Ccrown(k) = {0.1};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 10; 

        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg]  
        POI.II(k) =  {[0       0        1        0        0      0         0      0    ]>0};  

        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));

        otherwise
        disp('INDEX FOR VEGETATION PARAMETER INCONSISTENT')
        %return
        end
        
    else % if I have an station, this is in bare soil.    
      
        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 1.0;
        POI.Ccrown(k) = {0.0};
        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.VEG_CLASS(k) = 10; 

        %categories   [fir     larch    grass_A  grass_B  shrub  BLever    BLdec  NoVeg]  
        POI.II(k) =  {[0       0        1        0        0      0         0      0    ]>0};  

        % cc_max
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));

    end


% Defining zatm in each point
if ~isempty(zatm_surface(cell2mat(POI.II(k))>0))
POI.zatm(k) = max(zatm_surface(cell2mat(POI.II(k))>0)); %choose correct atmospheric reference height
else
POI.zatm(k) = 2;
end

end

%% Memory in
%profile on -memory

%% Single point Launcher
run_Point_Pro(root, outlocation, run_folder, names(78), POI, ksv, date_start, date_end);

%% Memory out
%profile off
%profile viewer

%{

%% Workers - Not needed for the cluster
% Specify the desired number of workers

numWorkers = 4; % Example: Set to 4 workers

% Create a parallel pool with the specified number of workers
poolobj = gcp('nocreate'); % Check if a pool already exists
if isempty(poolobj)
    parpool(numWorkers); % Create a new pool if one doesn't exist
elseif poolobj.NumWorkers ~= numWorkers
    % If a pool exists with a different size, delete it and create a new one
    delete(poolobj);
    parpool(numWorkers);
end



%% Computational time
%Time counter
tic;

%% Specific tasks - Here define the runs
% Define the specific indices you want to run

specific_indices = 2:length(names); % for all the runs
%specific_indices = [1, 23];
%specific_indices = [5, 11, 14, 22, 28, 30, 36, 38, 41, 55, 59, 65, 68, 76, 77, 82, 84, 90, 92, 95];

% Create a new cell array of names based on the specific indices
selected_names = names(specific_indices);

%% Parallel computing launch for the Model
parfor k = 1:length(selected_names) %length(names)  % 3    
    try
        run_Point_Pro(root, outlocation, run_folder, selected_names(k), POI, ksv, date_start, date_end)
    catch ME
        warning('Error occurred on worker %d: %s', k, ME.message);       
    end
end

%% Computational time
Computational_Time =toc;
disp('COMPUTATIONAL TIME PARFOR [h] ')
disp(Computational_Time/3600)


%}