%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%         MULTIPOINT SCALE MODEL                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
%{
Created on April 16, 2025
Author: MAXIMILIANO RODRIGUEZ
Region: Apennines

Code explanation: 
    This code launches the Point Scale version of TC model.
    
%}
%==========================================================================

%% CLEAR ALL
clc; clear;

%% INITIAL CONDITIONS
%==========================================================================

% Names of catchment and clusters
%--------------------------------------------------------------------------
IniCond.SITE = 'Tiber';
IniCond.Clusters = '501points';
IniCond.FORCING = "ERA5Land";
IniCond.DeltaGMT= 1; % for Italy
IniCond.versionDEM= '202595_250m';
% Name of the folder to save results
%--------------------------------------------------------------------------
IniCond.run_folder = 'Run_36';

% Modelling period
%--------------------------------------------------------------------------
dateRun.start = "01-Oct-1999 00:00:00"; % Starting point of the simulation
dateRun.end = "30-Sep-2000 23:00:00"; % Last timestep of the simulation

% Folder with forcings
%--------------------------------------------------------------------------
%{
20250806_A for 50points in Velino
20250723   for 100points in Velino
20250806_B for 500points in Velino
20250806_C for 999points in Velino
20250818 for 25points in Pianello
20250829 for 500points in Tiber
%}
forc_in = '20250903';

% Folder with Terrain inputs
%--------------------------------------------------------------------------
%{
1_Velino
2_Pianello_in_Chiascio
3_Tiber
%}
terrain_in = '3_Tiber';

% Choose point or multipoint
%--------------------------------------------------------------------------
IniCond.mode = "point";
% if point mode, then select the point
selected_point = "TiberCluster1";

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
Directories.save = [Directories.root '1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/4_Outputs/' IniCond.run_folder '/']; 

% Sub-path for forcings
%--------------------------------------------------------------------------
Directories.forc = [Directories.root '2_Forcing/3_Downscalling_ERA5/5_Radiation_Partition/1_Bias_corrected/' forc_in '/']; 

% Sub-path for inputs - DEM with version
%--------------------------------------------------------------------------
Directories.DEMinputs = [Directories.root '1_TC/3_Model_Source/2_MegaWat/4_Preparation_files/3_Distributed_GeoInputs_Apennines/8_Output/' IniCond.SITE '/' IniCond.versionDEM '/']; 

% Sub-path for terrain
%--------------------------------------------------------------------------
Directories.terrain = [Directories.model,'4_Preparation_files/4_GeoTerrain_MultiPoint/2_Results/' terrain_in '/1_ByCluster/' IniCond.Clusters '/']; 

% Dependencies
%--------------------------------------------------------------------------
addpath(genpath([Directories.model,'1_Functions'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([Directories.model,'5_Common_inputs'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([Directories.model,'3_Pyrenees_PointScale/2_Forcing'])); % Where is located the meteorological forcing and Shading matrix 
addpath(genpath([Directories.model,'3_Pyrenees_PointScale/3_Inputs'])); % Add path to Ca_Data


%% LOCATION OF OUTPUTS - CREATION OF FOLDERS
%==========================================================================

if ~exist(Directories.save, 'dir') 
    disp('Folders do not exist for outputs. Creating new folders')
    mkdir(Directories.save); 
    mkdir([Directories.save '1_Hourly/']); 
    mkdir([Directories.save '2_Daily/']); 
    mkdir([Directories.save '3_Params/']); 
    mkdir([Directories.save '4_INIT_COND/']); 
    mkdir([Directories.save '5_Env/']); 
    mkdir([Directories.save '6_POI_table/']); 
addpath(genpath(Directories.save)); 
end


%% MODEL INPUTS
%==========================================================================

% Points of interest
%--------------------------------------------------------------------------
POI = readtable([Directories.model '3_PointScale_version/3_Inputs/2_Apennine/' IniCond.SITE '_' IniCond.Clusters '.txt']); %import table with points info
UTM_zone = 33; % for Italy
[POI.LAT, POI.LON] = utm2ll(POI.UTM_X, POI.UTM_Y, UTM_zone);

% names of points
names = string(POI.Name);

% Parameters for vegetation
%--------------------------------------------------------------------------
%{
Here are the parameters of the model for vegetation.
opts check the format of the columns. 
Use opts to force all the columns with values to be numeric.
%}
%--------------------------------------------------------------------------

opts = detectImportOptions([Directories.root '1_TC/3_Model_Source/2_MegaWat/5_Common_inputs/Parameters_TC.xlsx']);
opts = setvartype(opts, [7:length(opts.VariableTypes)], 'double');

TT_par = readtable([Directories.root '1_TC/3_Model_Source/2_MegaWat/5_Common_inputs/Parameters_TC.xlsx'], opts);

%% DEM AND EXTRACTION OF MATRICES
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
%--------------------------------------------------------------------------
dtm_file = ['dtm_' IniCond.SITE '_250m.mat'];
MAIN_MATs = load([Directories.DEMinputs dtm_file]); % Distributed maps pre-processing. Useful here to get the DTM and initial snow depth

% Extraction of DEM and features
%--------------------------------------------------------------------------
DTM = MAIN_MATs.DTM_orig; % Use the full DEM in case running POI outside of mask
DTM(isnan(DTM)) = 0; %%??

[m_cell,n_cell]=size(DTM);
num_cell=numel(DTM);

xllcorner = MAIN_MATs.xllcorner; %bottom corner
yllcorner = MAIN_MATs.yllcorner; %bottom corner
cellsize = MAIN_MATs.cellsize;
VEG_CODE = MAIN_MATs.VEG_CODE; 

%figure()
%imagesc(mm.DTM)

%% PIXELS
%==========================================================================
% Construction of POI matrix with details of points.
%==========================================================================
k = 14;
for k = 1:height(POI)
id_location = char(string(POI.Name(k))); %id

y_coord = POI.UTM_Y(k);
x_coord = POI.UTM_X(k);

pos_col = floor((x_coord - xllcorner) / cellsize) + 1;
pos_row = floor((y_coord - yllcorner) / cellsize) + 1;

%ij = POI.idx(loc);
ij = sub2ind(size(DTM),pos_row, pos_col); % Location
[j, i] = ind2sub(size(DTM), ij); % Location

POI.ij(k) = ij;
POI.i(k) = i;
POI.j(k) = j;

POI.Zbas(k) = DTM(j,i); % Altitude
end

%Check
%Mask = mm.MASK;
%DTMp = DTM;
%DTMp(Mask == 0) = 0;
%figure()
%imagesc(DTMp)
%[r, c] = ind2sub(size(DTM), ij);
%hold on; % This keeps the image on the screen while you add the point
%plot(c, r, 'r*', 'MarkerSize', 10, 'LineWidth', 2);

%% LAND COVER PARTITION
%==========================================================================
% Here the model decides what means each land cover category.
%==========================================================================

%% Matrix
ksv=reshape(VEG_CODE,num_cell,1);

%% Surface elevation for vegetation
%--------------------------------------------------------------------------
% Here, note that it is needed to define the starting point of the data in
% the excel file. In this case it is from column 5 to the end. 
%--------------------------------------------------------------------------
zatm_surface = TT_par(strcmp(TT_par.Parameters,'zatm_surface'),7:size(TT_par,2));
zatm_hourly_on = 0;

%% Representation and definition of each vegetation class
%--------------------------------------------------------------------------
POIveg = table('Size', [10, 7], ...
          'VariableTypes', {'double', 'string', 'double', 'double', 'double', 'double', 'double'}, ...
          'VariableNames', {'Class', 'Veg_type', 'Ccrowns','Cwat','Curb','Crock','Cbare'});

POIveg.Class = (1:10)';

POIveg.Veg_type = {["BLdec_lowElev_ro2" "grass_Bondone"];                                      % Class 1:  Decidious Broad-leaved forest
                   ["grass_Bondone"];                                                          % Class 2:  Grassland/pasture
                   ["Crops_WW"  "Crops_WB"  "Crops_S"  "Crops_R" "grass_Bondone"];             % Class 3:  Crops
                   ["EvGreen_NeedLeaves_LeBray" "BLdec_highElev_Collelongo"];                  % Class 4:  Evergreen needleaves at high elevation (>700m)
                   ["shrub_Balsablanca_A" "shrub_Balsablanca_B" "BLdec_highElev_Collelongo"];  % Class 5:  Mediterranean shrublands
                   ["crops_Negrisia" "grass_Bondone"];                                         % Class 6:  Olives
                   ["NoVeg"];                                                                  % Class 7:  Urban
                   ["NoVeg"];                                                                  % Class 8:  Rock
                   ["NoVeg"];                                                                  % Class 9:  Water
                   ["grass_Alinya"]                                                            % Class 10: Bare soils
                  };

POIveg.Ccrowns =  {[0.8 0.2];
                   [0.9];
                   [0.1 0.1 0.2 0.1 0.3];
                   [0.5 0.2];
                   [0.6 0.05 0.2];
                   [0.4 0.2];
                   [0.0];
                   [0.0];
                   [0.0];
                   [0.1]};

% From class 1 to 10
%                 1     2     3     4     5      6     7     8     9     10
POIveg.Cwat  = ({ 0.0, 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,  0.0,  1.0,  0.0})';
POIveg.Curb  = ({ 0.0, 0.0,  0.0,   0.0,  0.0,   0.0,  1.0,  0.0,  0.0,  0.0})';
POIveg.Crock = ({ 0.0, 0.0,  0.0,   0.0,  0.0,   0.0,  0.0,  1.0,  0.0,  0.0})';
POIveg.Cbare = ({ 0.0, 0.1,  0.2,   0.3,  0.15,  0.4,  0.0,  0.0,  0.0,  0.9})';


% Assigning vegetation classification in POI table
%--------------------------------------------------------------------------
k=1
for k = 1:height(POI)
%disp(k)
    if ~strcmp(POI.Feature{k}, 'Snow_station') %If the POI is a snow depth station, then it is just bare soil

    vegCode = ksv(POI.ij(k));

    POI.VEG_CLASS(k) = vegCode;

    POI.Cwat(k) = POIveg.Cwat{vegCode};
    POI.Curb(k) = POIveg.Curb{vegCode};
    POI.Crock(k) = POIveg.Crock{vegCode};
    POI.Cbare(k) = POIveg.Cbare{vegCode};

    POI.II(k) =  {POIveg.Veg_type{vegCode}};   
    POI.Ccrown(k) = {POIveg.Ccrowns{vegCode}};

    POI.NCrown(k) = length(POIveg.Ccrowns{vegCode});
    POI.cc_max(k) = length(POIveg.Ccrowns{vegCode}); 
        
    else % if I have an station, this is in bare soil.    
      
    POI.VEG_CLASS(k) = 10;

    POI.Cwat(k) = 0.0;
    POI.Curb(k) = 0.0;
    POI.Crock(k) = 0.0;
    POI.Cbare(k) = 1.0;
             
    POI.II(k) =  {["grass_Alinya"]};  
    POI.Ccrown(k) = {[0.0]};

    POI.NCrown(k) = 1;      
    POI.cc_max(k) = 1;

    end


    % Defining zatm in each point
    if ~isempty(zatm_surface(:,POI.II{k}))
    POI.zatm(k) = {table2array(zatm_surface(:,POI.II{k}))}; %choose correct atmospheric reference height
    else
    POI.zatm(k) = {2};
    end

end

%% SAVE POI TABLE FOR REFERENCE
%==========================================================================
writetable(POI, [Directories.save '6_POI_table/POI_' IniCond.run_folder '.csv']);


%% SINGLE POINT LAUNCHER FOR PERSONAL COMPUTER
%==========================================================================

if IniCond.mode == "point"

% Time counter
%--------------------------------------------------------------------------
tic;

%Main function
%--------------------------------------------------------------------------
run_Point_Pro_BC(Directories, IniCond, selected_point, POI, ksv, dateRun, TT_par, zatm_surface, MAIN_MATs);

% Computational time
%--------------------------------------------------------------------------
Computational_Time =toc;
disp('Computation time - Single Point')
disp([num2str(round(Computational_Time/60,1)) ' mins'])

%% Memory use - Windows only
if ~contains(Directories.root,"nfs") 
[user, sys] = memory; % Windows only
disp(['Mem used in the Starter: ', num2str(user.MemUsedMATLAB/1e6), ' MB'])
end

% Memory out
%--------------------------------------------------------------------------
%profile off
%profile viewer


%% MULTIPLE POINT LAUNCHER 
%==========================================================================
% Specify the desired number of workers
%==========================================================================

elseif IniCond.mode == "multipoint"

% Workers - Only for personal computer
%--------------------------------------------------------------------------

%if I am working on my personal computer define this
if ~contains(Directories.root,"nfs") 

    numWorkers = 3; % Example: Set to 4 workers
    
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


% Computational time - Time counter - For personal computer or HPC
%--------------------------------------------------------------------------
tic;

% Specific tasks - Here define the runs
%--------------------------------------------------------------------------

idx = contains(names, "Cluster");
names = names(idx);

specific_indices = 1:length(names); % for all the runs
%specific_indices = [1, 23];
%specific_indices = [5, 11, 14, 22, 28, 30, 36, 38, 41, 55, 59, 65, 68, 76, 77, 82, 84, 90, 92, 95];
specific_indices = 1:20;

% Create a new cell array of names based on the specific indices
selected_names = names(specific_indices);

% Parallel computing launch for the Model
%--------------------------------------------------------------------------
parfor k = 1:length(selected_names) %length(names)  % 3    
    try
        run_Point_Pro_BC(Directories, IniCond, selected_names(k), POI, ksv, dateRun, TT_par, zatm_surface, MAIN_MATs);
    catch ME
        warning('Error occurred on worker %d: %s', k, ME.message);       
    end
end

% Computational time - End of processing
%--------------------------------------------------------------------------
Computational_Time =toc;
disp('Computation time - MultiPoint ')
disp([num2str(round(Computational_Time/60,1)) ' mins'])

%% Memory use - Windows only
if ~contains(Directories.root,"nfs") 
[user, sys] = memory; % Windows only
disp(['Mem used by the Starter: ', num2str(user.MemUsedMATLAB/1e6), ' MB'])
end

else %if the label is not point or multipoint, define it well
disp('The label is not right. Choose "point" or "multipoint"')
end