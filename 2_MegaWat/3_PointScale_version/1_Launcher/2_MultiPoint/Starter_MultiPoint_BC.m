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

%% CLEAR ALL
clc; clear;

%% INITIAL CONDITIONS
%==========================================================================

% Names of catchment and clusters
%--------------------------------------------------------------------------
SITE = 'Velino';
Clusters = '500points';

% Name of the folder to save results
%--------------------------------------------------------------------------
run_folder = 'Run_24';

% Modelling period
%--------------------------------------------------------------------------
date_start = ["01-Oct-1999 00:00:00"]; % Starting point of the simulation
date_end = ["30-Sep-2000 23:00:00"]; % Last timestep of the simulation

% Folder with forcings
%--------------------------------------------------------------------------
forc_in = '20250806_B';

%% DIRECTORIES
%==========================================================================
% All paths are here
%==========================================================================

% Main roots
%--------------------------------------------------------------------------
%root = '/nfs/scistore18/pelligrp/mrodrigu/' %HPC
Directories.root = 'C:/Users/mrodrigu/Desktop/19_ISTA/1_Science_MegaWat/'; %Personal computer

% Sub-path for the model
%--------------------------------------------------------------------------
Directories.model = [Directories.root '1_TC/3_Model_Source/2_MegaWat/'];

% Sub-path for outputs
%--------------------------------------------------------------------------
Directories.save = [Directories.root '1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/4_Outputs/' run_folder '/']; 

% Sub-path for forcings
%--------------------------------------------------------------------------
Directories.forc = [Directories.root '2_Forcing/3_Downscalling_ERA5/5_Radiation_Partition/1_Bias_corrected/' forc_in '/']; % Put here the path of where you downloaded the repository;

% Dependencies
%--------------------------------------------------------------------------
addpath(genpath([Directories.model,'1_Functions'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([Directories.model,'5_Common_inputs'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([Directories.model,'3_Pyrenees_PointScale/2_Forcing'])); % Where is located the meteorological forcing and Shading matrix 
addpath(genpath([Directories.model,'3_Pyrenees_PointScale/3_Inputs'])); % Add path to Ca_Data


%% LOCATION OF OUTPUTS - CREATION OF FOLDERS
%==========================================================================

if ~exist(Directories.save, 'dir'); 
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
POI = readtable([Directories.model '3_PointScale_version/3_Inputs/2_Apennine/Velino_' Clusters '.txt']); %import table with points info
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

opts = detectImportOptions([Directories.root '1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/3_Inputs/2_Apennine/Parameters_TC.xlsx']);
opts = setvartype(opts, [7:length(opts.VariableTypes)], 'double');

TT_par = readtable([Directories.root '1_TC/3_Model_Source/2_MegaWat/3_PointScale_version/3_Inputs/2_Apennine/Parameters_TC.xlsx'], opts);

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
dtm_file = 'dtm_Velino_250m.mat';
mm = load([Directories.model,'5_Common_inputs/',SITE,'/',dtm_file]); % Distributed maps pre-processing. Useful here to get the DTM and initial snow depth

% Extraction of DEM and features
%--------------------------------------------------------------------------
DTM = mm.DTM_orig; % Use the full DEM in case running POI outside of mask
DTM(isnan(DTM)) = 0; %%??

[m_cell,n_cell]=size(DTM);
num_cell=numel(DTM);

xllcorner = mm.xllcorner; %bottom corner
yllcorner = mm.yllcorner; %bottom corner
cellsize = mm.cellsize;
VEG_CODE = mm.VEG_CODE; 

%figure()
%imagesc(mm.DTM)

%% PIXELS
%==========================================================================
% Construction of POI matrix with details of points.
%==========================================================================
% k = 10;
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

%% Loop for vegetation classification
%k=3
for k = 1:height(POI)
%disp(k)
    if ~strcmp(POI.Feature{k}, 'Snow_station') %If the POI is a snow depth station, then it is just bare soil

    switch ksv(POI.ij(k))

    case 1 % Decidious Broad-leaved forest %
        % CORINE INFORMATION FOR TIBER BASIN 
        %------------------------------------------------------------------
        %   1) Broad-leaved forest (32.5%)
        %   2) Mixed forest (1.2%)
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
        % EXAMPLES OF SHRUBS
        %------------------------------------------------------------------
        %   1) Pixel 2, 3
        %
        % OPTIONS FOR CALIBRATION
        %------------------------------------------------------------------
        %   1) Collelongo         
        %      Original: 
        %           ["BLdec_highElev_Collelongo"]
        %           Ccrown = [1.0];
        %      Used: 
        %           ["BLdec_highElev_Collelongo" "grass_Bondone"]
        %           Ccrown = [0.8 0.2];         
        %      Comment for pixel 2:  LAI continues strongly into the winter. 
        %      Start ok.
        %
        %   2) Rocca Respampani        
        %      Original: 
        %           ["BLdec_lowElev_ro2"]
        %           Ccrown = [1.0];
        %      Used: 
        %           ["BLdec_lowElev_ro2"]
        %           Ccrown = [1.0];
        %      Comment for pixel 2: Sharp shape. LAI continues strongly into the winter.
        %      Start ok.
        %
        %   3) NGreece
        %      Original: 
        %          ["BLdec_lowElev_NGreece" "BLdec_lowElev_grass_NGreece"]
        %           Ccrown = [0.95 0.05];  
        %      Used:
        %          ["BLdec_lowElev_NGreece" "BLdec_lowElev_grass_NGreece"]
        %          Ccrown = [0.9 0.1];
        %      Comment for pixel 2: Mostly ok, but LAI continues well
        %      advanced the winter. Start ok.
        %% CODE
        POI.VEG_CLASS(k) = 1;

        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;

        POI.II(k) =  {["BLdec_lowElev_ro2" "grass_Bondone"]};   
        POI.Ccrown(k) = {[0.8 0.2]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));

 
    case 2 % Grassland/pasture %
        % CORINE INFORMATION FOR TIBER BASIN 
        %------------------------------------------------------------------
        %    1) Natural grasslands (5.0%)
        %    2) Pastures (1.3%)       
        %    3) Green urban areas (0.1%)
        %    4) Inland marshes (>0.05%)
        %
        %    It also includes, but not in Tiber basin:
        %    1) Peat bogs
        %    2) Salt marshes
        %    3) Salines
        %    4) Intertaidal flats

        % EXAMPLES OF SHRUBS
        %------------------------------------------------------------------
        %     1) Pixel 45

        % OPTIONS FOR CALIBRATION
        %------------------------------------------------------------------
        %   1) Vall d'Alinyà
        %     
        %      Original: 
        %           ["grass_Alinya"]
        %           Ccrown = [0.85];
        %           Cbare = 0.15;
        %      Used: 
        %           ["grass_Alinya"]
        %           Ccrown = [0.85]; 
        %           Cbare = 0.15;
        %      Comment for pixel 45: LAI keeps growing at first. Unknown reason. 
        %
        %   2) Monte Bondone
        %    
        %      Original:
        %           ["grass_Bondone"]
        %           Ccrown = [1.0];
        %      Used: 
        %           ["grass_Bondone"]
        %           Ccrown = [0.9]; 
        %           Cbare = 0.1;
        %      Comment for pixel 45: Ok, but LAI continues well advanced
        %      the winter.
        %
        %   3) Las Majadas del Tietar
        %      Original:
        %           ["grass_Tietar_A" "grass_Tietar_B"]
        %           Ccrown = [0.2 0.8];
        %      Used: 
        %           ["grass_Tietar_A" "grass_Tietar_B"]
        %           Ccrown = [0.2 0.8];
        %      Comment for pixel 45: Unrecognize variable phi1 from module
        %      Canopy_Radiative_Transfer. Parameter seems not to get
        %      defined when the set of parameters is used.
        %
        %% CODE
        POI.VEG_CLASS(k) = 2; 

        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.1;
                
        POI.II(k) =  {["grass_Bondone"]};
        POI.Ccrown(k) = {[0.9]};  

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 


    case 3 % Crops %
        % CORINE INFORMATION FOR TIBER BASIN 
        %------------------------------------------------------------------
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
        % EXAMPLES OF SHRUBS
        %------------------------------------------------------------------
        %     1) Pixel 55
        %
        % OPTIONS FOR CALIBRATION
        %------------------------------------------------------------------
        %   1) Aurade
        %      Original:
        %         ["Crops_WW" "Crops_WB" "Crops_S" "Crops_R"]
        %         Original: Ccrown = [0.0 0.0 0.0 0.0]; (Probably not well set)
        %         Cbare = 1.0;
        %      Used:
        %         ["Crops_WW" "Crops_WB" "Crops_S" "Crops_R" "grass_Bondone"]
        %         Used: Ccrown = [0.1 0.1 0.2 0.1 0.3]; 
        %         Cbare = 0.2; 
        %      Comment for pixel 55: 
        %
        %% CODE
        POI.VEG_CLASS(k) = 3;
        
        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.2;
                
        POI.II(k) =  {["Crops_WW"  "Crops_WB"  "Crops_S"  "Crops_R" "grass_Bondone"]};
        POI.Ccrown(k) = {[0.1 0.1 0.2 0.1 0.3]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 
    

    case 4 % Evergreen needleaves at high elevation (>700m) %   
        % CORINE INFORMATION FOR TIBER BASIN 
        %------------------------------------------------------------------
        %     1) Coniferous forest (1%)
        %
        % EXAMPLES OF EVERGREEN NEEDLEAVES
        %------------------------------------------------------------------
        %     1) Pixel 72
        %
        % OPTIONS FOR CALIBRATION
        %------------------------------------------------------------------
        %     1) Renon
        %        Original: 
        %           ["EvGreen_NeedLeaves_Renon"]
        %           Ccrown = [1.0];
        %        Used:
        %           ["EvGreen_NeedLeaves_Renon" "BLdec_highElev_Collelongo"] 
        %           Ccrown = [0.5 0.2];
        %           Cbare = 0.3;
        %        Comment for pixel 72: MODIS data seems to not behave like
        %        a evergreen forest. Check!
        %
        %     2) Le Bray
        %        Original: 
        %           ["EvGreen_NeedLeaves_LeBray"]
        %           Ccrown = [1.0];
        %        Used:
        %           ["EvGreen_NeedLeaves_LeBray" "BLdec_highElev_Collelongo"]
        %           Ccrown = [0.5 0.2];
        %           Cbare = 0.3;
        %        Comment for pixel 72: MODIS data seems to not behave like
        %        a evergreen forest. Check!
        %% CODE
        POI.VEG_CLASS(k) = 4;
        
        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.3;
                 
        POI.II(k) =  {["EvGreen_NeedLeaves_LeBray" "BLdec_highElev_Collelongo"]};  
        POI.Ccrown(k) = {[0.5 0.2]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));
    
    case 5 % Mediterranean shrublands %        
    % CORINE INFORMATION FOR TIBER BASIN  
    %----------------------------------------------------------------------
    %   1)  Transitional woodland-shrub (5%)    
    %   2)  Sclerophyllous vegetation (0.4%)
    %   3)  Moors and heathland (>0.05%)
    %
    % EXAMPLES OF SHRUBS
    %----------------------------------------------------------------------
    %   1) Pixel 89
    %
    % OPTIONS FOR CALIBRATION 
    %----------------------------------------------------------------------
    %   1) Garraf 
    %      Original:
    %           ["shrub_Garraf_A" "shrub_Garraf_B"]
    %           Ccrown = [0.01 0.70]; 
    %           Cbare = 0.29;
    %      Comment for pixel 89: Vegetation do not grow
    %
    %   2) Noe 
    %      Original: 
    %           ["shrub_Noe_A" "shrub_Noe_B"]
    %           Ccrown = [0.56 0.24];
    %           Cbare = 0.20;
    %      Comment for pixel 89: Not good. Constant increase
    %
    %   3) Balsablanca 
    %      (shrub_Balsablanca_A, shrub_Balsablanca_B)
    %      Original:
    %           ["shrub_Balsablanca_A" "shrub_Balsablanca_B"]
    %           Ccrown = [0.58 0.05]; 
    %           Cbare = 0.37; 
    %      Used:  
    %           ["shrub_Balsablanca_A" "shrub_Balsablanca_B" "BLdec_highElev_Collelongo"]
    %           Ccrown = [0.6 0.05 0.2]; 
    %           Cbare = 0.15;
    %      Comment for pixel 89: Relatevely ok. Peak displaced to the summer.
    %
    %   4) Achille
    %      Original: 
    %           ["shrub_Achille"]
    %           Ccrown = [0.8]; 
    %           Cbare = 0.2; 
    %      Comment for pixel 89: For mountains.  Relatevely ok. Peak displaced to the summer.
    %
    %% CODE
        POI.VEG_CLASS(k) = 5;

        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.15; %0.15
                
        POI.II(k) =  {["shrub_Balsablanca_A" "shrub_Balsablanca_B" "BLdec_highElev_Collelongo"]}; %"shrub_Balsablanca_A" "shrub_Balsablanca_B" "BLdec_highElev_Collelongo"
        POI.Ccrown(k) = {[0.6 0.05 0.2]};  %[0.6 0.05 0.2]

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k))); 
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));

    case 6 % Olives %
    % CORINE INFORMATION FOR TIBER BASIN   
    %----------------------------------------------------------------------
    %    1) Olive groves (4%)
    % EXAMPLES OF OLIVES
    %----------------------------------------------------------------------
    %    1) Pixel 94
    %
    % OPTIONS FOR CALIBRATION 
    %----------------------------------------------------------------------
    %    1) Negrisia 
    %       Original: 
    %           ["crops_Negrisia"]
    %           Ccrown = [0.75];
    %           Cbare = 0.25;
    %       Used:          
    %           ["crops_Negrisia" "grass_Bondone"]
    %           Ccrown = [0.6 0.2]; 
    %           Cbare = 0.25;
    %       Comment for pixel 94: Vegetation do not grow
    %
    %    2) Aurade
    %     (Crops_WW, Crops_WB, Crops_S, Crops_R)
    %      Original: 
    %           ["Crops_WW", "Crops_WB", "Crops_S", "Crops_R"]
    %           Ccrown = [0.0 0.0 0.0 0.0]; (Probably not well set)
    %           Cbare =  1.0; 
    %      Used: 
    %           Ccrown = [];
    %      Comment for pixel 94: 
    %% CODE
        POI.VEG_CLASS(k) = 6;

        POI.Cwat(k) = 0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.4;
                
        POI.II(k) =  {["crops_Negrisia" "grass_Bondone"]};
        POI.Ccrown(k) = {[0.4 0.2]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));                  
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 

    case 7 % Urban %
    % CORINE INFORMATION FOR TIBER BASIN    
    %----------------------------------------------------------------------
    %         1) Discontinuous urban fabric (3%)
    %         2) Industrial or commercial units (0.7%)
    %         3) Continuous urban fabric (0.6%)
    %         4) Sport and leisure facilities (0.1%)
    %         5) Road and rail networks and associated land (0.1%)
    %         6) Construction sites (0.1%)
    %         7) Port areas (>0.05%)   
    %         8) Airports (0.1%)
    %% CODE

        POI.VEG_CLASS(k) = 7;

        POI.Cwat(k) = 0.0; POI.Curb(k) = 1.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
        
        POI.II(k) =  {["NoVeg"]};
        POI.Ccrown(k) = {[0.0]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));
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
    % CORINE INFORMATION FOR TIBER BASIN   
    %----------------------------------------------------------------------
    %        1) Bare rocks (0.2%)
    %        2) Glaciers and perpetual snow (0%)
    %% CODE
        POI.VEG_CLASS(k) = 8;

        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 1.0; POI.Cbare(k) = 0.0;
        
        POI.II(k) =  {["NoVeg"]}; 
        POI.Ccrown(k) = {[0.0]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));         
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));

    case 9 % Water %
    % CORINE INFORMATION FOR TIBER BASIN
    %----------------------------------------------------------------------
    %        1) Water bodies (0.3%)
    %        2) Water courses (0.2%)

    %     It also includes, but not in Tiber basin:
    %       1) Coastal lagoons
    %       2) Estuaries
    %       3) Sea and Ocean
    %% CODE

        POI.VEG_CLASS(k) = 9;

        POI.Cwat(k) = 1.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.0;
                 
        POI.II(k) =  {["NoVeg"]};
        POI.Ccrown(k) = {[0.0]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));  
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k))); 

    case 10 % Bare soils %
    % CORINE INFORMATION FOR TIBER BASIN   
    %----------------------------------------------------------------------
    %        1) Mineral extraction sites (0.2%)
    %        2) Burnt areas (0.1%)
    %        3) Beaches - dunes - sands (>0.05%)
    %        4) Sparsely vegetated areas (0.7%)
    %% CODE
        POI.VEG_CLASS(k) = 10;

        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 0.9;
               
        POI.II(k) =  {["grass_Alinya"]}; 
        POI.Ccrown(k) = {[0.1]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));              
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));

        otherwise
        disp('INDEX FOR VEGETATION PARAMETER INCONSISTENT')
        %return
        end
        
    else % if I have an station, this is in bare soil.    
      
        POI.VEG_CLASS(k) = 10;

        POI.Cwat(k) = 0.0; POI.Curb(k) = 0.0 ; POI.Crock(k) = 0.0; POI.Cbare(k) = 1.0;
                 
        POI.II(k) =  {["grass_Alinya"]};  
        POI.Ccrown(k) = {[0.0]};

        POI.NCrown(k) = length(cell2mat(POI.Ccrown(k)));      
        POI.cc_max(k) = length(cell2mat(POI.Ccrown(k)));

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
writetable(POI, [Directories.save '6_POI_table/POI_' run_folder '.csv']);


%% SINGLE POINT LAUNCHER FOR PERSONAL COMPUTER
%==========================================================================

% Time counter
%--------------------------------------------------------------------------
tic;

%Main function
%--------------------------------------------------------------------------
run_Point_Pro_BC(Directories, run_folder, "VelinoCluster100", POI, ksv, date_start, date_end, TT_par, zatm_surface, Clusters);

% Computational time
%--------------------------------------------------------------------------
Computational_Time =toc;
disp('Computation time - Single Point')
disp([num2str(round(Computational_Time/60,1)) ' mins'])

% Memory out
%--------------------------------------------------------------------------
%profile off
%profile viewer


%% MULTIPLE POINT LAUNCHER 
%==========================================================================
% Specify the desired number of workers
%==========================================================================

% Workers - Only for personal computer
%--------------------------------------------------------------------------

%{
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



% Computational time - Time counter - For personal computer or HPC
%--------------------------------------------------------------------------
tic;

% Specific tasks - Here define the runs
%--------------------------------------------------------------------------

idx = contains(names, "Cluster");
names = names(idx)

specific_indices = 1:length(names); % for all the runs
%specific_indices = [1, 23];
%specific_indices = [5, 11, 14, 22, 28, 30, 36, 38, 41, 55, 59, 65, 68, 76, 77, 82, 84, 90, 92, 95];
%specific_indices = [4, 10, 14, 25, 36, 80];

% Create a new cell array of names based on the specific indices
selected_names = names(specific_indices);

% Parallel computing launch for the Model
%--------------------------------------------------------------------------
parfor k = 1:length(selected_names) %length(names)  % 3    
    try
        run_Point_Pro_BC(root, outlocation, run_folder, selected_names(k), POI, ksv, date_start, date_end, TT_par, zatm_surface, Clusters);
    catch ME
        warning('Error occurred on worker %d: %s', k, ME.message);       
    end
end

% Computational time - End of processing
%--------------------------------------------------------------------------
Computational_Time =toc;
disp('COMPUTATIONAL TIME PARFOR [h] ')
disp(Computational_Time/3600)


%}

%% Memory use
[user, sys] = memory; % Windows only
disp(['Mem used by worker: ', num2str(user.MemUsedMATLAB/1e6), ' MB'])
