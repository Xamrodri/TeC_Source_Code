%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE DEM-BASED DISTRIBUTED INPUTS FOR T&C
%
% (initially adapted from Mike's code, 24 Mar 2020, 
%  "tethys_chloris_dem_products.m")

% PB 2021/08/13
% Modified by AJ and SF since 02.11.2021
%
% NEW:  -inside buffer (1 cell) applied to basin mask to avoid
%        inconsistencies at border regions
%       -write basinmask to shapefile
%
% TODO:  -consistently use topotoolbox for streams, accumarea?
%       
% Pascal Buri | High Mountain Glaciers and Hydrology |
%  Swiss Federal Institute for Forest, Snow and Landscape Research, WSL |
%  Office MG E 35 | Zürcherstrasse 111, 8903 Birmensdorf |
%  pascal.buri@wsl.ch
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all

%==========================================================================
% SETTINGS
%==========================================================================

%define catchment name and target spatial resolution of all rasters [m]
targres = 100;
catchment_name = 'Kyzylsu';
username = 'jouberto';  %fugger

% Define path to T&C folder

if strcmp(username,'fugger')
    path_tc_git = 'C:\Users\fugger\TC_Setups\'; % Path to git repository 
    path_tc = path_tc_git; % 
    fnPOI = 'POI/POIs_KYZ_input.csv';
elseif strcmp(username,'jouberto')
    path_tc_git = 'C:\Users\jouberto\Desktop\T&C\TC_setups\'; % Path to git repository 
    path_tc = path_tc_git; %   
    fnPOI = 'POI/Kyzylsu_POIs.csv';
    root ='C:\Users\jouberto\Desktop';
end 

path_sp =  ['B:\home\Study_sites\' ...
    catchment_name '\Spatial_data\']; % Path of folder where spatial data is stored

%define site name
SIMnm =   [catchment_name '_' num2str(targres) 'm'];

%define name of outlet coordinate (in POI-file)
outlet_nm = 'Gauge4';
gauge_nm = ["Gauge1", "Gauge2", "Gauge3"];%define table of POI for T&C-output

%define DEM filename
 fnDem = 'DEMs\KYZYLSU_DEM_T_C.tif';

%define glaciers & debris shapefiles
fnGlacShp = 'Shapefiles\glacier_Kyzylsu_region_RGI60_UTM42N.shp';
fnDebShp = 'Shapefiles\manual_deliniated\DCG.shp';

%define debris thickness file
fnDebTh = 'DebrisThickness\Kyzylsu_modDT_BT_elev_utm.tif';

%define glacier thickness file
 fnGlaTh = 'GlacierThickness\Thickness_composite_KyzylsuRegion_Farinotti2019_utm.tif'; % Farinotti 2019

%define vegetation map file
fnVeg = 'Vegetation\Landcover_PROBA-V\PROBAV_LC100_2019_KYZ_UTM42N.tif';

%define soil map files
fnPSAN = 'Soil\SOILGRIDS\AVGPERC_Sand_0-100cm_KyzylsuRegion_UTM42N.tif';
fnPCLA = 'Soil\SOILGRIDS\AVGPERC_ClayContent_0-100cm_KyzylsuRegion_UTM42N.tif';
fnPORG = 'Soil\SOILGRIDS\AVGPERC_SoilOrganicCarbon_0-100cm_KyzylsuRegion_UTM42N.tif';

%define snow depth map file
% fnSND = 'SnowDepth\SND_KYZ_2018-09-09_v1.tif';
fnSND = 'SnowDepth\Landsat_SnowCoverHybrid_Kyzylsu_2014-08-30.tif';

%==========================================================================
% PATHS
%==========================================================================
%path to GIS-data
path_gisdata = path_sp;
cd(path_gisdata)

%path to functions and scripts
path_func = [path_tc_git 'Functions'];
addpath(path_func)
addpath([path_tc_git '\Scripts'])

%path to topotoolbox
addpath(genpath([path_tc_git 'Functions\topotoolbox']))

%path to store output fil
path_out = [root '\T&C\TC_preprocessing\' catchment_name '\Preprocessing\OUTPUTS\Outputs_' datestr(datetime("today")) '_' num2str(targres) 'm'];
path_out1 = [path_tc catchment_name '\RUNS\INPUTS'];

%path to filled Hugonnet 2021 GMB maps
path_GMB = 'GMB\Hugonnet_2021\Filled';

%path to figures
path_fig = [path_out '\Figures'];

% Create output folder if it doesn't exist

if ~exist(path_out, 'dir')
    mkdir(path_out); addpath(path_out)
    mkdir(path_fig)
    mkdir([path_out '\TIF'])
    mkdir([path_out '\Mask'])
end 

%define x, y limits to crop DEM to [xmin,xmax,ymin,ymax]
% (has to be the same limits as in "prepare_INPUT_PrecipitaionMap.r")

xylims_wgs = [71.38, 71.56, 38.98, 39.20];
xylims([1,3]) = ll2utm(xylims_wgs(3),xylims_wgs(1),42);
xylims([2,4]) = ll2utm(xylims_wgs(4),xylims_wgs(2),42);

%%
%==========================================================================
% LOAD, RESAMPLE & CROP DEM
%==========================================================================
DEM = GRIDobj(fnDem);
% DEM = resample(DEM,targres,'nearest'); No need since already 100m resolution
DEM = crop(DEM,xylims(1:2),xylims(3:4));
DEM.Z(DEM.Z == 0) = NaN;

Zone = string(DEM.georef.GeoKeyDirectoryTag.GTCitationGeoKey);
Zone_num = extractAfter(Zone,'WGS 84 / UTM zone ');
Zone_num = str2double(regexp(Zone_num,'\d*','Match'));
EPSG = string(DEM.georef.GeoKeyDirectoryTag.ProjectedCSTypeGeoKey);

%Plot DEM
figure()
imageschs(DEM)
colorbar
title('Elevation (m)')
saveas(gcf,[path_fig '\' SIMnm '_dem.png'])
close(gcf)

GRIDobj2geotiff(DEM,[path_out  '\TIF\' SIMnm '_DEM_' char(EPSG) '.tif'])

%==========================================================================
% FILL SINKS, DERIVE FLOW DIRECTION & CALCULATE UPSLOPE AREA
%==========================================================================
DEM.Z(isnan(DEM.Z)) = 0;
DEMmat = fill_sinks(DEM.Z);
flowDirection = dem_flow(DEMmat);

%Plot flow direction
figure()
imagesc(flowDirection)
colorbar;
title('Flow direction (rad)')
saveas(gcf,[path_fig '/' SIMnm '_flow_direction.png'])
close(gcf)

%Get sparse matrix of flow direction
flowMatrix = flow_matrix(DEMmat,flowDirection);
clear flowDirection %%recalculated below (upside down)

%Plot sparse matrix of flow direction
figure()
spy(flowMatrix)
title('Flow matrix')
saveas(gcf,[path_fig '/' SIMnm '_flow_matrix.png'])
close(gcf)

%Calculate upslope area
upslopeArea = upslope_area(DEMmat,flowMatrix);

%Plot upslope area
figure()
imagesc(log(upslopeArea))
colorbar
title('ln(upslope area) (pixels)')
saveas(gcf,[path_fig '/' SIMnm '_upslope_area.png'])
close(gcf)


%==========================================================================
% DERIVE STREAM NETWORK
%==========================================================================
%Get stream network such that some percentage of DEM is stream
maxUA = max(max(upslopeArea));
minUA = min(min(upslopeArea));
uAThresholds = linspace(maxUA,minUA,500);
nUAThresholds = length(uAThresholds);
totalPixels = numel(upslopeArea);
streamPerc = nan(1,nUAThresholds);
for iUAThreshold = 1:nUAThresholds
    streamPixels = upslopeArea > uAThresholds(iUAThreshold); %  & mask == 1
    streamPerc(iUAThreshold) = 100*sum(sum(streamPixels))/totalPixels;
end
optPerc = 2;
streamPercMinus5 = abs(streamPerc-optPerc);
[~,optUAThresholdInd] = min(streamPercMinus5);
optUAThreshold = uAThresholds(optUAThresholdInd);
optStreamPerc = streamPerc(optUAThresholdInd);

%Plot threshold optimisation
figure()
plot(uAThresholds,streamPerc); hold on
scatter(optUAThreshold,optStreamPerc)
text(0.95,0.9,['UAT = ' num2str(optUAThreshold,1) ', S% = ' ...
    num2str(optStreamPerc,2)],'Units','Normalized',...
    'HorizontalAlignment','right');
xlabel('Upslope area threshold (pixels)')
ylabel('Area (%)')
title('Optimise stream network')
saveas(gcf,[path_fig '/' SIMnm '_optimise_stream_network.png'])
close(gcf)

%Define stream network
streamNetwork = upslopeArea > optUAThreshold;
streamNetwork_tif = DEM;
streamNetwork_tif.Z = streamNetwork;

%Plot stream network
figure()
imagesc(streamNetwork)
colorbar
title('Stream network')
saveas(gcf,[path_fig '/' SIMnm '_stream_network.png'])
close(gcf)

GRIDobj2geotiff(streamNetwork_tif,[path_out  '\TIF\' SIMnm '_stream_network.tif'])

%%
%==========================================================================
% READ POI COORDINATES (DEGREES) & CONVERT TO UTM
%==========================================================================
poi = readtable(fnPOI); %read table with coordinates
if poi.x > 180
    poi_utm = poi;
    % convert to DEG
    poi_deg = poi;
    poi_deg(:,2:3) = num2cell(utm2ll(poi.x,poi.y,Zone_num,'wgs84'));
end
if poi.x < 180
    poi_deg = poi;
    % convert to UTM
    X_v = double(poi.x);
    Y_v = double(poi.y);
    [X_v,Y_v,~] = deg2utm(Y_v,X_v);  % deg2utm(LAT,LONG) !
    poi_utm = poi;
    poi_utm.x = X_v;
    poi_utm.y = Y_v;
end

%==========================================================================
% DROP POI COORDINATES OUTSIDE CROPPED DEM
%==========================================================================
n = height(poi_utm);
in_idx = zeros(0,n);
for i = 1:n
    x = poi_utm.x(i);
    y = poi_utm.y(i);
    if (x > xylims(1)) && (x < xylims(2)) && (y > xylims(3)) && (y < xylims(4))
        in_idx(i) = 1;
    else
        in_idx(i) = 0;
    end
end
in_idx = find(in_idx == 1);
poi = poi(in_idx,:) ;
poi_utm = poi_utm(in_idx,:);
poi_deg = poi_deg(in_idx,:);


%==========================================================================
% FIND POI AT OUTLET FOR BASIN DELINEATION
%==========================================================================
out_idx = find(string(poi.Name) == outlet_nm);
for i=1:numel(gauge_nm) %for some reason the indexing doesn't work without the loop
gauge_idx(i) = find(string(poi.Name) == gauge_nm(i));
end
clear('X_v','Y_v','poi');

%%
%==========================================================================
% DERIVE GLACIER MASK and IDs
%==========================================================================
glaciers = shaperead(fnGlacShp);

for i = 1:size(glaciers,1)
    RGI_Id_str = glaciers(i).RGIId;
    RGI_Id(i) = str2num(RGI_Id_str([7:8 10:end]));
end
C = num2cell(RGI_Id');
[glaciers.('RGIID_num')] = C{:};

glaciermaskid = polygon2GRIDobj(DEM,glaciers,'RGIID_num');  %expensive
glaciermaskid.Z(isnan(glaciermaskid.Z)) = 0;

glaciermask = polygon2GRIDobj(DEM,glaciers);  %expensive

%Plot glaciers
figure()
imageschs(DEM,glaciermask)
colorbar
title('Glacier mask')
saveas(gcf,[path_fig '/' SIMnm '_glaciers.png'])
close(gcf)

%==========================================================================
% READ DEBRIS MASK
%==========================================================================
debrismask = shaperead(fnDebShp);
debrismask = polygon2GRIDobj(DEM,debrismask);
debrismask = resample(debrismask,DEM,'bilinear');
% mask debris mask with glacier mask to avoid "Error: Debris without ice" in T&C
debrismask.Z(glaciermask.Z < 1) = 0; 
%==========================================================================
% GET GLACIER THICKNESS MAP
%==========================================================================
glThick = GRIDobj(fnGlaTh);
glThick = resample(glThick,DEM,'bilinear');
glThick = crop(glThick,xylims(1:2),xylims(3:4));

%mask glacier thickness map with glacier mask
glThick.Z(glaciermask.Z < 1) = NaN;
%add value of 10m ice thickness where glThick = 0 but debris-covered
glThick.Z(isnan(glThick.Z) & debrismask.Z == 1) = 10;
% set non-glacier to 0 (NaN cause errors in T&C (PARAMETERS_SOIL_Gletsch.m)
glThick.Z(isnan(glThick.Z)) = 0; 

%Plot glacier thickness
figure()
imagesc(glThick); caxis([0 300])
colorbar
title('Ice thickness (m)')
saveas(gcf,[path_fig '/' SIMnm '_ice_thickness.png'])
close(gcf)
%Output final debris thickness map als geotif
GRIDobj2geotiff(glThick,[path_out  '\TIF/' SIMnm '_ice_thickness.tif'])

% %==========================================================================
% % GET HUSS PARAMETERS
% %==========================================================================

invert_Huss_parameters

%==========================================================================
% GET DEBRIS THICKNESS MAP
%==========================================================================
dbThick_raw = GRIDobj(fnDebTh);
% %resample first to higher resolution than target resolution to avoid
% %inaccuracies with glaciermask ("Error: Debris without Ice" in T&C)
% %crop with buffer

% dbThick = crop(dbThick,xylims(1:2)+1000,xylims(3:4)+1000); 
% if targres >= 50
%     dbThick = resample(dbThick,10,'nearest');  
% end

dbThick_raw = resample(dbThick_raw,DEM,'bilinear');
% dbThick = resample(dbThick,DEM,'nearest');
dbThick_raw = crop(dbThick_raw,xylims(1:2),xylims(3:4));
dbThick_raw = dbThick_raw*10; %convert from cm to mm
dbThick = dbThick_raw;

%swap in measured debris thicknesses
dtMeas = readtable('DebrisThickness\DT_measurements.csv');
R = DEM.georef.SpatialRef;
res = 100;
nPts = size(dtMeas,1); %249
[x,y] = ll2utm(dtMeas.lat,dtMeas.lon,42);

r = nan(nPts,1);
c = nan(nPts,1);
[heightDTM,widthDTM] = size(DEMmat);

repl_tab = array2table(nan(nPts,2));
repl_tab.Properties.VariableNames = {'repl','ij'};
for iPt = 1:nPts
    
    % Convert to row and column
    rc = getrowcolfromdem(R,res,[x(iPt),y(iPt)]);
    r(iPt) = rc(1);
    c(iPt) = rc(2);
    
    % Get ij
    ij = sub2ind([heightDTM widthDTM],r(iPt),c(iPt));
    repl_tab.repl(iPt) = dtMeas.thickness_cm_(iPt)*10;
    repl_tab.ij(iPt) = ij;
end
%find median measurement per cell and swap it in
ijs = unique(repl_tab.ij);
for i = 1:length(ijs)
repls = repl_tab.repl(repl_tab.ij == ijs(i));
med = median(repls);% Swap in debris thicknesses
dbThick.Z(ijs(i)) = med;
end

%check debris before/after
figure()
subplot(1,2,1)
imagesc(dbThick_raw)
colorbar
caxis([0 1000])
subplot(1,2,2)
imagesc(dbThick)
colorbar
caxis([0 1000])

%Cap upper and lower limit debris thickness
dbThick.Z(dbThick.Z > 1000) = 1000;
dbThick.Z(dbThick.Z < 30 & debrismask.Z ==1 ) = 30;

% idx = find(glaciermask.Z == 0 & debrismask.Z == 1);

%mask debris thickness map with glacier masks
dbThick.Z(debrismask.Z < 1) = NaN;
%dbThick.Z(glaciermask.Z < 1) = NaN;  
dbThick.Z(glThick.Z==0) = 0; 
% set non-debris to 0 (NaN cause errors in T&C (PARAMETERS_SOIL_Gletsch.m)
dbThick.Z(isnan(dbThick.Z)) = 0; 

% %%%CHECKING
% imagesc(glaciermask.Z), colorbar
% imagesc(dbThick.Z), colorbar
% tt=glaciermask.Z*0;
% tt(glaciermask.Z == 1 & dbThick.Z == 0)=1;
% imagesc(tt), colorbar

%Plot debris thickness
figure()
imagesc(dbThick)
colorbar
hold on
title('Debris thickness (mm)')
saveas(gcf,[path_fig '/' SIMnm '_debris_thickness.png'])
close(gcf)

%Output final debris thickness map als geotif
GRIDobj2geotiff(dbThick,[path_out  '\TIF/' SIMnm '_debris_thickness.tif'])

%==========================================================================
% DERIVE DRAINAGE BASIN
%==========================================================================
%check location of outlet
DEM_bedrock = DEM;
DEM_bedrock.Z = DEM.Z - glThick.Z;

streams = (DEM*0 + streamNetwork);
streams.Z = logical(streams.Z);
figure()
[idx,d] = snap2stream(streams,table2array([poi_utm(out_idx,'x'),poi_utm(out_idx,'y')]),true);
title(['Outlet repositioning (d=',num2str(d),'m)'])
axis on
saveas(gcf,[path_fig '/' SIMnm '_outpos_check.png'])
close(gcf)
figure()
[idx_gauge,d_gauge] = snap2stream(streams,table2array([poi_utm(gauge_idx,'x'),poi_utm(gauge_idx,'y')]),true);
title(['Gauge repositioning (d=',num2str(d),'m)'])
axis on
saveas(gcf,[path_fig '/' SIMnm '_gaugepos_check.png'])
close(gcf)
[r,c] = size(streams.Z);
[Y_Outlet,X_Outlet] = ind2sub([r,c],idx);
[Y_Gauge,X_Gauge] = ind2sub([r,c],idx_gauge);
X_Outlet = X_Outlet;  % Shift so that the outlet is on a single pixel width flow line
Y_Outlet = Y_Outlet;
clear idx idx_gauge


DEM = fillsinks(DEM); %fill sinks in DEM before basin delineation
DEM_bedrock = fillsinks(DEM_bedrock); %fill sinks in DEM before basin delineation
FD = FLOWobj(DEM_bedrock,'preprocess','c'); %flow direction object
[X,Y] = sub2coord(DEM_bedrock,Y_Outlet,X_Outlet);
D = drainagebasins(FD,X,Y);
A  = flowacc(FD);
clear('X','Y');

% apply inside buffer to basin mask (to avoid inconsistencies at basin-borders)
D2=D;
for i = 2:(size(D2.Z,1)-1)
    for j = 2:(size(D2.Z,2)-1)
        % for i = 2:150
        %     for j = 2:150
        if double(D2.Z(i,j)) == 1
            vect=[D.Z(i-1,j) ...
                D.Z(i+1,j) ...
                D.Z(i,j-1) ...
                D.Z(i,j+1)];
            idx=find(vect == 0);
            if length(idx) > 0
                % D2.Z(i,j)=2; %%for checking only
                D2.Z(i,j)=0;
            else
                D2.Z(i,j)=1;
            end
            clear vect idx
        end
    end
end
%%checking
% imagesc(D2.Z); figure(2); imagesc(D.Z)
D=D2;
clear D2
D.Z(Y_Outlet,X_Outlet)=1; %make sure outlet cell is = 1

% % Load basin shapefile
% Cat_outline = shaperead(fnCat);
% Dnew = polygon2GRIDobj(DEM, Cat_outline);

% create basin shapefile & write to file
[Bas,bas_x,bas_y]=GRIDobj2polygon(D,'Geometry','Polygon','simplify',true,'tol',0);
% [Bas,bas_x,bas_y]=GRIDobj2polygon(D,'Geometry','Line','simplify',true,'tol',0);
Bas.Name='Basin';
% tt=polybuffer([bas_x,bas_y],'lines',100);
fn = [path_out '/Basin_' SIMnm '.shp'];
shapewrite(Bas,fn);    % Bas=shaperead(fn);
clear fn

%Plot basin
figure()
imageschs(DEM)
hold on
plot(bas_x,bas_y,'-r')
colorbar('off')
title('Drainage basin')
saveas(gcf,[path_fig '/' SIMnm '_basin.png'])
close(gcf)

%Plot flow accumulation
figure()
imageschs(DEM,log(A))
hold on
plot(bas_x,bas_y,'-r')
title('Flow accumulation (pixel)')
saveas(gcf,[path_fig '/' SIMnm '_flowacc.png'])
close(gcf)

%mask DEM with basin
D.Z = double(D.Z);
DEM_basin = DEM * D;
DEM_basin.Z(DEM_basin.Z == 0) = NaN; %zeros to NaNs (for T&C)

%create mask
MASK = double(GRIDobj2mat(D)); %if simulation should be for watershed (=1)

%==========================================================================
% VEGETATION MAP
%==========================================================================
% doc topotoolbox
ext = getextent(DEM);
Veg_map = GRIDobj(fnVeg);
Veg_map = resample(Veg_map,DEM,'nearest'); %'nearest' important for discrete data!
Veg_map = crop(Veg_map,ext(1:2),ext(3:4));
Veg_map.Z(isnan(Veg_map.Z)) = 0;
clear('ext');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%T&C VEGETATION CODE Kyzylsu %%%
%%% 0 Outside mask
%%% 1 Fir  (evergr.)
%%% 2 Larch (decid.)
%%% 3 Grass C3
%%% 4 Shrub (decid.)
%%% 5 Broadleaf evergreen
%%% 6 Broadleaf deciduous
%%% 7 Rock
%%% 8 Ice
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%PROBA-V LANDCOVER MAP %%%
%%%      20 = Shrubs -> 4
%%%      30 = Herbaceous vegetation -> 3
%%%      40 = Cultivated and managed vegetation/agriculture (cropland) -> 3
%%%      50 = Urban / built up -> 7
%%%      60 = Bare / sparse vegetation (exposed soil, sand, or rocks and
%%%      never has more than 10 % vegetated cover during any time of the
%%%      year) -> 7
%%%      70 = Snow and Ice (throughout the year) -> 8
%%%      80 = Permanent water bodies -> 3 (for now)
%%%      90 = Herbaceous wetland -> 3
%%%     100 = Moss and lichen -> 3
%%%     110 = closed forest, deciduous needle leaf -> 2
%%%     111 = closed forest, evergreen needle leaf -> 1
%%%     112 = closed forest, evergreen broad leaf -> 5
%%%     114 = closed deciduous broad leaf forest -> 6
%%%     115 = closed forest, mixed -> 6
%%%     116 = closed forest, unknown -> 6
%%%     121 = open forest, evergreen needle leaf -> 1
%%% 	122 = open forest, evergreen broad leaf -> 5
%%%     123 = open forest, deciduous needle leaf -> 2
%%%     124 = open forest, deciduous broad leaf -> 6
%%%     125 = open forest, mixed -> 6
%%%     126 = Open forest, unknown -> 1 

Veg_map.Z(Veg_map.Z == 111  | Veg_map.Z == 121 ) = 1;  %Fir
Veg_map.Z(Veg_map.Z == 110 | Veg_map.Z == 113 |  Veg_map.Z == 115 | Veg_map.Z == 123 | Veg_map.Z == 126) = 2;  %Larch
Veg_map.Z(Veg_map.Z == 30 | Veg_map.Z == 40 | Veg_map.Z == 80 | Veg_map.Z == 90 | Veg_map.Z == 100 ) = 3;  %Grass
Veg_map.Z(Veg_map.Z == 20) = 4;  %Shrub (decid.)
Veg_map.Z(Veg_map.Z == 112 | Veg_map.Z == 122 ) = 5;  %Broadleaf (evergreen.)
Veg_map.Z(Veg_map.Z == 114 | Veg_map.Z == 115 | Veg_map.Z == 116 | Veg_map.Z == 124 | Veg_map.Z == 125) = 6;  %Broadleaf (decid.)
Veg_map.Z(Veg_map.Z == 50 | Veg_map.Z == 60) = 7;  %Rock
Veg_map.Z(Veg_map.Z == 70) = 8;  %Ice

%%%Corrections of vegetation map
%%% 1) convert ice (=8) into rock (=7) and only glaciermask into ice again
Veg_map.Z(Veg_map.Z == 8) = 7; 
Veg_map.Z(glaciermask.Z == 1) = 8;
% %%% 2) convert grass/shrub/forest (<7) above 4500m a.s.l. into rock (=7) 
% Veg_map.Z(Veg_map.Z < 7 & DEM_basin.Z > 4500) = 7;
%%% 3) convert no data values (=0) within basin to rock (=7)
Veg_map.Z(Veg_map.Z == 0) = 7; 


%Plot vegetation map
figure()
imagesc(Veg_map)
colorbar
hold on
plot(bas_x,bas_y,'-r')
title('Vegetation code')
saveas(gcf,[path_fig '/' SIMnm '_vegetation_code.png'])
close(gcf)


%==========================================================================
% SOIL PROPERTIES MAP
%==========================================================================
ext = getextent(DEM);
%Map with percent SAND in soil profile
PSAN_map = GRIDobj(fnPSAN);   
PSAN_map = resample(PSAN_map,DEM,'nearest'); %'nearest' important for discrete data!
PSAN_map = crop(PSAN_map,ext(1:2),ext(3:4));
%Map with percent CLAY in soil profile
PCLA_map = GRIDobj(fnPCLA);   
PCLA_map = resample(PCLA_map,DEM,'nearest'); 
PCLA_map = crop(PCLA_map,ext(1:2),ext(3:4));
%Map with percent ORGANIC in soil profile
PORG_map = GRIDobj(fnPORG);   
PORG_map = resample(PORG_map,DEM,'nearest');
PORG_map = crop(PORG_map,ext(1:2),ext(3:4));
clear ext;

%%%Corrections of soil maps (consistency with glacier mask and vegetation map
%%% 1) convert NaNs to 0 (no soil)
PSAN_map.Z(isnan(PSAN_map.Z)) = 0;
PCLA_map.Z(isnan(PCLA_map.Z)) = 0; 
PORG_map.Z(isnan(PORG_map.Z)) = 0;

%%% 2) convert rock (VEG=7) & glaciermask to 0 (no soil)
PSAN_map.Z(Veg_map.Z == 7) = 0; PSAN_map.Z(glaciermask.Z == 1) = 0;
PCLA_map.Z(Veg_map.Z == 7) = 0; PCLA_map.Z(glaciermask.Z == 1) = 0;
PORG_map.Z(Veg_map.Z == 7) = 0; PORG_map.Z(glaciermask.Z == 1) = 0;

%%% 3) convert vegetated cells (VEG<7) without soil (SOIL=0) to soil cell
%%% with properties from nearest soil cell >0

% Index for cells with vegetation BUT NO soil 
ix_torepl = find(and((PSAN_map.Z+PCLA_map.Z+PORG_map) == 0, and(Veg_map.Z >0, Veg_map.Z <7)));
% check: unique(PSAN_map.Z(i3)); unique(PCLA_map.Z(i3)); unique(PORG_map.Z(i3)); unique(Veg_map.Z(i3))
%Index for cells with vegetation AND soil
ix_todupl = find((PSAN_map.Z+PCLA_map.Z+PORG_map) > 0);

[x_torepl,y_torepl] = ind2coord(PSAN_map,ix_torepl);
[x_todupl,y_todupl] = ind2coord(PSAN_map,ix_todupl);
coord_todupl = [x_todupl,y_todupl];

% find nearest cell with soil > 0 and duplicate to cell with soil=0
for i=1:length(ix_torepl)
    coord_torepl = [x_torepl(i) y_torepl(i)];
    %compute Euclidean distances
    distances = sqrt(sum(bsxfun(@minus, coord_todupl, coord_torepl).^2,2)); 
    %find the smallest distance and use that as an index:
    closest = coord_todupl(find(distances==min(distances)),:);
    ix = coord2ind(PSAN_map, closest(1,1), closest(1,2)); %take only first coordinate
    PSAN_map.Z(ix_torepl(i)) = PSAN_map.Z(ix);
    PCLA_map.Z(ix_torepl(i)) = PCLA_map.Z(ix);
    PORG_map.Z(ix_torepl(i)) = PORG_map.Z(ix);
    clear coord_torepl distances closest ix;
end
% check: min(PSAN_map.Z(ix_torepl)); min(PCLA_map.Z(ix_torepl)); min(PORG_map.Z(ix_torepl))
% add soil under glacier, 100% sand
PSAN_map.Z(glaciermask.Z==1) = 50;
PCLA_map.Z(glaciermask.Z==1) = 10;

%Plot soil maps
figure()
imagesc(PSAN_map)
colorbar
hold on
plot(bas_x,bas_y,'-r')
title('Sand content [%]')
saveas(gcf,[path_fig '/' SIMnm '_soil_sand.png'])
close(gcf)

figure()
imagesc(PCLA_map)
colorbar
hold on
plot(bas_x,bas_y,'-r')
title('Clay content [%]')
saveas(gcf,[path_fig '/' SIMnm '_soil_clay.png'])
close(gcf)

figure()
imagesc(PORG_map)
colorbar
hold on
plot(bas_x,bas_y,'-r')
title('Organic content [%]')
saveas(gcf,[path_fig '/' SIMnm '_soil_org.png'])
close(gcf)


%==========================================================================
% SOIL THICKNESS MAP IN %, f(SLOPE)
%==========================================================================
SLP = arcslope(DEM, 'deg'); %% imageschs(DEM,SLP);
ix_soil = find((PSAN_map.Z+PCLA_map.Z+PORG_map) > 0);
ix_nosoil = find((PSAN_map.Z+PCLA_map.Z+PORG_map) == 0);
% hist(SLP.Z(ix_soil));  %%show elevations with soil

%%% linear decrease of soil depth (absolute soil depth is defined in T&C >> PARAMETERS_SOIL_XX): 
%%%  >=60°: 1% soil depth   |   0°: 100% soil depth
mx_slp = 60;
mx_slth = 2.5;
SLP.Z(SLP.Z > mx_slp) = mx_slp;
SOIL_TH = (mx_slth-SLP/mx_slp)*100; %% imageschs(DEM,SOIL_TH);

%%% no soil where glacier
%SOIL_TH.Z(glaciermask.Z == 1) = 0; 
%%%
%SOIL_TH.Z(ix_nosoil) = 0;  %% can not be NaN, otherwise problems in T&C 
%%%

figure()
imageschs(DEM,SOIL_TH);
colorbar
hold on
plot(bas_x,bas_y,'-r')
title('Soil thickness [%]')
saveas(gcf,[path_fig '/' SIMnm '_soil_thickness.png'])
close(gcf)


%==========================================================================
% SNOW DEPTH MAP (FOR INITIAL CONDITIONS ONLY)
%==========================================================================
ext = getextent(DEM);
SNC_map = GRIDobj(fnSND);
SNC_map = resample(SNC_map,DEM,'nearest'); 
SNC_map = crop(SNC_map,ext(1:2),ext(3:4));
SNC_map.Z = SNC_map.Z>0;

SND_map = DEM; %GRIDobj(fnSND);  
SND_map.Z = DEM.Z.*0 + 0.40;
SND_map.Z(SNC_map.Z == 0) = 0;

figure()
imageschs(DEM,SND_map);
colorbar
hold on
plot(bas_x,bas_y,'-r')
title('Initial snow depth [m]')
saveas(gcf,[path_fig '/' SIMnm '_snow_depth.png'])
close(gcf)

%==========================================================================
% FLOW MATRIX
%==========================================================================
% recalculate flow direction, now upside down
flowDirection_ud = dem_flow(flipud(DEMmat)); %%[rad] from North
% ouside basin = NaN
DEMmat(MASK == 0) = NaN; 
flowDirection_ud(flipud(MASK) == 0) = NaN;
% imagesc(flowDirection*180/pi), axis equal, colorbar

%   T = flow_matrix(E, R) computes a sparse linear system representing flow from
%   pixel to pixel in the DEM represented by the matrix of height values, E.  R
%   is the matrix of pixel flow directions as computed by dem_flow.
%   T is  numel(E)-by-numel(E).  The value T(k,l) is the negative of the fractional
%   flow from pixel l to pixel k, where pixels are numbered columnwise. For
%   example, if E is 15-by-10, then T is 150-by-150, and T(17,18) is the
%   negative of the fractional flow from pixel 18 (row 3, column 2) to pixel 17
%   (row 2, column 2).
%   [i,j] = ind2sub([m,n],l); flow into --->  [i,j] = ind2sub([m,n],k);
T = flow_matrix(flipud(DEMmat),flowDirection_ud);


%==========================================================================
% PREPARE MAP-COORDINATES
%==========================================================================
[~,x_v,y_v] = GRIDobj2mat(DEM);    % x y coordinates
y_v = flip(y_v);    %y has to go from smallest to largest coordinate (S-N)


%==========================================================================
% PREPARE MAP-CELL INDICES FOR POI
%==========================================================================
% %x & y indices to outlet cell (geographic matrix, not flipped)
% [~,X_Outlet] = min(abs(x_v - table2array(poi_utm(out_idx,'x'))));
% [~,Y_Outlet] = min(abs(y_v - table2array(poi_utm(out_idx,'y'))));

%Y-coordinate has to be flipped to be in T&C format
% correct cell coordinate can be checked here:
% fldem = flipud(DEM.Z)
% imagesc(fldem);
% axis xy
Y_Outlet = r - Y_Outlet + 1;
Y_Gauge = r - Y_Gauge + 1;
clear streams idx d

%x & y indices to all POI cells ("tracked points")
X_v = table2array(poi_utm(:,'x'));
Y_v = table2array(poi_utm(:,'y'));
[npoi,~] = size(poi_utm);
X_POI = zeros(npoi,1); Y_POI = zeros(npoi,1);
for i = 1:npoi
    if i == out_idx
        X_POI(i) = X_Outlet;
        Y_POI(i) = Y_Outlet;
    elseif any(gauge_idx==i)
        idx_g = find(gauge_idx == i);
        X_POI(i) = X_Gauge(idx_g);
        Y_POI(i) = Y_Gauge(idx_g);
    else
    [~,X_POI(i)] = min(abs(x_v - X_v(i)));
    [~,Y_POI(i)] = min(abs(y_v - Y_v(i)));
    end
end
clear('X_v','Y_v');

% figure()
% imageschs(DEM)
% hold on
% plot(bas_x,bas_y,'-r');
% scatter(x_v(X_POI), y_v(Y_POI),20,'+r','LineWidth',2)
% text(x_v(X_POI), y_v(Y_POI),string(poi_deg.Name))
% colorbar('off')
% title('Drainage basin')


figure()
imageschs(DEM,log(A))
hold on
plot(bas_x,bas_y,'-r')
title('Flow accumulation (pixel)')
scatter(x_v(X_POI), y_v(Y_POI),20,'+r','LineWidth',2)
text(x_v(X_POI), y_v(Y_POI),string(poi_deg.Name))
colorbar('off')
title('Drainage basin')
%==========================================================================
% FILL ALL VARIABLES FOR T&C INPUT
%==========================================================================
% (MASK: defined earlier already)
Aacc = upslopeArea; %Upstream area
Aacc(MASK == 0) = NaN;
basin_ln = Bas;      %Basin outline
cellsize = round(DEM.cellsize); %Pixel size of DEM (m)
DEB_MAP = double(GRIDobj2mat(dbThick)); %Debris thickness map (mm)
DTM = double(GRIDobj2mat(DEM_basin)); %basin-masked DEM (m)
DTM_orig = double(GRIDobj2mat(DEM)); %entire DEM (m)
GLA_MAP2 = double(GRIDobj2mat(glaciermask)); %Glacier ID map (0=no gl, 1,.. for each gl.)
GLA_ID = double(GRIDobj2mat(glaciermaskid)); %Glacier ID map (0=no gl, 1,.. for each gl.)
GLH = double(GRIDobj2mat(glThick)); %Glacier thickness (m)
[m,n] = size(DEM.Z); %DEM rows and columns
POI_names = string(char(poi_utm.Name)); %names of POI
PSAN = double(GRIDobj2mat(PSAN_map)).*MASK; % Soil: sand percentage map
PCLA = double(GRIDobj2mat(PCLA_map)).*MASK; % Soil: clay percentage map
PORG = double(GRIDobj2mat(PORG_map)).*MASK; % Soil: organic material percentage map
SN = streamNetwork*1; % Stream network... check NaNs
SN(MASK == 0) = NaN;
SNOWD = double(GRIDobj2mat(SND_map)); %Snow depth map (m)
SNOWD(isnan(SNOWD)) = 0; % correct for NaNs, otherwise error in HYPERION!
SOIL_TH = double(GRIDobj2mat(SOIL_TH)).*MASK; % Soil: relative thickness based on slope (%)
SOIL_TH(isnan(SOIL_TH)) = 0; % correct for NaNs, otherwise error in HYPERION!
T_flow = T; %Sparse matrix of flow direction
VEG_CODE = double(GRIDobj2mat(Veg_map)); %Vegetation ID map
x = x_v;    % x coordinates
y = y_v;    % y coordinates
xllcorner = min(x_v)-cellsize/2; % Bottom left corner x coordinate
yllcorner = min(y_v)-cellsize/2; % Bottom left corner y coordinate
Xoutlet = X_Outlet; %X-index to outlet cell
Youtlet = Y_Outlet; %Y-index to outlet cell
Xout = X_POI; %X-indices to POI cells
Yout = Y_POI; %Y-indices to POI cells


%==========================================================================
% CONVERT RASTERS ETC. TO "UPSIDE DOWN"
%==========================================================================
Aacc = flipud(Aacc);
DEB_MAP = flipud(DEB_MAP);
DTM = flipud(DTM);   %check:   imagesc(DTM),axis equal,colorbar()
DTM_orig = flipud(DTM_orig);
GLA_MAP2 = flipud(GLA_MAP2);
GLA_ID = flipud(GLA_ID);
GLH = flipud(GLH);
MASK = flipud(MASK);
PSAN= flipud(PSAN);
PCLA= flipud(PCLA);
PORG= flipud(PORG);
SN = flipud(SN);
SNOWD = flipud(SNOWD);
SOIL_TH = flipud(SOIL_TH);
VEG_CODE = flipud(VEG_CODE);


% %%%CHECKING
% idx=find(POI_names == "HagensColSouth                     ");
% DTM(Yout(idx),Xout(idx))
%
% tt=GLA_MAP2*0;
% tt(GLA_MAP2 == 1 & DEB_MAP == 0)=1;
% tt=GLH*0;
% tt(GLH > 0 & DEB_MAP == 0)=1;
% imagesc(tt), axis xy, colorbar


%==========================================================================
% WRITE GEOREFERENCE TO TXT-FILE
%==========================================================================
T = table(Zone, EPSG);
writetable(T,[path_out '/GEOREF_dtm_' SIMnm '.txt']);
writetable(T,[path_out1 '/GEOREF_dtm_' SIMnm '.txt']);

clear('Zone','T');


%==========================================================================
% WRITE POI-TABLE TO TXT-FILE
%==========================================================================
idx = sub2ind(size(DTM),Yout,Xout);
T = poi_utm;
T.col = Xout;
T.row = Yout;
T.idx = idx;
T.lon = poi_deg.x;
T.lat = poi_deg.y;
writetable(T,[path_out '/POI_dtm_' SIMnm '.txt']);
writetable(T,[path_out1 '/POI_dtm_' SIMnm '.txt']);

clear('idx');

figure()
imageschs(DEM)
hold on
plot(bas_x,bas_y,'-r');
scatter(T.x, T.y,20,'+r','LineWidth',2)
colorbar('off')
title('Drainage basin')
saveas(gcf,[path_fig '/' SIMnm '_basin_POI.png'])
close(gcf)

%==========================================================================
% Add debris-cover to veg map and export to geotiff 
%==========================================================================

% Create a new class to distinguish clean ice and debris covered ice
VEG_classes = unique(VEG_CODE(:));
DEB_class = VEG_classes(end)+1;
VEG_CODE_full = VEG_CODE;
VEG_CODE_full(DEB_MAP>0) = DEB_class;

veg_grid = DEM;
veg_grid.Z = flipud(VEG_CODE_full);
GRIDobj2geotiff(veg_grid,[path_out  '\TIF\' catchment_name '_VEG_CODE_' num2str(round(targres,-1)) 'm_' char(EPSG)  '.tif'])

veg_grid = DEM;
veg_grid.Z = flipud(VEG_CODE_full);
veg_grid.Z(flipud(GLH==0)) = min(veg_grid.Z(flipud(GLH>0)),[],'all');
GRIDobj2geotiff(veg_grid,[path_out  '\TIF\' catchment_name '_VEG_CODE_glaciers_' num2str(round(targres,-1)) 'm_' char(EPSG)  '.tif'])

% Remove glacier for non-glacier TopoSUB
veg_grid_nogla = DEM;
veg_grid_nogla.Z = flipud(VEG_CODE_full);
veg_grid_nogla.Z(flipud(GLH>0)) = min(veg_grid_nogla.Z(flipud(GLH==0)),[],'all');
GRIDobj2geotiff(veg_grid_nogla,[path_out  '\TIF\' catchment_name '_VEG_CODE_NoGlaciers_' num2str(round(targres,-1)) 'm_' char(EPSG)  '.tif'])

%==========================================================================
% Export produced catchment mask to geotiff and shapefile
%==========================================================================

mask_grid = DEM;
mask_grid.Z = flipud(MASK);

GRIDobj2geotiff(mask_grid,[path_out  '\Mask\' catchment_name '_catchment_mask_' num2str(round(targres,-1)) 'm_' char(EPSG)  '.tif'])
GRIDobj2geotiff(mask_grid,[path_out  '\TIF\' catchment_name '_catchment_mask_' num2str(round(targres,-1)) 'm_' char(EPSG)  '.tif'])

mask_shp = GRIDobj2polygon(mask_grid);
shapewrite(mask_shp,[path_out  '\Mask\' catchment_name '_catchment_mask_' num2str(round(targres,-1)) 'm_' char(EPSG)  '.shp'])

% Without glacier mask for TopoSUB
mask_grid_nogla = DEM;
mask_grid_nogla.Z = flipud(MASK & GLH==0);
GRIDobj2geotiff(mask_grid_nogla,[path_out  '\TIF\' catchment_name '_catchment_mask_NoGla_' num2str(round(targres,-1)) 'm_' char(EPSG)  '.tif'])

%==========================================================================
% Export produced glacier mask to geotiff and shapefile
%==========================================================================

gla_grid = DEM;
gla_grid.Z = flipud(GLH>0 & MASK==1);
GRIDobj2geotiff(gla_grid,[path_out '\Mask\' catchment_name '_glacier_mask_' num2str(round(targres,-1)) 'm_' char(EPSG) '.tif'])

gla_shp = GRIDobj2polygon(gla_grid,'multipart',true);
shapewrite(gla_shp, [path_out '\Mask\' catchment_name '_glacier_mask_' num2str(round(targres,-1)) 'm_' char(EPSG) '.shp'])
%==========================================================================
% OVERWRITE DEBRIS THICKNESS MAP AT POIs WHERE OBSERVED THICKNESS IS
%  AVAILABLE
%==========================================================================
%find which POI has measured debris  thickness [m]
% & overwrite debris thickness map
% idx = find(isnan(poi_utm.hDebris) == false);
% for i = 1:length(idx)
%     IDX = idx(i);
%     %only apply meas. debris thickness if not 0 (=outside glaciermask)
%     if DEB_MAP(T.row(IDX),T.col(IDX)) > 0
%         DEB_MAP(T.row(IDX),T.col(IDX)) = poi_utm.hDebris(IDX)*1000;
%     end
%     clear IDX
% end
% clear idx

%==========================================================================
% SAVE EVERYTHING IN T&C-FORMAT
%==========================================================================
fn = [path_out '/dtm_' SIMnm '.mat'];
fn1 = [path_out1 '/dtm_' SIMnm '.mat'];

save(fn,'Aacc','basin_ln','cellsize',...
    'DEB_MAP','DTM','DTM_orig','GLA_MAP2','GLA_ID','GLH','m','MASK','POI_names',...
    'PSAN','PCLA','PORG','n','SN','SNOWD','SOIL_TH','T_flow','VEG_CODE','x','xllcorner',...
    'Xout','Xoutlet','y','yllcorner','Yout','Youtlet','Gla_nEl_nDH')

save(fn1,'Aacc','basin_ln','cellsize',...
    'DEB_MAP','DTM','DTM_orig','GLA_MAP2','GLA_ID','GLH','m','MASK','POI_names',...
    'PSAN','PCLA','PORG','n','SN','SNOWD','SOIL_TH','T_flow','VEG_CODE','x','xllcorner',...
    'Xout','Xoutlet','y','yllcorner','Yout','Youtlet','Gla_nEl_nDH')

clear fn fn1

% % Save this preparing script into the pre-processing output folder

currentfile = matlab.desktop.editor.getActiveFilename;
newbackup= [path_out '/prepare_INPUT_GeoData_' catchment_name '.mat'];
copyfile(currentfile,newbackup);

% % Copy the output folder in the as the OUTPUTS folder (such that OUTPUTS
% always has the most updated pre-processing outputs )

currentfile = path_out;
newbackup = [path_tc catchment_name '\Preprocessing\OUTPUTS'];
copyfile(currentfile,newbackup);

% %==========================================================================
% % CALCULATE SKY-VIEW-FACTOR
% %==========================================================================
% [Slo_top,Aspect]=Slope_Aspect_indexes(DEMmat,cellsize,'mste');
% [HZ,Zasp] = Horizon_Angle(DEMmat,cellsize);
% [SvF,Ct] = Sky_View_Factor(DEMmat,atan(Slo_top)*180/pi,Aspect,HZ,Zasp);
% 
% %Plot sky view factor
% figure()
% imagesc(SvF)
% colorbar
% title('Sky-view factors ()')
% saveas(gcf,[path_fig '/' SIMnm '_sky_view_factors.png'])
% close(gcf)


% %functions
% function rowCol = getrowcolfromdem(demR,pixelSize,xy)
% % Get row and column relative to DEM
% rowCol(1) = floor((demR.YWorldLimits(1,2)+0.5*pixelSize-xy(2))/...
%     pixelSize);
% rowCol(2) = floor((xy(1)-(demR.XWorldLimits(1,1)-0.5*pixelSize))/...
%     pixelSize);
% end