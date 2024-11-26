%% Analyze T&C fully distributed outputs

% Achille Jouberton, achille.jouberton@ista.ac.at
% May 2021 - last updated: 06-09-2024 

clear all; close all; clc

%%%%%  Computer %%%%%¨

Bdrive = 1;
machine=getenv('computername');

if strcmp(machine,'LAPTOP-P75HBCRE') %%% << Laptop >>
    root='C:\Users\jouberto\Desktop';
elseif isempty(machine) %%% << HYPERION >>
    root='\home\jouberto';
end

root_B = 'B:\home';

%==========================================================================
% SETTINGS
%==========================================================================
glacier = 'Kyzylsu';
glacier_NAME = 'KYZ';
gla_id = 1319847; % RGI ID of the main glacier of interest
sim_res = 100; % [m] Model resolution
PP_dataset = 'ERA5Land';

forcing_nm = 'ERA5Land';
sim_nm = '010824_1999_2023_PG00_Tmod0_2000mmSWEcap_a012_c145'; %
sim_nm = [forcing_nm '_' sim_nm];
outlet_nm = 'Gauge4';

monthly_LS_figure = 0; % 1 to plot all landsat scene against T&C output
monthy_modis_figure = 0; % 1 to plot one MODIS scene per month against T&C
modis_validation = 1; % Plot snow cover fraction vs MODIS

labelled_output = 1; % 0 if using old OUTPUT_MANAGER_DIST

% To double check air temperature downscaling
fn_LR ='B:\group\pelligrp\Field_Data_Curated\Kyzylsu\Air_temp\LR_MH_Kyzylsu_2021-07-01_2023-09-13.mat';

% Define path to T&C folder and postprocessing functions
path_tc = [root '\T&C\TC_setups\']; 
addpath(genpath('C:\Users\jouberto\Desktop\T&C\TC'));
path_func = [path_tc 'Functions'];
addpath(path_func)
addpath(genpath([path_tc 'Functions\topotoolbox']))
addpath(genpath('C:\Users\jouberto\Desktop\T&C\Post-processing\Matlab_scripts\Postprocessing_functions'))

load('colorbrewer.mat') % Nice color palette

%path to Landsat/Sentinel-8 maps
path_L8S2 = 'B:\home\Remote_sensing\Landsat\Processed_data\Snow_cover\KYZ_filtered\Kept';

% Path to T&C outputs
if Bdrive == 1
    dir_tcout = [root_B '\TC_outputs\' glacier '\Distributed\' sim_nm];
else 
    dir_tcout = [root '\T&C\TC_outputs\' glacier '\Distributed\' sim_nm];
end 

%Where to store the figures
dir_fig = [root '\T&C\Post-processing\Figures\' glacier '\' sim_nm];

if ~exist(dir_fig, 'dir') || ~exist([dir_fig '\Snow_cover'], 'dir')
    mkdir(dir_fig)
    mkdir([dir_fig '\GMB'])
    mkdir([dir_fig '\Avalanching'])
    mkdir([dir_fig '\Snow_cover'])
end 

% Define paths to T&C pre-processing inputs

path_tcpre = [dir_tcout '\' 'dtm_' glacier '_' num2str(sim_res) 'm.mat']; % dtm file saved with T&C output run 
if exist(path_tcpre,'file') == 0 % In case a previous version of the launcher was used (not saving dtm.file)
    path_tcpre = [root '\T&C\TC_setups\' glacier '\RUNS\INPUTS\dtm_' glacier '_' num2str(sim_res) 'm.mat'];
end

load(path_tcpre,'DTM','x','y','GLH','GLA_ID','POI_names','Xout','Yout','DEB_MAP','VEG_CODE')
MASK = ~isnan(DTM); % catchment mask
mask = ~isnan(DTM(:));
mask_gla = GLH > 0; % glacier mask
[nRows,nCols] = size(DTM);
nCatchPix = sum(~isnan(DTM),'all'); % Numbe of modelled pixel
[xmin_map, xmax_map, ymin_map, ymax_map] = Define_maps_limits(x,y,MASK); % For figures with maps or spatial outputs

%%%%  Load catchment average output table  to extract main datetime vector %%%% 

SPAVG = load_catchment_average(dir_tcout,glacier,labelled_output);

startDate = SPAVG{1,1};
SPAVG_dm = retime(table2timetable(SPAVG),'daily',@nanmean);
SPAVG_ds = retime(table2timetable(SPAVG),'daily',@nansum);
SPAVG_mm = retime(table2timetable(SPAVG),'monthly',@nanmean);
save([dir_fig '\SPAVG_dm'],'SPAVG_dm')
Date = SPAVG.Date;

month_labels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

%%% Load DEM and glacier mask used for T&C and produce glacier mask shapefile 
DEM = GRIDobj([root '\Study_sites\' glacier '\Spatial_data\DEMs\' upper(glacier) '_DEM_T_C.tif']);
demInfo = geotiffinfo([root '\Study_sites\' glacier '\Spatial_data\DEMs\' upper(glacier) '_DEM_T_C.tif']);
utm_zone = demInfo.Zone;
[demLons,demLats] = pixcenters(demInfo.RefMatrix,size(DEM.Z),'makegrid');  

GLA = DEM; GLA.Z = flipud(GLA_ID); % Create a georeferenced glacier masp
gla_shp = GRIDobj2polygon(GLA,'multipart',1,'holes',1); % Create glacier polygons
GRADIENT = gradient8(DEM,'deg'); % Compute slope from DEM
ASPECT = aspect(DEM); % Compute aspect from DEM

%%%% Load catchment glacier outlines (native resolution, not super useful here)

run("Load_catchment_glacier_outlines.m")

%%% Load Hugonnet et al. 2021 GMB values and compute elevation profile
run("load_Hugonnet2021_dhmaps.m")

%%% Load MODIS one day snow cover product
run('import_MODIS_snow_cover.m')

%%% Load Landsat-8 / Sentinel-2 fsc and SLA

opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["dateTime", "se_hybrid", "fsnow_hybrid", "fData"];
opts.VariableTypes = ["datetime", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "dateTime", "InputFormat", "dd-MMM-yyyy");
LS_sno = readtable([path_L8S2 '\' glacier '_fsc_L8S2.csv'], opts);

% Load Monthly T&C spatial  output
run("load_TC_spatial_outputs.m")

nYears = size(SWEm_map,3)./12; % number of years
Years_no = unique(year(Date)); % years

% Extract daily maps of snow depth, SWE, snow density and snow melt

run("load_TC_snow_spatial_outputs.m")

%% Compute simulated daily snow line elevation 

se = NaN(length(date_h),1);
tic
for d = 1:length(date_h)
    se(d) = Snowline_Krajci(DTM,snow_pres(:,:,d));
end 
toc   

%% Compute glacier mass balance for mass balance measurements

load('C:\Users\jouberto\Desktop\T&C\Post-processing\Figures\Kyzylsu\MultiPoints\ERA5Land_090724_2015_2023_newcluster_Ding_Tmod0_Pmod00\SMB_full_period.mat');

SMB_full_period_mp = SMB_full_period;

gla_id_smbs = [1319847, 1319824, 1318354];
gla_names_smb = {'Kyzylsu','Koshkul','Muzgazy'};

y_lim_smb = [4320200 4335333]; x_lim_smb =[707705 718617]; % SMB map limits
date_start_smb = datetime(2018,10,1);
date_end_smb = datetime(2019,9,1);

% Two periods for SMB profile comparison 
date_start_smb1 = datetime(1999,11,1);
date_end_smb1 = datetime(2018,9,1);

date_start_smb2 = datetime(2018,10,1);
date_end_smb2 = datetime(2023,9,1);

 run("Distributed_SMB_validation.m")

%% Compute mean annual SMB, pre-monsoon, monsoon, post-monsoon air temperature/precipitation

% gla_id_smb = 1319847;% Kyzylsu = 1319847, Koshkul = 1319824, Muzgazy = 1318354
% run("SMB_component_analysis.m")

%% Plot one snow cover map per month of MODIS vs T&C

if monthy_modis_figure == 1
    run("MODIS_spatial_comparison.m");
end

%% Monthly SWE maps

dir_fig_SWE = [dir_fig '\Monthly_SWE'];
if ~exist(dir_fig_SWE, 'dir')mkdir(dir_fig_SWE); end 

c4 = colormap(getPyPlot_cMap('GnBu'));

fi3 = figure('Renderer', 'painters', 'Position', [25.6667 87 733.3333 584.6667]);
for yy = 1:length(Years_no)
    
tiledlayout(3,4,'TileSpacing','compact');
month_of_that_year = unique(month(date_m(year(date_m) == Years_no(yy))));
for m = 1:numel(month_of_that_year)
id_m = find(month(date_m) == month_of_that_year(m) & year(date_m) == Years_no(yy),1);
if ~isempty(id_m)
nexttile
imagesc(x,y,SWEm_map(:,:,id_m),'AlphaData',DTM > 0); %clim([0 1]); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.9 0.9 0.9])
title(datestr(date_m(id_m),'dd-mmm-yyyy'))
colormap([1,1,1; c4])
if m == 4 || m == 8 || m == 12; cb = colorbar; set(cb,'Position',get(cb,'Position')+[-0.003 0 0 0]); 
ylabel(cb,'SWE [mm w.e.]','FontSize',11); end
set(gca,'YDir','normal'); clim([0 3000]);
if m == 9; xlabel(forcing_nm); end
end  
end
exportgraphics(fi3,[dir_fig_SWE '\' glacier '_monthly_SWE_' forcing_nm '_' num2str(Years_no(yy)) ... 
    '.png'],'Resolution',300,'BackgroundColor','none')
end 
close(gcf)

%% Avalanches contribution

gla_id_ava = 1319847;% Kyzylsu = 1319847, Koshkul = 1319824, Muzgazy = 1318354
run('Avalanches_quantification.m')

%% MODIS and L8/S2 fsc and SLA comparison 

%Needs
% - scas (daily fsc from T&C outputs)
% - date_h ( daily time vector)
% - MODIS and Landsat fsc/sla tables
% - Daily mean temperature and its corresponding time vector (SPAVG_dm.Date, SPAVG_dm.Ta_tg)
% - Daily sum of solid and total precipitation and time vectors(SPAVG_ds.Pr_sno_tg, SPAVG_ds.Pr_sno)  

Date_d = date_h-hours(12); % Corresponds to daily map timestamp

run("MODIS_LSS2_fsc_validation.m")

date_seas = date_m; 
SWEm_map = SWE_map;
run("Snow_cover_trends.m")

%% Map of final snow depth
c4 = colormap(getPyPlot_cMap('GnBu',10));
final_SND = flipud(snow_depth(:,:,end)); final_SND(isnan(flipud(DTM))) = NaN;

fi3 = figure('Renderer', 'painters', 'Position', [25.6667 87 733.3333 584.6667]);

imageschs(DEM,final_SND,'colormap',[1,1,1; c4]); %clim([0 1]); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.9 0.9 0.9])
title(datestr(date_m(end),'dd-mmm-yyyy'))
% colormap([1,1,1; c4]); 
cb = colorbar; %clim([0 50])
ylabel(cb,'Snow depth [m]','FontSize',11); 
set(gca,'YDir','normal');
exportgraphics(fi3,[dir_fig '\' glacier '_final_SND_' datestr(Date(end),'yyyy-mmm-dd') ...
    '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% High-resolution snow cover map validation

% Import landsat 5/7/8/9 and Sentinel-2 snow cover maps
run("Import_l8s2_snowcover.m")

%Analyse difference in snow cove between two periods, with the middle date defined below:

middle_date = datetime(2018,9,30);
run("Landsat_S2_validation.m")

%% Initial snow depth and comp with closest L8S2 image

%Needs:
% -Datestamp of the first timestep (startDate(1))
% -path to Landsat scenes (path_L8S2)
% -daily snow depth from T&C output (snow_depth)
% -DTM

run("Check_init_snow_depth.m")

%% Snowfall-Snowmelt-sublimation map

run("Snowfall_sublim_map.m")

% EGU abstract result

Snowfall_GLA = nanmean(Tot_snowfall(GLA_ID>0)./perYear);
Sublim_GLA = nanmean(Tot_sublim(GLA_ID>0)./perYear);

disp(['Annual on-glacier snowfall: ' num2str(round(Snowfall_GLA)) ' (mm)'])
disp(['Snowfall removed by sublimation: ' num2str(round(100*Sublim_GLA./Snowfall_GLA)) ' (%)'])

 %% Hydrograph figure
 
%  Needs:
% - the name of the outlet pixel (outlet_nm)
% - monthly or daily maps of liquid precip, snowmelt and icemelt

run("Hydrograph_figure.m"); 

%% Compute monthly-hourly temperature lapse-rates and check it corresponds to the ones provided to the downscaling

% Load air temperature LR used for downscaling
load(fn_LR)

for mm = 1:size(TA_map,3)
    T = TA_map(:,:,mm); T = T(:);
    E = flipud(DEM.Z); E(isnan(T)) = NaN; E = E(:);
    E(isnan(E)) = [];
    T(isnan(T)) = [];
    PLY = polyfit(E,T,1);
    LR_monthly(mm) = PLY(1);
end 

for m = 1:12
    LR_monthly_m(m) = nanmean(LR_monthly(month(date_m-days(15))==m));
end 

fi4 = figure('Renderer', 'painters', 'Position', [185 308.3333 578.6667 406.6667]);
plot(1:12, LR_monthly_m,'--r','LineWidth',1.1); grid on; hold on;
plot(1:12,nanmean(LR_hour_month,2),'k','LineWidth',1.1); grid on; hold on;
legend('LR in T&C','LR provided for downscaling','Location','NorthWest')
xlabel('Month of the year'); ylabel('Temperature lapse-rate [°C/m]')
title([glacier ' catchment'])
exportgraphics(fi4,[dir_fig '\T&C_TaLapseRates_check.png'],'Resolution',300,'BackgroundColor','none')

%% Figure 2. equivalent of Jouberton et al. 2022

% needs:
% - Monthly or daily maps of: air temperature, sold and total
% precipitation, avalanches, snowmelt, icemelt, evaporation (Ice and snow)
% - Monthly or daily catchment discharge
% - Corresponding time vector (date_m)

date_f2 = date_m;
hydro_year_f2_i = 1; % Use hydrological years (1, default), use calendar years (0)
show_ratio = 1; % Show monsoon snowfall ratio (1, default)

run("Fig2_Jouberton2022.m");

%% Figure 3 of Jouberton et al. 2022

% needs:
% - Monthly or daily maps of: air temperature, sold and total
% precipitation, avalanches, snowmelt, icemelt, evaporation (Ice and snow)
% - Corresponding time vector (date_m)

gla_fig3_num = 1; 
gla_fig3_id = gla_id_smbs(gla_fig3_num);
gla_fig3_label = gla_names_smb{gla_fig3_num};

date_f3 = date_m;

run("Fig3_Jouberton2022.m")

%% Catchment trend 

% needs:
% - Monthly or daily maps of: sublimation, precipitation, SWE, evaporation,
% snow melt, ice melt
% - Corresponding time vector (date_m)

trend_area = 'catchment'; % 'catchment','glaciers';

if length(Years_no) > 20
    date_d = Date_d;
    run("Trend_analysis.m")
end

%% Elevation profiles   

date_el = date_m(year(date_m)>1969);
volumes = 0;
middle_year = 2018;
run('Altitudinal_analysis.m')

%% Spatial analysis of runoff fluxes

date_start_flux_p1 = datetime(2000,1,1);
date_end_flux_p1 = datetime(2018,9,1);

date_start_flux_p2 = datetime(2018,10,1);
date_end_flux_p2 = datetime(2023,9,1);

run("Spatial_fluxes.m");

%% Seasonal cycle analysis

date_start_runoff_p1 = datetime(1999,10,1);
date_end_runoff_p1 = datetime(2012,9,1);

date_start_runoff_p2 = datetime(2012,10,1);
date_end_runoff_p2 = datetime(2023,9,1);

gla_id_smb = gla_id;
run('Seasonal_runoff.m')

% Heatmaps of anomalies
run('Heatmap_monthly_analysis.m')

%% ADDITIONAL FROM THAT POINT (site-specific)
%% Mean monthly variables 

for m = 1:12
    month_m = (month(date_m)-1); month_m(month_m==0)=12;
    ind_m = month_m == m;

    Psnow_m_i = nanmean(PSNOW_map(:,:,ind_m),3);
    Psnow_m(:,:,m) = Psnow_m_i;
    Psnow_m_GLA(m) = nanmean(Psnow_m_i(GLA_ID>0));

    Prain_m_i = nanmean(PRAIN_map(:,:,ind_m),3);
    Prain_m(:,:,m) = Prain_m_i;
    Prain_m_GLA(m) = nanmean(Prain_m_i(GLA_ID>0));

    % Compute mean monthly snow and rain per elevation bands

    dEL=100; % width of elevation bins
    ELs_GLA = nanmin(DTM(GLA_ID > 0),[],'all'):dEL:nanmax(DTM(GLA_ID > 0),[],'all');

    for iel = 1:numel(ELs_GLA)
        cur=(DTM<(ELs_GLA(iel)+dEL/2))&(DTM>=(ELs_GLA(iel)-dEL/2)); %current section of DEM
        pPrain_el_m(iel,m) = nanmean(Prain_m_i(cur & GLA_ID > 0)); 
        pPsnow_el_m(iel,m) = nanmean(Psnow_m_i(cur & GLA_ID > 0)); 
    end 
end 

% fi4 = figure('Renderer', 'painters', 'Position', [185 127.6667 724 587.3333]);
% tiledlayout(4,3,'TileSpacing','compact')
% for m = 1:12
% nexttile
% bar(ELs_GLA, [pPrain_el_m(:,m) pPsnow_el_m(:,m)],'stacked');
% if m == 1; legend('Rain','Snowfall'); end 
% view(90,-90); grid on; ylim([0 300]);
% if ~ismember(m, 10:12); set(gca,'YTickLabels',[]); else ylabel('Precip [mm]','FontSize',11); end
% if ~ismember(m, [1 4 7 10]); set(gca,'XTickLabels',[]); else xlabel('Elevation [m a.s.l.]','FontSize',11); end
% text(3500, 230, month_labels{m},'FontWeight','bold')
% end 

nansum(Psnow_m_GLA(3:6)./nansum(Psnow_m_GLA))
nansum(Psnow_m_GLA([12 1 2 3])./nansum(Psnow_m_GLA))
nansum(Psnow_m_GLA([6 7 8 9])./nansum(Psnow_m_GLA))

%% Compare Pleiades snow depth with T&C snow depth

Pleiades_path = 'B:\group\pelligrp\Pleiades_stereo_images\Differencing\Kyzylsu\Coreg_Kyzylsu_2022_09_24-2023_05_17_2m_32642';
Pleiades_dh = GRIDobj([Pleiades_path '\Kyzylsu_2022_09_24-2023_05_17_diff_NK-DR_de-unduling_10m.tif']);
Pleiades_dh_100m = resample(Pleiades_dh,DEM,'bilinear');
Pleiades_ortho = GRIDobj([Pleiades_path '\Kyzylsu_2023_05_17_coregistered_truecolor_10m.tif']);
Pleiades_info = georasterinfo([Pleiades_path '\Kyzylsu_2023_05_17_coregistered_truecolor_10m.tif']);
Pleiades_dem = GRIDobj([Pleiades_path '\Kyzylsu_2022_09_24-2023_05_17_dem_Marin_de-unduling_10m.tif']);
Pleiades_rgb = double(cat(3,Pleiades_ortho.Z(:,:,1),Pleiades_ortho.Z(:,:,2),Pleiades_ortho.Z(:,:,3))./4000);
[x_ortho, y_ortho] = worldGrid(Pleiades_info.RasterReference);
Pleiades_ortho.Z = Pleiades_rgb;
res_pleiades = Pleiades_ortho.cellsize;

Pleiades_slope = gradient8(Pleiades_dem,'degree');

Pleiades_dh_filt = Pleiades_dh.Z;
Pleiades_dh_filt(Pleiades_slope.Z > 30) = NaN;

%Load AW3D dem
DEM_aw3d_orig = GRIDobj('B:\group\pelligrp\Remote_Sensing_Data\Regional_Global\AW3D\N039E071\N039E071_AVE_DSM.tif');
DEM_aw3d = reproject2utm(DEM_aw3d_orig,30);
DEM_aw3d = crop(DEM_aw3d, [demLons(1,1) demLons(end,end)], [demLats(1,1) demLats(end,end)]);
res_aval = 100;
DEM_aw3d_res = resample(DEM_aw3d,res_aval);
DEM_aval = DEM_aw3d_res; DEM_aval.Z = double(DEM_aval.Z);
cellsize = DEM_aval.cellsize;
MASK_aval = ~isnan(DEM_aval.Z);

% Load avalanche deposits manually delineated from Pleiades differencing
Aval_deposit = shaperead([root '\Study_sites\Kyzylsu\Pleiades_differencing\Avalanche_deposits_Kyzylsu_2023_05_17_utm.shp']);
AVAL_deposit = polygon2GRIDobj(DEM_aval,Aval_deposit);
AVAL_deposit_pleaides = polygon2GRIDobj(Pleiades_dem,Aval_deposit);

% Extract T&C snow depth map at the Pleiades acquisition date

Date_acq_start = datetime(2022,9,24);
Date_acq = datetime(2023,5,17);
Date_acq_tc_id = find(date_d==Date_acq,1); Date_acq_tc_id_mm = find(date_m==Date_acq,1); 
SND_pleiades_tc = DEM;
SND_pleiades_tc.Z = flipud(snow_depth(:,:,Date_acq_tc_id));

% Compute avalanched SWE from T&C corresponding to the acquisition dates
diff_start = abs(date_m - Date_acq_start);
[~, Date_acq_tc_id_start_mm] = min(diff_start);
diff_end = abs(date_m - Date_acq);
[~, Date_acq_tc_id_mm] = min(diff_end);
AVAL_pleiades_tc = DEM;
AVAL_pleiades_tc.Z = flipud(nansum(AVA_map(:,:,Date_acq_tc_id_start_mm:Date_acq_tc_id_mm),3));

% Compute slope and aspect of AW3D resampled DEM
[Slo_top,Aspect]=Slope_Aspect_indexes(DEM_aval.Z,cellsize,'mste');
Asur = (1./cos(atan(Slo_top)));
Area = (cellsize^2)*sum(sum(MASK_aval)); %% Projected area [m^2]

% Load deposit detected from Sentinel-1 (Kneib et al. 2024)
S1_ASC = resample(GRIDobj('C:\Users\jouberto\Desktop\Study_sites\Kyzylsu\Spatial_data\Avalanches_S1\Heat_map_ASC_022017-122023.tif'), DEM);
S1_DESC = resample(GRIDobj('C:\Users\jouberto\Desktop\Study_sites\Kyzylsu\Spatial_data\Avalanches_S1\Heat_map_DESC_022017-122023.tif'), DEM);
S1_ALL = S1_ASC; S1_ALL.Z = S1_ASC.Z + S1_DESC.Z;

S1_BRIGHT_ASC = resample(GRIDobj('C:\Users\jouberto\Desktop\Study_sites\Kyzylsu\Spatial_data\Avalanches_S1\Brightness_ASC.tif'), DEM);
S1_BRIGHT_DESC = resample(GRIDobj('C:\Users\jouberto\Desktop\Study_sites\Kyzylsu\Spatial_data\Avalanches_S1\Brightness_DESC.tif'), DEM);

S1_mask = S1_BRIGHT_ASC.Z<0.82 | S1_BRIGHT_DESC.Z < 0.82 & (atan(Slo_top).*180/pi < 35);

%% Pleiades avalanche deposit - catchment validation

xmin_aval = [707705 7.115*10^5 7.09*10^5];%xmin_map+2000;
xmax_aval = [718617 7.143*10^5 7.123*10^5]; %xmax_map-2000
ymin_aval = [4322500 4.328*10^6 4.3253*10^6];
ymax_aval = [4333333 4.3308*10^6 4.32857*10^6];
labels_zoom = {'Catchment','Koshkul','Kyzylsu'};
depo_line_width = 0.6;

cmap_redblue = flipud(redblue); cmap_redblue = cmap_redblue(1:255,:);

fi3 =figure('Renderer', 'painters', 'Position', [104.3333 55.6667 704.6667 636.6667]);
tiledlayout(3,3,"TileSpacing","compact")

for zz = 1:length(labels_zoom)

n1 = nexttile;
imagesc(Pleiades_ortho); hold on;
set(gca,'YDir','normal')
plot(catchShp_utm.X, catchShp_utm.Y,'k','LineWidth',0.8); 
if zz == 1
    plot([xmin_aval(2) xmax_aval(2) xmax_aval(2) xmin_aval(2) xmin_aval(2)], [ymin_aval(2) ymin_aval(2) ymax_aval(2) ymax_aval(2) ymin_aval(2)],'r','LineWidth',0.9)
    plot([xmin_aval(3) xmax_aval(3) xmax_aval(3) xmin_aval(3) xmin_aval(3)], [ymin_aval(3) ymin_aval(3) ymax_aval(3) ymax_aval(3) ymin_aval(3)],'b','LineWidth',0.9)
end
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.9290 0.6940 0.1250]); end 
for yy = 1:size(Aval_deposit,1); plot(Aval_deposit(yy).X,Aval_deposit(yy).Y,'Color',[255, 105, 180]./255,'LineWidth',depo_line_width); end
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); 
xlim([xmin_aval(zz) xmax_aval(zz)]); ylim([ymin_aval(zz) ymax_aval(zz)]);
if zz == 1; title('Pleiades ortho-image','FontSize',9); end 
if zz == 2; set(n1,'YColor',[1 0 0]); set(n1,'XColor',[1 0 0]); end
if zz == 3; set(n1,'YColor',[0 0 1]); set(n1,'XColor',[0 0 1]); xlabel(datestr(Date_acq),'Color','k'); end 

nexttile
imageschs(Pleiades_dem,Pleiades_dh,'caxis',[-10 10],'colormap',cmap_redblue,'colorbar',0); hold on; %clim([-10 10]); colormap(flipud(redblue)); hold on
plot(catchShp_utm.X, catchShp_utm.Y,'k','LineWidth',0.8);
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); 
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.9290 0.6940 0.1250]); end 
for yy = 1:size(Aval_deposit,1); plot(Aval_deposit(yy).X,Aval_deposit(yy).Y,'Color',[255, 105, 180]./255,'LineWidth',depo_line_width); end
% cb1 = colorbar('westoutside'); ylabel(cb1, 'Elevation change [m]')
 xlim([xmin_aval(zz) xmax_aval(zz)]); ylim([ymin_aval(zz) ymax_aval(zz)]);
if zz == 1; title('Pleiades elevation change [m]','FontSize',9); end
if zz == 3; xlabel([datestr(Date_acq) ' - ' datestr(Date_acq_start)]); end 

n3 = nexttile;
imageschs(DEM,SND_pleiades_tc,'caxis',[-10 10],'colormap',cmap_redblue,'colorbar',1,'brighten',0)%,'AlphaData',flipud(MASK))
set(gca,'Color',[0.8 0.8 0.8])
%clim([-10 10]); colormap(flipud(redblue)); 
hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.9290 0.6940 0.1250]); end 
for yy = 1:size(Aval_deposit,1); plot(Aval_deposit(yy).X,Aval_deposit(yy).Y,'Color',[255, 105, 180]./255,'LineWidth',depo_line_width); end
 xlim([xmin_aval(zz) xmax_aval(zz)]); ylim([ymin_aval(zz) ymax_aval(zz)]);
 plot(catchShp_utm.X, catchShp_utm.Y,'k','LineWidth',0.8)
% cb1 = colorbar('westoutside'); ylabel(cb1, 'Snow depth [m]','FontSize',11)
if zz ==1; title('Simulated snow depth [m]','FontSize',9); end 
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
if zz == 3; xlabel(datestr(Date_acq)); end 
end
exportgraphics(fi3,[dir_fig '\Avalanching\Pleiades_aval_TC_evaluation.png'],'Resolution',300,'BackgroundColor','none')

%% Compare Sentinel-1 avalanche deposits with simulated deposits:

xmin_aval = 707705;
xmax_aval = 718617; 
ymin_aval = 4321000;
ymax_aval = 4333333;
SWE_thres_deposit = 0;
period_S1_aval = ismember(year(date_m), 2017:2023);
mask_perf = S1_mask==1 & flipud(MASK==1); % Mask over which compare Sentinel-1 and T&C outputs

SWE_avalanched_i = flipud(nansum(AVA_map(:,:,period_S1_aval),3));
SND_show = DEM;
SND_show.Z = DEM.Z.*NaN;
SND_show.Z((SWE_avalanched_i >SWE_thres_deposit) & (S1_ALL.Z>0 & S1_mask==1)) = 1;
SND_show.Z((SWE_avalanched_i >SWE_thres_deposit) & (S1_ALL.Z==0 & S1_mask==1)) = 2;
SND_show.Z((SWE_avalanched_i <=SWE_thres_deposit) & (S1_ALL.Z>0 & S1_mask==1)) = 3;

nPix_TruPos_S1_SnowSlide = nansum(SWE_avalanched_i(mask_perf)>SWE_thres_deposit & S1_ALL.Z(mask_perf)>0,'all');
nPix_TruNeg_S1_SnowSlide = nansum(SWE_avalanched_i(mask_perf)<SWE_thres_deposit & S1_ALL.Z(mask_perf) ==0,'all');
nPix_FalPos_S1_SnowSlide = nansum(SWE_avalanched_i(mask_perf)>SWE_thres_deposit & S1_ALL.Z(mask_perf) ==0,'all');
nPix_FalNeg_S1_SnowSlide = nansum(SWE_avalanched_i(mask_perf)<=SWE_thres_deposit & S1_ALL.Z(mask_perf) >0,'all');
nPix_catch_affected_SnowSlide = nansum(SWE_avalanched_i(mask_perf)>SWE_thres_deposit)./nansum(mask_perf,'all');
S1_catch_desposit_affected = nansum(S1_ALL.Z(mask_perf) >0,'all')./nansum(mask_perf,'all');
Mean_catchment_deposit_error = nPix_catch_affected_SnowSlide - S1_catch_desposit_affected;
F1_score = 2.*nPix_TruPos_S1_SnowSlide./(2.*nPix_TruPos_S1_SnowSlide + nPix_FalPos_S1_SnowSlide + nPix_FalNeg_S1_SnowSlide);

fi3 =figure('Renderer', 'painters', 'Position', [179 97.6667 637.3333 596]);
imageschs(DEM,SND_show.Z,'colormap',flipud(cbrewer('qual','Set1',3)));%, 'AlphaData', flipud(MASK) & S1_mask)
cb = colorbar; set(cb,'XTick',[1.33 2 2.66]); set(cb,'XTickLabel',{'True Positive','False positive','False negative'}); 
%clim([-2 2]); 
set(gca,'ZTick',[0.66 1.33 2 ],'CLim',[1 3])
hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k')
% cb1 = colorbar('westoutside'); ylabel(cb1,'Snow displaced [m w.e.]','FontSize',11); 
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.9290 0.6940 0.1250]); end 
for yy = 1:size(Aval_deposit,1); plot(Aval_deposit(yy).X,Aval_deposit(yy).Y,'Color',[255, 105, 180]./255,'LineWidth',0.4); end
set(gca,'YDir','normal'); set(gca,'Color',[0.8 0.8 0.8])
% set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
% xlabel(['V_{deposit} = ' num2str(Aval_deposit_SnowSlide_volume(id).*10^-6,3) '* 10^6 m^3'])
xlim([xmin_aval xmax_aval]); ylim([ymin_aval ymax_aval]); 
xlabel(['Deposit cover fraction: ' num2str(S1_catch_desposit_affected,2) ', MeanErr: ' num2str(Mean_catchment_deposit_error,2) ', Fscore: ' num2str(F1_score,2)],'FontSize',10)
title('Sentinel-1 avalanche deposit comparison (2017-2023)')
exportgraphics(fi3,[dir_fig '\Sentinel_aval_TC_validation.png'],'Resolution',300,'BackgroundColor','none')

%%%% Avalanches desposits comparison with Sentinel-1 data  (AVAL)

SWE_avalanched_show = SWE_avalanched_i;
SWE_avalanched_show(flipud(MASK) ~= 1) = NaN;
id = 1;
%%
fi3 =figure('Renderer', 'painters', 'Position', [123 281.6667 520.6667 430]);
t = tiledlayout(1,1,"TileSpacing","compact");
ax1 = axes(t);
SND_show = DEM;
SND_show.Z = SWE_avalanched_show./(nansum(period_S1_aval)./12);
imagesc(SND_show, 'AlphaData', flipud(MASK))% & S1_mask)
ax1.XTick = []; ax1.XTickLabels = [];
ax1.YTick = []; ax1.YTickLabels = [];
ax1.Color = [0.9 0.9 0.9];
% ax1.Layout.Tile=1; % <--- key new line of code
clim([-1000 1000]); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k')
if id == 1; cb1 = colorbar('westoutside'); ylabel(cb1,'Snow displaced [mm w.e./yr]','FontSize',12); end
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.8 0.8 0.8]); end 
for yy = 1:size(Aval_deposit,1); plot(Aval_deposit(yy).X,Aval_deposit(yy).Y,'Color',[255, 105, 180]./255,'LineWidth',0.4); end
ax2 = axes(t);
im2 = imagesc(S1_ALL,'AlphaData',(S1_ALL.Z>0 & S1_mask).*0.6);
linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = []; ax2.XTickLabels = [];
ax2.YTick = []; ax2.YTickLabels = [];
% ax2.Layout.Tile=1; % <--- key new line of code
colormap(ax1,flipud(redblue));
colormap(ax2,'copper');
set(ax2,'color','none','visible','off');
cb1 = colorbar('eastoutside'); ylabel(cb1,'Deposits occurence [-]','FontSize',12); 
set(gca,'YDir','normal'); set(gca,'Color',[0.8 0.8 0.8])
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
% xlabel(['V_{deposit} = ' num2str(Aval_deposit_SnowSlide_volume(id).*10^-6,3) '* 10^6 m^3'])
% xlim([xmin_aval xmax_aval]); ylim([ymin_aval ymax_aval]); 
axis image
exportgraphics(fi3,[dir_fig '\Avalanching\Sentinel_aval_TC_overlay.png'],'Resolution',300,'BackgroundColor','none')

%%
% % Solve monthly maps which can be used for SnowSlide calibration
% 
% SS_input_maps.SWE = SWEm_map;
% SS_input_maps.SMS = SMSm_map;
% SS_input_maps.SSN = SSN_map;
% SS_input_maps.SND = SNDm_map;
% SS_input_maps.AVA_map = AVA_map;
% SS_input_maps.PSNOW = PSNOW_map;
% SS_input_maps.date = date_m;
% SS_input_maps.Asur = Asur;
% SS_input_maps.mask = MASK;
% SS_input_maps.DEM = DEM;
% SS_input_maps.Slo_top = Slo_top;
% 
% save([dir_fig '/SS_input_maps'],'SS_input_maps')