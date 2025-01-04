%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPATIAL PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CLEANING
clc; clear all;

%% DIRECTORIES WHERE RESULTS ARE STORED
path_output = 'M:/19_ISTA/1_TC/3_Model_Source/2_MegaWat/2_Pyrenees_distributed/5_Output_files/1_Model_outputs/';
result = 'Region_(Distributed_Ebro)_SubBasin_(Cinca_Mid)_Forcing_(Chelsa)_RunningDate_(031224_1301)_RunningPeriod_(2021_2021)/';
Directory=[path_output result];
cd(Directory)

load("Spatial_data/Osat_1.mat");
load("Spatial_data/Ohy_1.mat");
SNOW349 = load("Spatial_data_intermediate/OUTPUT_Cinca_Mid_SNOWMAP_349.mat");
SNOW445 = load("Spatial_data_intermediate/OUTPUT_Cinca_Mid_SNOWMAP_445.mat");

%% Base map
dtm_shape =flipud(SNOW349.DTM);
dtm_shape(dtm_shape >= 0) = 1;

%% DEM
figure(1)
dtm_flipped =flipud(SNOW349.DTM)
imagesc(dtm_flipped);
caxis([min(dtm_flipped(:)) max(dtm_flipped(:))]);
cmap = [1 1 1; copper]; % White for NaN, gray for other values
colormap(cmap);
colorbar;


%% Spatial SWE
figure(2)
SWE_MAT = reshape(Smelt_spatial_daily, length(x_cell), length(y_cell));
imagesc(SWE_MAT)

%Spatial SWE
figure(3)
SSN_MAT = reshape(SSN_spatial_daily, length(x_cell), length(y_cell));
imagesc(SSN_MAT)

%Osat
figure(4)
imagesc(Osat)

%Osat
figure(5)
imagesc(Ohy)

%Matrix
figure(6)
[M, c]=contour(DTM, [1000 1500 2000 2500 3000 3500 4000], "ShowText",true, ...
    "LabelFormat","%d m");
c.LineWidth = 2;

%Contour and heatmap overlapped
figure(7)
Ohy_flipped = flipud(Ohy);
imagesc(Ohy_flipped);
caxis([min(Ohy_flipped(:)) max(Ohy_flipped(:))]);
colormap(sky);
colorbar;
set(gca, 'XTick', [], 'YTick', []);

hold on;
DTM_flipped = flipud(DTM);
[M2, c2]=contour(DTM_flipped, [1000 1500 2000 2500 3000 3500 4000 5000], "ShowText",true, ...
    "LabelFormat","%d m", 'LineColor', 'k');
hold off;

%Labels
xlabel({'EAST'});
ylabel({'NORTH'});

%xl = xlim; yl = ylim;

%X ticks
xcord_x = 1:length(x_cell);
xcord_y = 240*ones(1, length(x_cell));
xcord_text = strsplit(num2str(round(x_cell,0)));
xc_selected= 1:20:length(x_cell);

xticks(xc_selected);
xticklabels(xcord_text(xc_selected));

%Y ticks
new_y_cell = flip(y_cell);
ycord_x = 1:length(new_y_cell);
ycord_y = 240*ones(1, length(new_y_cell));
ycord_text = strsplit(num2str(round(new_y_cell,0)));
yc_selected= 1:20:length(new_y_cell);

yticks(yc_selected);
yticklabels(ycord_text(yc_selected));
%{
text(xcord_x(xc_selected), xcord_y(xc_selected),xcord_text(xc_selected) , ...
    'Rotation',45, 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
%}

%% ICE IN THE CATCHMENT
figure(8)

imagesc(dtm_shape);
myColorMap = [0 0 0; 1 1 1]; % two colors
colormap(myColorMap);

hold all
A_ICE_resh = reshape(A.Ice_spatial_daily,length(A.y_cell),length(A.x_cell));
A_ICE_resh_flipped = flipud(A_ICE_resh);

% Create an alpha matrix
A_ICE_resh_flipped_trans = zeros(size(A_ICE_resh_flipped)); % Start with all opaque
%A_ICE_resh_flipped_trans(A_ICE_resh_flipped == 0) = 0; % Set transparent for zeros

h2 = imagesc(A_ICE_resh_flipped);
h2.AlphaData = 0.5*A_ICE_resh_flipped;
set(h2, 'AlphaData', A_ICE_resh_flipped_trans); % Adjust transparency as needed
myColorMap = [0 0 0, 'jet']
colormap(myColorMap); % Set colormap if needed