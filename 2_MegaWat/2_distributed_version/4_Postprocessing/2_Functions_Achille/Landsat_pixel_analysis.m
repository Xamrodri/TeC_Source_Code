function [POI_obs, POI_obs_mean, POI_obs_std, Date_l8s2] = Landsat_pixel_analysis(Poi_name,var_obs,POI_table,...
    path_L8S2, glacier, Date)

% Poi_name = 'Pluvio';  % POI at which to extract the Landsat/Sentinel-2 time-series
% var_obs = 'Albedo';    % Which observed variable to extract from Landsat/Sentinel

Xutm_tc = POI_table.LON_UTM(strcmp(Poi_name,POI_table.Name));
Yutm_tc = POI_table.LAT_UTM(strcmp(Poi_name,POI_table.Name));

fn_maps = dir([path_L8S2 '/TIF/*' var_obs '*.tif']);

% Extract the dates for which landsat/sentinel-2 snow cover maps are available
for ii = 1:size(fn_maps,1)
    fn_scene = fn_maps(ii).name;
    date_l8s2_char = char(extractBetween(fn_scene,[glacier '_'],'.tif'));
    sat_l8s2_char = extractBefore(fn_scene,'_Albedo');
    sat_l8s2{ii} = sat_l8s2_char;
    date_l8s2(ii) = datetime(str2num(date_l8s2_char(1:4)),str2num(date_l8s2_char(6:7)),str2num(date_l8s2_char(9:10)));
end

% Extract the indices for which landsat/sentinel-2 maps we have T&C outputs
for ii = 1:numel(Date)
    find_id_L8S2 = find(datetime(year(Date(ii)),month(Date(ii)),day(Date(ii))) == date_l8s2',1);

    if ~isempty(find_id_L8S2)
        id_l8s2(ii,1) = find_id_L8S2;
        id_l8s2(ii,2) = ii;
    else
        id_l8s2(ii,1) = 0;
        id_l8s2(ii,2) = ii;
    end 
end
id_l8s2(id_l8s2(:,1) == 0,:) = [];

% Load geospatial information from the first scene
[~, SCENE1_info] = readgeoraster([path_L8S2 '\TIF\' fn_maps(1).name]);
SCENE1_INFO = georasterinfo([path_L8S2 '\TIF\' fn_maps(1).name]);
[x_l8s2,y_l8s2] = worldGrid(SCENE1_info);
% hls_GRIDobj = GRIDobj([path_L8S2 '\TIF\' fn_maps(1).name]);
% DEM_l8s2 = project(DEM,hls_GRIDobj);

% Find the coordinate of the Landsat images which is closest to our POI our interest
% [l8s2_obs] = NaN(SCENE1_INFO.RasterSize(1),SCENE1_INFO.RasterSize(2), size(id_l8s2,1));
DIST = sqrt((x_l8s2-Xutm_tc).^2 + (y_l8s2-Yutm_tc).^2);

[xi_l8s2, yi_l8s2] = find(DIST == min(DIST,[],'all'));

[POI_obs, POI_obs_mean, POI_obs_std] = deal(NaN(size(id_l8s2,1),1));

for ii = 1:size(id_l8s2,1)
   l8s2_obs_load =  imread([path_L8S2 '\TIF\' fn_maps(id_l8s2(ii,1)).name]);
   POI_obs(ii) =  l8s2_obs_load(xi_l8s2,yi_l8s2);
   POI_obs_mean(ii) =  nanmean(l8s2_obs_load((xi_l8s2-1):1:(xi_l8s2+1),(yi_l8s2-1):1:(yi_l8s2+1)),'all');
   POI_obs_std(ii) =  nanstd(l8s2_obs_load((xi_l8s2-1):1:(xi_l8s2+1),(yi_l8s2-1):1:(yi_l8s2+1)),[],'all');
end       

Date_l8s2 = date_l8s2(id_l8s2(:,1));

% figure; imagesc(x_l8s2(1,:), y_l8s2(:,1), l8s2_obs_load); hold on; %scatter(yi_l8s2,xi_l8s2,'r+')
% plot(catchShp_utm.X, catchShp_utm.Y,'r')
% scatter(x_l8s2(1,yi_l8s2),y_l8s2(xi_l8s2,1),'r+')
% set(gca,'YDir','normal')
% 
figure
scatter(date_l8s2(id_l8s2(:,1)),POI_obs); grid on; hold on;
scatter(date_l8s2(id_l8s2(:,1)),POI_obs_mean,'r','filled')
plot(Pluvio.Date(hour(Pluvio.Date)==12), Pluvio.SWout(hour(Pluvio.Date)==12)./Pluvio.SWin(hour(Pluvio.Date)==12),'Color',[0 0 0 0.6],'LineWidth',1.1)

% % plot(AWS.AWS5850_DATE(ind_noon_meas), AWS.Albedo(ind_noon_meas),'LineWidth',1.1);