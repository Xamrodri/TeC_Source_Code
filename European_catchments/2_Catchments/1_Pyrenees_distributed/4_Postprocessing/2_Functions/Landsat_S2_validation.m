

%% Maps of snow cover frequency per season

ind_win = month(Date_l8s2) > 11 | month(Date_l8s2) < 3;
ind_spr = month(Date_l8s2) > 2 & month(Date_l8s2) < 6;
ind_sum = month(Date_l8s2) > 5 & month(Date_l8s2) < 9;
ind_aut = month(Date_l8s2) > 8 & month(Date_l8s2) < 12;

Snow_freq_win = nansum(l8s2_sno(:,:,ind_win),3)./nansum(~isnan(l8s2_sno(:,:,ind_win)),3);
Snow_freq_spr = nansum(l8s2_sno(:,:,ind_spr),3)./nansum(~isnan(l8s2_sno(:,:,ind_spr)),3);
Snow_freq_sum = nansum(l8s2_sno(:,:,ind_sum),3)./nansum(~isnan(l8s2_sno(:,:,ind_sum)),3);
Snow_freq_aut = nansum(l8s2_sno(:,:,ind_aut),3)./nansum(~isnan(l8s2_sno(:,:,ind_aut)),3);

ind_seas = [ind_aut, ind_win, ind_spr, ind_sum];
season_labels = {'Sep - Nov','Dec - Feb','Mar - May','Jun - Aug'};

cmap_snow_freq = cbrewer('seq','Blues',80); cmap_snow_freq(cmap_snow_freq<0) = 0; cmap_snow_freq(cmap_snow_freq>1) = 1;
cmap_snow_std = cbrewer('seq','Reds',80); cmap_snow_freq(cmap_snow_freq<0) = 0; cmap_snow_freq(cmap_snow_freq>1) = 1;
%
fi2 = figure('Renderer', 'painters', 'Position',[203 88.3333 497.3333 562.6667]) ;
tiledlayout(2,2,'TileSpacing','tight')
for ii = 1:size(ind_seas,2)
n1 = nexttile;
imagesc(demLons(1,:), demLats(:,1),nanmean(l8s2_sno(:,:,ind_seas(:,ii)),3),'AlphaData',MASK_hls & ~isnan(nanmean(l8s2_sno(:,:,ind_seas(:,ii)),3))); 
title([season_labels{ii} ' (' num2str(sum(ind_seas(:,ii))) ': scenes)']); set(gca,'YDir','normal')
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.8 0.8 0.8]) 
clim([0 1]); axis image; colormap(n1,cmap_snow_freq); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons(:)) max(demLons(:))]); ylim([min(demLats(:)) max(demLats(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
if ismember(ii,[2 4]); cb1 = colorbar; ylabel(cb1,'Mean snow presence [-]','FontSize',11); end
end
exportgraphics(fi2,[dir_fig '\' glacier '_all4season_SnowFreqMaps.png'],'Resolution',300,'BackgroundColor','none')

%% Compute snow cover frequency per month and per year  [TO DELETE IF NO PROBLEM OCCURS, commented on 03-11-2024

% Years_l8s2 = unique(year(Date_l8s2));
% L8s2_ym_sp = deal(nan(size(DEM_l8s2_orig.Z,1),size(DEM_l8s2_orig.Z,2),numel(Years_l8s2).*12));
% 
% for yy = 1:numel(Years_l8s2)
% for mm = 1:12
% 
% ind_scene = month(Date_l8s2) == mm & year(Date_l8s2) == Years_l8s2(yy);
% 
% if nansum(ind_scene) == 0 % skip in case there is one month without any l8s2 data
%     continue; 
% end 
% 
% Snow_occurence_m = nansum(l8s2_sno(:,:,ind_scene),3);
% Data_occurence_m = nansum(~isnan(l8s2_sno(:,:,ind_scene)),3);
% 
% Data_occurence_m(Data_occurence_m==0) = NaN; % Put as NaN pixels for which there are no data for the whole given month
% Snow_occurence_m(Data_occurence_m==0) = NaN; % Put as NaN pixels for which there are no data for the whole given month
% 
% L8s2_ym_sp(:,:,12*(yy-1)+mm) = Snow_occurence_m./Data_occurence_m;
% end 
% end 
% 
% date_l8s2_m = (datetime(Years_l8s2(1),1,15,0,0,0):calmonths(1):datetime(Years_l8s2(end),12,15,0,0,0))'; 
% date_l8s2_m_data =  squeeze(nansum(L8s2_ym_sp,[1 2]) > 0);
% L8s2_ym_sp = L8s2_ym_sp(:,:,date_l8s2_m_data);
% date_l8s2_m = date_l8s2_m(date_l8s2_m_data); clear date_l8s2_m_data

%% Maps of snow cover frequency per season between two periods

% Compute monthly mean snow frequency

month_labels = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};

for ii = 1:12

ind_mm_p1 = (month(Date_l8s2) == ii)  & Date_l8s2 < middle_date;
ind_mm_p2 = (month(Date_l8s2) == ii)  & Date_l8s2 >= middle_date;
    
ind_nan_p1 = sum(isnan(l8s2_sno(:,:,ind_mm_p1)),3)<(0.33*size(l8s2_sno(:,:,ind_mm_p1),3));    
ind_nan_p2 = sum(isnan(l8s2_sno(:,:,ind_mm_p2)),3)<(0.33*size(l8s2_sno(:,:,ind_mm_p2),3));    

Monthly_freq_map_p1_i = nanmean(l8s2_sno(:,:,ind_mm_p1),3);
Monthly_freq_map_p1_i(~MASK_hls | ~ind_nan_p1) = NaN;
Monthly_freq_map_p1(:,:,ii) = Monthly_freq_map_p1_i;

Monthly_freq_map_p2_i = nanmean(l8s2_sno(:,:,ind_mm_p2),3);
Monthly_freq_map_p2_i(~MASK_hls | ~ind_nan_p2) = NaN;
Monthly_freq_map_p2(:,:,ii) = Monthly_freq_map_p2_i;

end

ind_seas_p1 = [ismember(month(Date_l8s2),10:12)  & Date_l8s2 < middle_date ... 
    ismember(month(Date_l8s2),1:3)  & Date_l8s2 < middle_date  ...
    ismember(month(Date_l8s2),4:6)  & Date_l8s2 < middle_date ...
    ismember(month(Date_l8s2),7:9)  & Date_l8s2 < middle_date];

ind_seas_p2 = [ismember(month(Date_l8s2),10:12)  & Date_l8s2 >= middle_date ... 
    ismember(month(Date_l8s2),1:3)  & Date_l8s2 >= middle_date  ...
    ismember(month(Date_l8s2),4:6)  & Date_l8s2 >= middle_date ...
    ismember(month(Date_l8s2),7:9)  & Date_l8s2 >= middle_date];

Season_freq_map_p1(:,:,1) = nanmean(Monthly_freq_map_p1(:,:,10:12),3);
Season_freq_map_p1(:,:,2) = nanmean(Monthly_freq_map_p1(:,:,1:3),3);
Season_freq_map_p1(:,:,3) = nanmean(Monthly_freq_map_p1(:,:,4:6),3);
Season_freq_map_p1(:,:,4) = nanmean(Monthly_freq_map_p1(:,:,7:9),3);

Season_freq_map_p2(:,:,1) = nanmean(Monthly_freq_map_p2(:,:,10:12),3);
Season_freq_map_p2(:,:,2) = nanmean(Monthly_freq_map_p2(:,:,1:3),3);
Season_freq_map_p2(:,:,3) = nanmean(Monthly_freq_map_p2(:,:,4:6),3);
Season_freq_map_p2(:,:,4) = nanmean(Monthly_freq_map_p2(:,:,7:9),3);

season_labels = {'Oct - Dec','Jan - Mar','Apr - Jun','Jul - Sep'};

cmap_snow_freq = cbrewer('seq','Blues',80); cmap_snow_freq(cmap_snow_freq<0) = 0; cmap_snow_freq(cmap_snow_freq>1) = 1;
cmap_snow_std = cbrewer('seq','Reds',80); cmap_snow_freq(cmap_snow_freq<0) = 0; cmap_snow_freq(cmap_snow_freq>1) = 1;

%% Analyze seasonal changes in Sentinel-2/ Landsat_8 snow cover

for ii = 1:size(ind_seas_p1,2)

fi2 = figure('Renderer', 'painters', 'Position',[89.6667 275.6667 812.0000 377.9999]) ;
tiledlayout(1,3,'TileSpacing','tight')
n1 = nexttile;
imagesc(demLons_l8s2(1,:), demLats_l8s2(:,1),Season_freq_map_p1(:,:,ii)); 
title([season_labels{ii} ' (' num2str(sum(ind_seas_p1(:,ii))) ': scenes)']); set(gca,'YDir','normal')
%title([season_labels{ii}]); set(gca,'YDir','normal')
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.8 0.8 0.8]) 
clim([0 1]); axis image; cb1 = colorbar('Location','westoutside'); ylabel(cb1,'Mean snow presence [-]','FontSize',11); colormap(n1,cmap_snow_freq); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons(:)) max(demLons(:))]); ylim([min(demLats(:)) max(demLats(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
xlabel([num2str(min(year(Date_l8s2(ind_seas_p1(:,ii))))) '-'  num2str(max(year(Date_l8s2(ind_seas_p1(:,ii)))))])
n2 = nexttile;
imagesc(demLons_l8s2(1,:), demLats_l8s2(:,1),Season_freq_map_p2(:,:,ii)); 
title([season_labels{ii} ' (' num2str(sum(ind_seas_p2(:,ii))) ': scenes)']); set(gca,'YDir','normal')
%title([season_labels{ii}]); set(gca,'YDir','normal')
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.8 0.8 0.8]) 
clim([0 1]); axis image; %cb1 = colorbar; ylabel(cb1,'Mean snow presence [-]','FontSize',11); 
colormap(n2,cmap_snow_freq); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons(:)) max(demLons(:))]); ylim([min(demLats(:)) max(demLats(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
xlabel([num2str(min(year(Date_l8s2(ind_seas_p2(:,ii))))) '-'  num2str(max(year(Date_l8s2(ind_seas_p2(:,ii)))))])

n3 = nexttile;
imagesc(demLons_l8s2(1,:), demLats_l8s2(:,1),Season_freq_map_p2(:,:,ii) - Season_freq_map_p1(:,:,ii),'AlphaData',~isnan(Season_freq_map_p2(:,:,ii) - Season_freq_map_p1(:,:,ii))); 
title([season_labels{ii}]); set(gca,'YDir','normal')
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.8 0.8 0.8]) 
clim([-0.4 0.4]); axis image; cb1 = colorbar; ylabel(cb1,'Snow presence anomaly [-]','FontSize',11); colormap(n3,flipud(redblue)); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons(:)) max(demLons(:))]); ylim([min(demLats(:)) max(demLats(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
xlabel([num2str(min(year(Date_l8s2(ind_seas_p2(:,ii))))) '-'  num2str(max(year(Date_l8s2(ind_seas_p2(:,ii))))) ' vs ' ...
    num2str(min(year(Date_l8s2(ind_seas_p1(:,ii))))) '-'  num2str(max(year(Date_l8s2(ind_seas_p1(:,ii)))))])

exportgraphics(fi2,[dir_fig '\Snow_cover\' glacier '_' season_labels{ii} '_SnowFreqMaps_2periods_comp.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
end

%% Analyze monthly changes in Sentinel-2/ Landsat_8 snow cover

for ii = 1:12

ind_mm_p1 = (month(Date_l8s2) == ii)  & Date_l8s2 < middle_date;
ind_mm_p2 = (month(Date_l8s2) == ii)  & Date_l8s2 >= middle_date;
    
ind_nan_p1 = sum(isnan(l8s2_sno(:,:,ind_mm_p1)),3)<(0.33*size(l8s2_sno(:,:,ind_mm_p1),3));    
ind_nan_p2 = sum(isnan(l8s2_sno(:,:,ind_mm_p2)),3)<(0.33*size(l8s2_sno(:,:,ind_mm_p2),3));    

fi2 = figure('Renderer', 'painters', 'Position',[89.6667 279 978.6666 374.6666]) ;
tiledlayout(1,3,'TileSpacing','compact')
n1 = nexttile;
imagesc(demLons_l8s2(1,:), demLats_l8s2(:,1),nanmean(l8s2_sno(:,:,ind_mm_p1),3),'AlphaData',MASK_hls & ind_nan_p1); 
title([month_labels{ii} ' (' num2str(sum(ind_mm_p1)) ': scenes)']); set(gca,'YDir','normal')
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.8 0.8 0.8]) 
clim([0 1]); axis image; cb1 = colorbar; ylabel(cb1,'Mean snow presence [-]','FontSize',11); colormap(n1,cmap_snow_freq); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons(:)) max(demLons(:))]); ylim([min(demLats(:)) max(demLats(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
xlabel([num2str(min(year(Date_l8s2(ind_mm_p1)))) '-'  num2str(max(year(Date_l8s2(ind_mm_p1))))])

n2 = nexttile;
imagesc(demLons_l8s2(1,:), demLats_l8s2(:,1), nanmean(l8s2_sno(:,:,ind_mm_p2),3),'AlphaData',MASK_hls & ind_nan_p2); 
title([month_labels{ii} ' (' num2str(sum(ind_mm_p2)) ': scenes)']); set(gca,'YDir','normal')
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.8 0.8 0.8]) 
clim([0 1]); axis image; cb1 = colorbar; ylabel(cb1,'Mean snow presence [-]','FontSize',11); colormap(n2,cmap_snow_freq); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons(:)) max(demLons(:))]); ylim([min(demLats(:)) max(demLats(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
xlabel([num2str(min(year(Date_l8s2(ind_mm_p2)))) '-'  num2str(max(year(Date_l8s2(ind_mm_p2))))])

n3 = nexttile;
imagesc(demLons_l8s2(1,:), demLats_l8s2(:,1),nanmean(l8s2_sno(:,:,ind_mm_p2),3) - nanmean(l8s2_sno(:,:,ind_mm_p1),3),'AlphaData',MASK_hls & ind_nan_p1 & ind_nan_p2); 
title([month_labels{ii} ' (' num2str(sum(ind_mm_p2)) ': scenes)']); set(gca,'YDir','normal')
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.8 0.8 0.8]) 
clim([-0.5 0.5]); axis image; cb1 = colorbar; ylabel(cb1,'Mean snow presence [-]','FontSize',11); colormap(n3,flipud(redblue)); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons(:)) max(demLons(:))]); ylim([min(demLats(:)) max(demLats(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
xlabel([num2str(min(year(Date_l8s2(ind_mm_p2)))) '-'  num2str(max(year(Date_l8s2(ind_mm_p2)))) ' vs ' ...
    num2str(min(year(Date_l8s2(ind_mm_p1)))) '-'  num2str(max(year(Date_l8s2(ind_mm_p1))))])

exportgraphics(fi2,[dir_fig '\Snow_cover\' glacier '_' month_labels{ii} '_SnowFreqMaps_2periods_comp.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
end

%% Load MODIS daily snow cover for the same extent as Sentinel-2 and Landsat-8 maps

if strcmp(glacier,'Kyzylsu')

site_name = 'Tajikistan';
foSC = [root_B '\Remote_sensing\MODIS\Data\' site_name]; %raw data
foPSC = [root_B '\Remote_sensing\MODIS\Processed_data\' site_name]; % processed data

MODIS_folders =  dir([foPSC '\20*']);
MODIS_years = str2double(string({MODIS_folders.name}));

% Load first scene to get geospatial information

dataset = 'MOD10A1';
fns = dir([foSC '\' dataset '.061_NDSI_Snow_Cover_doy*.tif']);
[SCENE1, SCENE1_info] = readgeoraster([foSC '\' fns(30).name]);
SCENE1_INFO = georasterinfo([foSC '\' fns(30).name]);

imInfo = geotiffinfo([foSC '\' fns(30).name]);
imRM = imInfo.RefMatrix;
imSize = [imInfo.Height,imInfo.Width];
[imLons,imLats] = pixcenters(imRM,imSize);
[imLons,imLats] = meshgrid(imLons,imLats);
minLon = min(imLons(:)); maxLon = max(imLons(:)); 
minLat = min(imLats(:)); maxLat = max(imLats(:)); 

SCENE1_gridobj = GRIDobj(imLons,imLats,SCENE1(:,:,1));  % convert it to a GRIDobj to reproject the AW3D on it.

% Load all MODIS scenes corresponding to the landsat and sentinel-2 dates

dataset = 'MOD10A1';
fns = dir([foSC '\' dataset '.061_NDSI_Snow_Cover_doy*.tif']);
ndsiT = 0.4;

% Get the datestamps of all images
dateTime = NaT(size(fns,1),1);

for iIm = 1:size(fns,1)
    fn = fns(iIm).name;
    
    yr = fn(end-18:end-15);
    doy = fn(end-14:end-12);
    yr = str2num(yr);
    doy = str2num(doy);
    dateTime(iIm) = datetime(yr,1,doy);
end 

L8s2_mod_dates= intersect(dateTime,Date_l8s2_d,'rows');
i_mod = ismember(dateTime, L8s2_mod_dates);

fns = fns(i_mod);
date_mod = dateTime(i_mod);
nIms = length(fns);
[ims] = deal(nan(size(SCENE1,1),size(SCENE1,2),nIms));

for iIm = 1:nIms
    iIm
    % Get filename
    fn = fns(iIm).name;

    % Check file size
    s=dir([foSC '\' fn]);
    the_size=s.bytes;

    if the_size == 0
        continue
    end

    % Get image
    im = geotiffread([foSC '\' fn]);
    im = double(im);
    im(im > 100) = NaN; % Unknown
    im(im < ndsiT) = 0; % No snow
    im(im >= ndsiT) = 1; % Snow
    ims(:,:,iIm) = im;
end 

% find intersect between MODIS dates, T&C simulation and LS/S2 dates

L8s2_mod_dates= intersect(date_mod,Date_l8s2_d,'rows');
L8s2_TC_dates = intersect(Date_d,Date_l8s2_d,'rows');

i_l8s2 = ismember(Date_l8s2_d, L8s2_mod_dates);
i_mod = ismember(date_mod, L8s2_mod_dates);
i_tc = ismember(Date_d, L8s2_TC_dates); % Indices of TC dates for which we have landsat-Sentinel-2 snow cover maps

% Compare two different periods

middle_date = datetime(2018,9,30,0,0,0);

for mm = 1:12
    ind_p1 = (date_mod < middle_date) & (ismember(month(date_mod), mm));
    ind_p2 = (date_mod >= middle_date) & (ismember(month(date_mod), mm));

    ind_p1_l8s2 = (Date_l8s2_d < middle_date) & (ismember(month(Date_l8s2_d), mm));
    ind_p2_l8s2 = (Date_l8s2_d >= middle_date) & (ismember(month(Date_l8s2_d), mm));

    ind_p1_tc = (Date_d < middle_date) & (ismember(month(Date_d), mm));
    ind_p2_tc = (Date_d >= middle_date) & (ismember(month(Date_d), mm));

    MODIS_sp_p1_m(:,:,mm) = nanmean(ims(:,:,ind_p1 & i_mod),3);
    MODIS_sp_p2_m(:,:,mm) = nanmean(ims(:,:,ind_p2 & i_mod),3);

    L8S2_sp_p1_m(:,:,mm) = nanmean(l8s2_sno(:,:,ind_p1_l8s2 & i_l8s2),3);              
    L8S2_sp_p2_m(:,:,mm) = nanmean(l8s2_sno(:,:,ind_p2_l8s2 & i_l8s2),3);     

    TC_sp_p1_m(:,:,mm) = nanmean(snow_pres(:,:,ind_p1_tc & i_tc),3);              
    TC_sp_p2_m(:,:,mm) = nanmean(snow_pres(:,:,ind_p2_tc & i_tc),3);  
end 

cmap_snow_freq = cbrewer('seq','Blues',80); cmap_snow_freq(cmap_snow_freq<0) = 0; cmap_snow_freq(cmap_snow_freq>1) = 1;
cmap_snow_freq(1,:) = [1, 1, 1];

season_comp = [3 4 5 6];
season_comp_label = 'Mar-Jun';

MODIS_obj = SCENE1_gridobj;
MODIS_obj.Z = nanmean(ims(:,:,ismember(month(date_mod), season_comp)),3);

Diff_obj = SCENE1_gridobj;
Diff_obj.Z = nanmean(MODIS_sp_p2_m(:,:,season_comp),3) - nanmean(MODIS_sp_p1_m(:,:,season_comp),3);
Diff_obj_utm = reproject2utm(Diff_obj,500,'zone','42N');

Diff_obj_tc = DEM;
Diff_obj_tc.Z = flipud(nanmean(TC_sp_p2_m(:,:,season_comp),3) - nanmean(TC_sp_p1_m(:,:,season_comp),3));

Diff_obj_l8s2 = DEM_l8s2_orig; ii = 3;
Diff_obj_l8s2.Z = nanmean(L8S2_sp_p2_m(:,:,season_comp),3) - nanmean(L8S2_sp_p1_m(:,:,season_comp),3);

%% Anomaly map, MODIS vs Sentinel vs T&C

cmap_fsc_anom = cbrewer('div', 'BrBG', 80); cmap_fsc_anom(cmap_fsc_anom<0)=0; cmap_fsc_anom(cmap_fsc_anom>1)=1;

fi2 = figure('Renderer', 'painters', 'Position',[50 57.6667 988.6667 410.0000]) ;
tiledlayout(1,3,"TileSpacing","tight")
n1 = nexttile;
imagesc(Diff_obj_l8s2,'AlphaData',~isnan(Diff_obj_l8s2.Z)); 
set(gca,'YDir','normal')
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.8 0.8 0.8]) 
clim([-0.3 0.3]); axis image; colormap(n1,flipud(redblue)); hold on;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
title('Landsat & Sentinel-2','FontWeight','normal','FontSize',12)

n2 = nexttile;
imagesc(Diff_obj_tc,'AlphaData',flipud(MASK == 1)); hold on; 
clim([-0.3 0.3]); axis image; set(gca,'Color',[0.8 0.8 0.8])
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
colormap(cmap_fsc_anom); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
% xlabel([num2str(min(year(Date_l8s2(ind_seas_p2(:,ii))))) '-'  num2str(max(year(Date_l8s2(ind_seas_p2(:,ii))))) ' vs ' ...
%     num2str(min(year(Date_l8s2(ind_seas_p1(:,ii))))) '-'  num2str(max(year(Date_l8s2(ind_seas_p1(:,ii)))))])
xlabel(['2018-2023 vs 2000-2018, ' season_comp_label])
title('Simulated','FontWeight','normal','FontSize',12)

n3 = nexttile; 
imagesc(Diff_obj_utm); hold on; 
clim([-0.3 0.3]); axis image;
plot(catchShp_utm.X, catchShp_utm.Y,'k'); xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
cb1 = colorbar; ylabel(cb1,'Snow persistence anomaly [-]','FontSize',12);
colormap(cmap_fsc_anom); hold on;
title('MODIS','FontWeight','normal','FontSize',12)

exportgraphics(fi2,[dir_fig '\Snow_cover\' glacier '_' season_comp_label '_SnowFreqMaps_2periods_comp_MODIS_' num2str(ndsiT) '_L8S2_brown.png'],'Resolution',300,'BackgroundColor','none')

end

%%

for mm = 1:12 % For all of T&C simulated days (not for strict comparison with MODIS or Landsat

    ind_p1_tc = (Date_d < middle_date) & (ismember(month(Date_d), mm));
    ind_p2_tc = (Date_d >= middle_date) & (ismember(month(Date_d), mm));
    
    TC_sp_p1_all_m(:,:,mm) = nanmean(snow_pres(:,:,ind_p1_tc),3);              
    TC_sp_p2_all_m(:,:,mm) = nanmean(snow_pres(:,:,ind_p2_tc),3);  

    TC_snd_p1_all_m(:,:,mm) = nanmean(snow_depth(:,:,ind_p1_tc),3);              
    TC_snd_p2_all_m(:,:,mm) = nanmean(snow_depth(:,:,ind_p2_tc),3);  
end 

%% Show all for seasons of simulated snow persistence anomaly

Season_months = [10, 11, 12; ...
                 1, 2, 3; ...
                 4, 5, 6; ...
                 7, 8, 9];

fi2 = figure('Renderer', 'painters', 'Position',[50 50 891.3333 923.3333]) ;
tiledlayout(3,4,"TileSpacing","compact")

for ss = 1:size(Season_months,1)
n1 = nexttile;
Diff_obj_tc_ss = DEM;
Diff_obj_tc_ss.Z = flipud(nanmean(TC_sp_p1_all_m(:,:,Season_months(ss,:)),3));

imagesc(Diff_obj_tc_ss,'AlphaData',flipud(MASK == 1)); hold on; 
clim([0 1]); axis image; set(gca,'Color',[0.8 0.8 0.8])
plot(catchShp_utm.X, catchShp_utm.Y,'k'); 
axis image;
xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
colormap(n1, cmap_snow_freq); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
title(season_labels{ss},'FontWeight','normal','FontSize',12); 
if ss==1; ylabel('2000-2018','FontWeight','bold','FontSize',12); end

if ss == 4
    cb1 = colorbar('Location','eastoutside'); ylabel(cb1,'Mean Snow presistence [-]','FontSize',12); 
end 
end

for ss = 1:size(Season_months,1)
n1 = nexttile;
Diff_obj_tc_ss = DEM;
Diff_obj_tc_ss.Z = flipud(nanmean(TC_sp_p2_all_m(:,:,Season_months(ss,:)),3));

imagesc(Diff_obj_tc_ss,'AlphaData',flipud(MASK == 1)); hold on; 
clim([0 1]); axis image; set(gca,'Color',[0.8 0.8 0.8])
plot(catchShp_utm.X, catchShp_utm.Y,'k'); 
axis image;
xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
colormap(n1, cmap_snow_freq); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
% title(season_labels{ss},'FontWeight','normal','FontSize',12); 
if ss==1; ylabel('2018-2023','FontWeight','bold','FontSize',12); end

if ss == 4
    cb2 = colorbar('Location','eastoutside'); ylabel(cb2,'Mean Snow presistence [-]','FontSize',12); 
end 
end


for ss = 1:size(Season_months,1)
n2 = nexttile;
Diff_obj_tc_ss = DEM;
Diff_obj_tc_ss.Z = flipud(nanmean(TC_sp_p2_all_m(:,:,Season_months(ss,:)),3) - nanmean(TC_sp_p1_all_m(:,:,Season_months(ss,:)),3));

imagesc(Diff_obj_tc_ss,'AlphaData',flipud(MASK == 1)); hold on; 
clim([-0.5 0.5]); axis image; set(gca,'Color',[0.8 0.8 0.8])
plot(catchShp_utm.X, catchShp_utm.Y,'k'); 
axis image;
xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
colormap(n2,flipud(redblue)); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
if ss== 1; ylabel('2018-2023 - 2000-2018','FontWeight','bold','FontSize',12); end

if ss == 4
    cb3 = colorbar('Location','eastoutside'); ylabel(cb3,'\DeltaSnow presistence [-]','FontSize',12); 
end 
end

exportgraphics(fi2,[dir_fig '\Snow_cover\' glacier '_seasonal_TC_SnowFreqMaps_2periods_anomaly.png'],'Resolution',300,'BackgroundColor','none')

%% Show all for seasons of simulated snow mean snow depth

Season_months = [10, 11, 12; ...
                 1, 2, 3; ...
                 4, 5, 6; ...
                 7, 8, 9];

fi2 = figure('Renderer', 'painters', 'Position',[50 50 891.3333 923.3333]) ;
tiledlayout(3,4,"TileSpacing","compact")

for ss = 1:size(Season_months,1)
n1 = nexttile;
Diff_obj_tc_ss = DEM;
Diff_obj_tc_ss.Z = flipud(nanmean(TC_snd_p1_all_m(:,:,Season_months(ss,:)),3));

imagesc(Diff_obj_tc_ss,'AlphaData',flipud(MASK == 1)); hold on; 
clim([0 2]); axis image; set(gca,'Color',[0.8 0.8 0.8])
plot(catchShp_utm.X, catchShp_utm.Y,'k'); 
axis image;
xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
colormap(n1, cmap_snow_freq); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
title(season_labels{ss},'FontWeight','normal','FontSize',12); 
if ss==1; ylabel('2000-2018','FontWeight','bold','FontSize',12); end

if ss == 4
    cb1 = colorbar('Location','eastoutside'); ylabel(cb1,'Mean Snow Height [m]','FontSize',12); 
end 
end

for ss = 1:size(Season_months,1)
n1 = nexttile;
Diff_obj_tc_ss = DEM;
Diff_obj_tc_ss.Z = flipud(nanmean(TC_snd_p2_all_m(:,:,Season_months(ss,:)),3));

imagesc(Diff_obj_tc_ss,'AlphaData',flipud(MASK == 1)); hold on; 
clim([0 2]); axis image; set(gca,'Color',[0.8 0.8 0.8])
plot(catchShp_utm.X, catchShp_utm.Y,'k'); 
axis image;
xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
colormap(n1, cmap_snow_freq); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
% title(season_labels{ss},'FontWeight','normal','FontSize',12); 
if ss==1; ylabel('2018-2023','FontWeight','bold','FontSize',12); end

if ss == 4
    cb2 = colorbar('Location','eastoutside'); ylabel(cb2,'Mean Snow Height [m]','FontSize',12); 
end 
end

for ss = 1:size(Season_months,1)
n2 = nexttile;
Diff_obj_tc_ss = DEM;
Diff_obj_tc_ss.Z = flipud(nanmean(TC_snd_p2_all_m(:,:,Season_months(ss,:)),3) - nanmean(TC_snd_p1_all_m(:,:,Season_months(ss,:)),3));

imagesc(Diff_obj_tc_ss,'AlphaData',flipud(MASK == 1)); hold on; 
clim([-1.5 1.5]); axis image; set(gca,'Color',[0.8 0.8 0.8])
plot(catchShp_utm.X, catchShp_utm.Y,'k'); 
axis image;
xlim([min(demLons_l8s2(:))+3000 max(demLons_l8s2(:))]); ylim([min(demLats_l8s2(:)) max(demLats_l8s2(:))])
for gg = 1:size(glaciers,1), plot(glaciers(gg).X,glaciers(gg).Y,'Color',[0.7 0.7 0.7]); end 
colormap(n2,flipud(redblue)); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
% title(season_labels{ss},'FontWeight','normal','FontSize',12); 
if ss== 1; ylabel('2018-2023 - 2000-2018','FontWeight','bold','FontSize',12); end

if ss == 4
    cb3 = colorbar('Location','eastoutside'); ylabel(cb3,'\DeltaSnow Height [m]','FontSize',12); 
end 
end

exportgraphics(fi2,[dir_fig '\Snow_cover\' glacier '_seasonal_TC_SnowHeight_2periods_anomaly.png'],'Resolution',300,'BackgroundColor','none')


%% Topographic parameters vs snow cover frequency  [Uncomment if useful ...]

% % Observed
% dEL = 100;
% min_EL = nanmin(DEM_l8s2_orig.Z(MASK_hls==1));
% max_EL = nanmax(DEM_l8s2_orig.Z(MASK_hls==1));
% ELs = (min_EL-100):dEL:(max_EL+100);
% for ii = 1:size(ind_seas,2)
% for iel = 1:numel(ELs)
%     cur = (DEM_l8s2_orig.Z<(ELs(iel)+dEL/2))&(DEM_l8s2_orig.Z>=(ELs(iel)-dEL/2)) & (MASK_hls==1); %current section of DEM
%     Seas_mean_sno = nanmean(l8s2_sno(:,:,ind_seas(:,ii)),3);    
%     Seas_std_sno = nanstd(l8s2_sno(:,:,ind_seas(:,ii)),0,3);
%     SF_seas_mean_el(iel,ii) = nanmean(Seas_mean_sno(cur));
%     SF_seas_std_el(iel,ii) = nanmean(Seas_std_sno(cur));
% end
% end
% 
% % Modelled
% 
% ind_win_tc = (month(Date_d) > 11 | month(Date_d) < 3) & ismember(Date_d,Date_l8s2);
% ind_spr_tc = (month(Date_d) > 2 & month(Date_d) < 6) & ismember(Date_d,Date_l8s2);
% ind_sum_tc = (month(Date_d) > 5 & month(Date_d) < 9) & ismember(Date_d,Date_l8s2);
% ind_aut_tc = (month(Date_d) > 8 & month(Date_d) < 12) & ismember(Date_d,Date_l8s2);
% 
% ind_seas_tc = [ind_aut_tc, ind_win_tc, ind_spr_tc, ind_sum_tc];
% 
% dEL = 100;
% ELs_tc = nanmin(DTM,[],'all'):dEL:nanmax(DTM,[],'all');
% 
% for ii = 1:size(ind_seas,2)
% for iel = 1:numel(ELs_tc)
%     cur = (DTM<(ELs_tc(iel)+dEL/2))&(DTM>=(ELs_tc(iel)-dEL/2)) & (DTM>0); %current section of DEM
%     Seas_mean_sno_tc = nanmean(snow_pres(:,:,ind_seas_tc(:,ii)),3);    
%     Seas_std_sno_tc = nanstd(snow_pres(:,:,ind_seas_tc(:,ii)),0,3);
%     SF_seas_mean_el_tc(iel,ii) = nanmean(Seas_mean_sno_tc(cur));
%     SF_seas_std_el_tc(iel,ii) = nanmean(Seas_std_sno_tc(cur));
% end
% end
% 
% %%
% fi2 = figure('Renderer', 'painters', 'Position',[253.6667 41.6667 163.3333 678.0000]) ;
% tiledlayout(4,1,"TileSpacing","compact")
% for ii = 1:size(ind_seas,2)
% nexttile
% errorbar(ELs./1000,SF_seas_mean_el(:,ii),SF_seas_std_el(:,ii),'r'); ylim([0 1]); hold on; grid on;
% plot(ELs./1000,SF_seas_mean_el(:,ii),'k','LineWidth',1.2)
% plot(ELs_tc./1000,SF_seas_mean_el_tc(:,ii),'b','LineWidth',1.2)
% if ii == 4; legend('1 \sigma','Mean','Modelled','Location','south'); end
% view(90,-90); hold on; grid on;
% xlabel('Elevation [km]','FontSize',10); 
% if ii ~= 4 && ii ~=1 ; set(gca,'YTickLabels',[]);end
% if ii == 4; ylabel('Mean snow presence [-]','FontSize',10); end 
% if ii == 1; ylabel('Mean snow presence [-]','FontSize',10); set(gca,'YAxisLocation','right'); end 
% % title(season_labels{ii})
% xlim([min(ELs)./1000 max(ELs)./1000])
% text(5.8,0.05,season_labels{ii},'FontSize',8,'FontWeight','bold');
% end
% 
% exportgraphics(fi2,[dir_fig '\Snow_cover\' 'SentinelLandsat_SnowFreqMaps_MeanSTD_T&Ccomo.png'],'Resolution',300,'BackgroundColor','none')
% 
% 
% fi2 = figure('Renderer', 'painters', 'Position',[253.6667 41.6667 163.3333 678.0000]) ;
% tiledlayout(4,1,"TileSpacing","compact")
% for ii = 1:size(ind_seas,2)
% nexttile
% errorbar(ELs_tc./1000,SWE_y_el_seas_mean(:,ii),SWE_y_el_seas_std(:,ii),'r'); hold on; grid on;
% plot(ELs_tc./1000,SWE_y_el_seas_mean(:,ii),'k','LineWidth',1.2)
% if ii == 1; legend('1 \sigma','Mean','Location','southeast'); end
% view(90,-90); hold on; grid on; ylim([0 700])
% xlabel('Elevation [km]'); yticks(0:150:600)
% if ii ~= 4 && ii ~= 1; set(gca,'YTickLabels',[]);end
% if ii == 4; ylabel('Mean SWE [mm w.e.]'); end 
% if ii == 1; ylabel('Mean SWE [mm w.e.]'); set(gca,'YAxisLocation','right'); end 
% text(5.8,100,season_labels{ii},'FontSize',8,'FontWeight','bold');
% xlim([min(ELs)./1000 max(ELs)./1000])
% end 
% exportgraphics(fi2,[dir_fig '\Snow_cover\Simulated_SWE_season_elevation_band_STD.png'],'Resolution',300,'BackgroundColor','none')


%% Reproject landsat/sentinel-2 snow cover maps onto the T&C grid, and
% compare with T&C outputs

l8s2_sno_proj = NaN(size(DEM.Z,1),size(DEM.Z,2),size(l8s2_sno,3));
l8s2_sno_comp = NaN(size(DEM.Z,1),size(DEM.Z,2),size(l8s2_sno,3));

for ii = 1:size(l8s2_sno,3)

   l8s2_sno_comp_i = NaN(size(DEM.Z,1),size(DEM.Z,2));

   l8s2_sno_pro_i = DEM_l8s2;
   l8s2_sno_pro_i.Z = l8s2_sno(:,:,ii);
   l8s2_sno_pro_i = resample(l8s2_sno_pro_i, DEM,'nearest');
   l8s2_sno_proj(:,:,ii) = l8s2_sno_pro_i.Z;

   l8s2_sno_comp_i((flipud(snow_pres(:,:,id_l8s2(ii,2))) == 0) & (l8s2_sno_proj(:,:,ii) == 0)) = 0;
   l8s2_sno_comp_i((flipud(snow_pres(:,:,id_l8s2(ii,2))) == 1) & (l8s2_sno_proj(:,:,ii) == 0)) = 1;
   l8s2_sno_comp_i((flipud(snow_pres(:,:,id_l8s2(ii,2))) == 0) & (l8s2_sno_proj(:,:,ii) == 1)) = 2;
   l8s2_sno_comp_i((flipud(snow_pres(:,:,id_l8s2(ii,2))) == 1) & (l8s2_sno_proj(:,:,ii) == 1)) = 3;
   l8s2_sno_comp_i(flipud(MASK) ~= 1) = NaN;
   l8s2_sno_comp(:,:,ii) = l8s2_sno_comp_i;

   % Compute DICE coefficient
    TP = sum(flipud(MASK)==1 & flipud(snow_pres(:,:,id_l8s2(ii,2)))==1 & l8s2_sno_proj(:,:,ii)==1,'all');
    FP = sum(flipud(MASK)==1 & flipud(snow_pres(:,:,id_l8s2(ii,2)))==1 & l8s2_sno_proj(:,:,ii)==0,'all');
    FN = sum(flipud(MASK)==1 & flipud(snow_pres(:,:,id_l8s2(ii,2)))==0 & l8s2_sno_proj(:,:,ii)==1,'all');
    DICE(ii) = 2*TP/(2*TP+FP+FN);
end 

%% Figure of Landsat-Sentinel snow cover vs T&C snow depth

if monthly_LS_figure == 1

c4 = colormap(getPyPlot_cMap('GnBu'));
c_class = [0.3594    0.2500    0.1992; 0 0 1; 1 0 0; 1 1 1];

for ii = 1:size(id_l8s2,1)
fi2 = figure('Renderer', 'painters', 'Position', [185.6667 309.6667 1.0027e+03 387.3333]) ;
tiledlayout(1,3,'TileSpacing','compact');
s3 = nexttile;
imagesc(x,y,snow_depth(:,:,id_l8s2(ii,2)),'AlphaData',DTM > 0); %clim([0 1]); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.9 0.9 0.9])
title(['T&C with ' forcing_nm ' : Snow depth'])
set(gca,'YDir','normal'); set(gca,'YAxisLocation','Left')
colormap(s3,[1,1,1; c4])
cb = colorbar('Location','westoutside'); 
ylabel(cb,'Snow depth [m]','FontSize',11)
clim([0 1.5]); xlim([min(x) max(x)]); ylim([min(y) max(y)])
xlabel(['fsc = ' num2str(round(scas(id_l8s2(ii,2)),2)) ', SLA = ' num2str(round(SPAVG_dm.SLE(id_l8s2(ii,2)+1),-1)) ' m'] )

nexttile
imagesc(x_l8s2(1,:),y_l8s2(:,1),cat(3,l8s2_r(:,:,ii), l8s2_g(:,:,ii), l8s2_b(:,:,ii)),'AlphaData',l8s2_r(:,:,ii) ~= -0.9999); clim([0 1]); hold on;
[C, h] = contour(x,y,DTM,[2000:500:6000],'Color',[1 0.6 0.6]);
clabel(C,h,[2000:500:6000],'LabelSpacing',1000,'Color',[1 0.6 0.6],'FontWeight','bold')
plot(catchShp_utm.X,catchShp_utm.Y,'r','LineWidth',1.1)
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.9 0.9 0.9])
text(4.4771*10^5, 3.08909*10^6,datestr(date_h(id_l8s2(ii,2)),'dd-mmm-yYYy'),'FontSize',11,'FontWeight','bold')
xlabel([datestr(date_h(id_l8s2(ii,2)),'dd-mmm-yYYy') ', fsc = ' num2str(round(LS_sno.fsnow_hybrid(id_l8s2(ii,1)),2)) ', SLA = ' num2str(round(LS_sno.se_hybrid(id_l8s2(ii,1)),-1)) ' m'] )
xlim([min(x) max(x)]); ylim([min(y) max(y)])
title([sat_l8s2{id_l8s2(ii,1)} ' - RGB'])
set(gca,'YDir','normal')

c3 = nexttile;
imagesc(x,y,flipud(l8s2_sno_comp(:,:,ii)),'AlphaData',flipud(~isnan(l8s2_sno_comp(:,:,ii)))); clim([0 3]); hold on;
plot(catchShp_utm.X,catchShp_utm.Y,'r','LineWidth',1.1); set(gca,'YDir','normal');
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.9 0.9 0.9])
title('Snow cover'); colormap(c3,c_class); cb = colorbar;
set(cb,'YTick',0:1:3); cb.TickLabels = {'No snow','Mod-NoObs','Obs-NoMod','Snow'};
xlabel(['DICE coefficient: ' num2str(round(DICE(ii),2))])
exportgraphics(fi2,[dir_fig '\Snow_cover\' glacier '_SnowCover_L8S2_validation_' datestr(date_h(id_l8s2(ii,2)),'yYYy-mm-dd') ... 
    '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
end 
end

%% Snow frequency map, all season

date_l8s2_tc = date_l8s2(id_l8s2(:,1));

[C ia ic] = unique([year(date_l8s2_tc)' month(date_l8s2_tc)'],'rows');
date_l8s2_monthly = date_l8s2_tc(ia);

Date_d = datetime(year(date_h),month(date_h),day(date_h));

for dd = 1:length(date_l8s2_monthly)
    ind_l8s2_tc_monthly(dd) = find(Date_d == date_l8s2_monthly(dd),1) ;
    ind_l8s2_monthly(dd) = find(date_l8s2_tc == date_l8s2_monthly(dd),1) ;
end 

Snow_frequency_hls_tc = nansum(snow_pres(:,:,ind_l8s2_tc_monthly),3)./size(snow_pres(:,:,ind_l8s2_tc_monthly),3);
Snow_frequency_hls = nansum(l8s2_sno(:,:,ind_l8s2_monthly) > 0,3)./nansum(~isnan(l8s2_sno(:,:,ind_l8s2_monthly)),3);

DEM_l8s2.size = size(DEM_l8s2.Z);
Snow_frequency_hls(nansum(l8s2_b,3)<-70) = NaN;
Snow_frequency_hls_tc(isnan(DTM))=NaN;

%Project obs snow cover onto T&C DEM
Snow_freq_obj = DEM_l8s2;
Snow_freq_obj.Z = Snow_frequency_hls;
Snow_freq_obj = resample(Snow_freq_obj,DEM);
Snow_freq_obj.Z(isnan(flipud(DTM)))=NaN;

cmap = redblue;

fi4 = figure('Renderer', 'painters', 'Position', [185 362.3333 865.3333 352.6667]);
tiledlayout(1,3,'TileSpacing','compact')
nexttile
imageschs(DEM,flipud(Snow_frequency_hls_tc),'colorbarylabel','Snow cover frequency [-]','caxis',[0 1],'colorbar',1)
cb = colorbar;
cb.Location = 'westoutside'; hold on;
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
xlim([min(x) max(x)]); ylim([min(y) max(y)])
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
title('T&C')
nexttile
imageschs(DEM_l8s2,Snow_frequency_hls,'colorbar',0,'caxis',[0 1]); hold on;
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
xlim([min(x) max(x)]); ylim([min(y) max(y)])
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
title('Landsat-8/Sentinel-2')
n3 = nexttile;
imageschs(DEM,flipud(Snow_frequency_hls_tc)-Snow_freq_obj.Z,'colorbar',1,'colormap',flipud(cmap(2:end,:)),'caxis',[-1 1],...
    'colorbarylabel','\DeltaSnow cover frequency [-]'); hold on;
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
xlim([min(x) max(x)]); ylim([min(y) max(y)])
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]); %
title('T&C - obs')
exportgraphics(fi4,[dir_fig '\T&C_SnowCoverFrequency_HLS.png'],'Resolution',300,'BackgroundColor','none')

%% Dice coefficient plot

fi4 = figure('Renderer', 'painters', 'Position', [185 318.3333 561.3333 396.6667]);
scatter(day(date_l8s2_tc, "dayofyear"),DICE,18,LS_sno.fsnow_hybrid(id_l8s2(:,1)),'filled');
cb = colorbar; ylabel(cb,'Observed fractional snow cover [-]','FontSize',11);
ylabel('DICE coefficient'); xlabel('Day of the year'); grid on; 
xlim([0 365])
exportgraphics(fi4,[dir_fig '\T&C_DICE_vs_fsc.png'],'Resolution',300,'BackgroundColor','none')

