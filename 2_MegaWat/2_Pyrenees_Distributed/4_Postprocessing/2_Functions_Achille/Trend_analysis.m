% needs:
% - Monthly or daily maps of: sublimation, precipitation, SWE, evaporation,
% snow melt, ice melt
% - Corresponding time vector (date_m)

if ~exist('trend_area','var')
    trend_area = 'catchment'; % 'catchment','glaciers';
end 

% Where to store spatially distributed trend figures

dir_fig_trends = [dir_fig '\Spatial_trends'];
if ~exist(dir_fig_trends, 'dir'); mkdir(dir_fig_trends); end 

SMB_map = PSNOW_map - SMSm_map - ESN_map ...
     - SMG_map;

Hydro_years = year(date_m(find(month(date_m) == 9 & day(date_m) == 1,1))):1:year(date_m(find(month(date_m) == 9 & day(date_m) == 1,1,'last')));

glaciers_id = unique(GLA_ID(MASK==1  & (GLA_ID>0)));
glaciers_numb = numel(glaciers_id);

[SWE_y, SWE_std, SSN_y, SSN_std, SMS_y, SMS_std, fsc_y, Precip_y, Precip_std, ...
   Psno_y, Psno_std,  Ta_y, Ta_std, WS_y, WS_std, RH_y, RH_std, SMB_y] = deal(NaN(numel(Hydro_years)-1,1));

[SWE_maps, SSN_maps, SMS_maps, Snow_freq_maps, Precip_maps, Psno_maps, ...
    Ta_maps, WS_maps, RH_maps, SMB_maps] = deal(NaN(size(SMS_map,1), size(SMS_map,2),numel(Hydro_years)-1));

% Mask used for the trend analysis
if strcmp(trend_area,'glacier') || strcmp(trend_area,'glaciers')
    mask_trend = MASK==1 &(GLA_ID>0);
else 
    mask_trend = MASK==1;
end 

for yy = 1:numel(Hydro_years)-1

    ind_start = find(year(date_m) == Hydro_years(yy) & month(date_m) == 9 & day(date_m) == 1,1);
    ind_end =   find(year(date_m) == Hydro_years(yy+1) & month(date_m) == 9 & day(date_m) == 1,1);
    ind_period = ind_start:ind_end;

    ind_start_d = find(year(date_d) == Hydro_years(yy) & month(date_d) == 9 & day(date_d) == 30,1);
    ind_end_d =   find(year(date_d) == Hydro_years(yy+1) & month(date_d) == 9 & day(date_d) == 30,1);
    ind_period_d = ind_start_d:ind_end_d;

    SWE_map_y = nanmean(SWE_map(:,:,ind_period),3); SWE_map_y(mask_trend==0)=NaN;
    SWE_y(yy) = nanmean(SWE_map_y,'all');
    SWE_std(yy) = nanstd(squeeze(nanmean(SWE_map(:,:,ind_period),[1 2])));
    SWE_maps(:,:,yy) = SWE_map_y;

    SSN_map_y = nansum(SSN_map(:,:,ind_period),3); SSN_map_y(mask_trend==0)=NaN;
    SSN_y(yy) = nanmean(SSN_map_y,'all');
    SSN_std(yy) = nanstd(SSN_map_y,[],'all');
    SSN_maps(:,:,yy) = SSN_map_y;

    SMS_map_y = nansum(SMS_map(:,:,ind_period),3); SMS_map_y(mask_trend==0)=NaN;
    SMS_y(yy) = nanmean(SMS_map_y,'all');
    SMS_std(yy) = nanstd(SMS_map_y,[],'all');
    SMS_maps(:,:,yy) = SMS_map_y;

    Snow_freq_y = nansum(snow_pres(:,:,ind_period_d),3)./numel(ind_period_d);
    Snow_freq_y(mask_trend==0)=NaN;
    Snow_freq_maps(:,:,yy) = Snow_freq_y;

    fsc_y(yy) = nanmean(squeeze(nansum(snow_pres(:,:,ind_period_d),[1 2]))./nCatchPix);

    Precip_map_y = nansum(PRECIP_map(:,:,ind_period),3); Precip_map_y(mask_trend==0)=NaN;
    Precip_y(yy) = nanmean(Precip_map_y,'all');
    Precip_std(yy) = nanstd(Precip_map_y,[],'all');
    Precip_maps(:,:,yy) = Precip_map_y;

    Psno_map_y = nansum(PSNOW_map(:,:,ind_period),3); Psno_map_y(mask_trend==0)=NaN;
    Psno_y(yy) = nanmean(Psno_map_y,'all');
    Psno_std(yy) = nanstd(Psno_map_y,[],'all');
    Psno_maps(:,:,yy) = Psno_map_y;

    Ta_map_y = nanmean(TA_map(:,:,ind_period),3); Ta_map_y(mask_trend==0)=NaN;
    Ta_y(yy) = nanmean(Ta_map_y,'all');
    Ta_std(yy) = nanstd(Ta_map_y,[],'all');
    Ta_maps(:,:,yy) = Ta_map_y;

    WS_map_y = nanmean(WS_map(:,:,ind_period),3); WS_map_y(mask_trend==0)=NaN;
    WS_y(yy) = nanmean(WS_map_y,'all');
    WS_std(yy) = nanstd(WS_map_y,[],'all');
    WS_maps(:,:,yy) = WS_map_y;

    RH_map_y = nanmean(RH_map(:,:,ind_period),3); RH_map_y(mask_trend==0)=NaN;
    RH_y(yy) = nanmean(RH_map_y,'all');
    RH_std(yy) = nanstd(RH_map_y,[],'all');
    RH_maps(:,:,yy) = RH_map_y;

    SMB_mapy = nansum(SMB_map(:,:,ind_period),3); SMB_mapy(GLA_ID == 0) = NaN;
    for gi = 1:numel(glaciers_id)
        SMB_y_gi(yy,gi) = nanmean(SMB_mapy(GLA_ID==glaciers_id(gi)),'all');
    end 
    SMB_y(yy) = nanmean(SMB_mapy,'all');
    SMB_maps(:,:,yy) = SMB_mapy;
    
%     SMB_y_std(yy) = nanstd(SMB_y_gi);
end 

%%
cmap_sub = cbrewer('seq','BuPu',20);
cmap_sub(cmap_sub<0) = 0;

fi2 = figure('Renderer', 'painters', 'Position', [199.6667 263 729.3333 442.6667]);
tiledlayout(2,2,'TileSpacing','compact')
nexttile
plot(Hydro_years(2:end), SMB_y,'k','LineWidth',1.1); grid on; hold on;
plot(Hydro_years(2:end),movmean(SMB_y,10),'r','LineWidth',1)
ylabel('Glacier SMB [mm w.e./yr]')  
legend('Individual years','10yr mov-avg','Location','NorthWest')
for gi = 1:numel(glaciers_id)
     plot(Hydro_years(2:end), SMB_y_gi(:,gi),'Color',[0.7 0.7 0.7 0.5],'HandleVisibility','off');
end 
% ax = gca; set(gca,'YAxisLocation','right')
xlim([Hydro_years(2) Hydro_years(end)])
title([glacier ' ' trend_area])

nexttile
plot(Hydro_years(2:end), Psno_y,'k','LineWidth',1.1); grid on; hold on;
plot(Hydro_years(2:end),movmean(Psno_y,10),'r','LineWidth',1,'HandleVisibility','off')
plot(Hydro_years(2:end), Precip_y,'-','Color',[0 0 0 0.5],'LineWidth',0.5)
plot(Hydro_years(2:end),movmean(Precip_y,10),'Color',[1 0 0 0.5],'LineWidth',0.5,'HandleVisibility','off')
ylabel('Precipitation [mm w.e.]'); ylim([0 1.15*max(Precip_y)])
legend('Solid','Total','Location','southwest')
yyaxis right
plot(Hydro_years(2:end),movmean(Psno_y./Precip_y,10),'Color','b','LineWidth',1,'HandleVisibility','off'); 
ax = gca; ax.YAxis(2).Color = 'b';
ylabel('Annual snowfall ratio [-]'); ylim([0 1])
xlim([Hydro_years(2) Hydro_years(end)])
title([glacier ' ' trend_area])

nexttile
plot(Hydro_years(2:end), SSN_y,'k','LineWidth',1.1); grid on; hold on;
plot(Hydro_years(2:end),movmean(SSN_y,10),'r','LineWidth',1)
scatter(Hydro_years(2:end), SSN_y, 20, SSN_y./(Psno_y),'filled'); colormap(cmap_sub);
clim([0.0 0.30]); if strcmp(glacier,'Mugagangqiong'); clim([0.0 0.50]); end
cb = colorbar; ylabel(cb,'Sublimation/precipitation ratio [-]')
ylabel('Sublimation [mm w.e.]')    
xlim([Hydro_years(2) Hydro_years(end)]); 
if strcmp(trend_area, 'glaciers'); ylim([0 200]); end
if strcmp(trend_area, 'catchment'); ylim([0 160]); end

nexttile
shadedErrorBar(Hydro_years(2:end), SWE_y,SWE_std); grid on; hold on;
plot(Hydro_years(2:end),movmean(SWE_y,10),'r','LineWidth',1)
ylabel('SWE [mm w.e.]'); ylim([0 1.3*max(SWE_y)])
xlim([Hydro_years(2) Hydro_years(end)])
yyaxis right
plot(Hydro_years(2:end), fsc_y,'Color',[0 0 1 0.6],'LineStyle','-','LineWidth',0.6); grid on; hold on;
plot(Hydro_years(2:end),movmean(fsc_y,10),'Color',[0 0 1 0.8],'LineStyle','-','LineWidth',1)
ylim([max(round(min(fsc_y),1)-0.1,0) min(1,round(max(fsc_y),1)+0.1)])
ax = gca; ax.YAxis(2).Color = 'b';
ylabel('Snow cover fraction [-]')

exportgraphics(fi2,[dir_fig '\T&C_decadal_trends_' trend_area '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

disp(['Period average (max) sublimation to snowfall ratio: ' num2str(nanmean(SSN_y./(Psno_y))) ...
    ' (' num2str(nanmax(SSN_y./(Psno_y))) ')'])

disp(['Period average (max) sublimation to snowpack ablation ratio: ' num2str(nanmean(SSN_y./(SSN_y+SMS_y))) ...
    ' (' num2str(nanmax(SSN_y./(SSN_y+SMS_y))) ')'])

disp(['Period average (max) sublimation [mm w.e.]: ' num2str(nanmean(SSN_y)) ...
    ' (' num2str(nanmax(SSN_y)) ')'])

%% Sublimation variability

% figure
% tiledlayout(2,2,'tilespacing','compact')
% nexttile
% scatter(SSN_y,Psno_y,20,'filled'); grid on;
% xlabel('Sublimation [mm w.e./yr]'); ylabel('Snowfall [mm w.e./yr]')
% nexttile
% scatter(SSN_y,Ta_y,20,'filled'); grid on;
% xlabel('Sublimation [mm w.e./yr]'); ylabel('Mean air temperature [°C]')
% nexttile
% scatter(SSN_y,WS_y,20,'filled'); grid on;
% xlabel('Sublimation [mm w.e./yr]'); ylabel('Mean wind speed [m/s]')
% nexttile
% scatter(SSN_y,RH_y,20,'filled'); grid on;
% xlabel('Sublimation [mm w.e./yr]'); ylabel('Mean relative humidity [-]')
% 
% 
% figure
% tiledlayout(2,2,'tilespacing','compact')
% nexttile
% scatter(SSN_y./(Psno_y),Psno_y,20,'filled'); grid on;
% xlabel('Sublimation [% of snowfall]'); ylabel('Snowfall [mm w.e./yr]')
% nexttile
% scatter(SSN_y./(Psno_y),Ta_y,20,'filled'); grid on;
% xlabel('Sublimation [% of snowfall]'); ylabel('Mean air temperature [°C]')
% nexttile
% scatter(SSN_y./(Psno_y),WS_y,20,'filled'); grid on;
% xlabel('Sublimation [% of snowfall]'); ylabel('Mean wind speed [m/s]')
% nexttile
% scatter(SSN_y./(Psno_y),RH_y,20,'filled'); grid on;
% xlabel('Sublimation [% of snowfall]'); ylabel('Mean relative humidity [-]')
% 
% 
%% Maps and trends !!
cmap = cbrewer('seq','YlGnBu',50); cmap(cmap>1) = 1; cmap(cmap<0) = 0;

fi2 = figure('Renderer', 'painters', 'Position', [199.6667 435.6667 622.6666 270]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile
imagesc(x,flipud(y),flipud(nanmean(Psno_maps./Precip_maps,3)),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Snowfall ratio [-]'); clim([0 1]); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap);
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
nexttile
imagesc(x,flipud(y),flipud(trend(Psno_maps./Precip_maps)).*10,'AlphaData',~isnan(flipud(trend(Psno_maps./Precip_maps)))); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Snowfall ratio trend [-/dec]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, flipud(redblue));
% clim([-max(abs(trend(Psno_maps./Precip_maps)),[],'all') +max(abs(trend(Psno_maps./Precip_maps)),[],'all')])
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
% nexttile
% imagesc(x,flipud(y),flipud(trend(Psno_maps)),'AlphaData',~isnan(flipud(trend(Psno_maps)))); hold on;
% for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
% set(gca,'YDir','normal');
% title('Snowfall ratio trend [-/yr]'); set(gca,'Color',[0.8 0.8 0.8])
% colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, flipud(redblue));
% clim([-max(abs(trend(Psno_maps)),[],'all') +max(abs(trend(Psno_maps)),[],'all')])
% xlim([xmin_map xmax_map]); ylim([ymin_map ymax_map]); 
% xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
exportgraphics(fi2,[dir_fig_trends '\T&C_SnowfallRatio_spatial_trends.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%%
fi2 = figure('Renderer', 'painters', 'Position', [199.6667 435.6667 622.6666 270]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile
imagesc(x,flipud(y),flipud(nanmean(SWE_maps,3)),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('SWE [mm w.e.]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap);
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
clim([0 max(nanmean(SWE_maps,3),[],'all')]);
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
nexttile
imagesc(x,flipud(y),flipud(trend(SWE_maps)),'AlphaData',~isnan(flipud(trend(SWE_maps)))); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('SWE trend [mm w.e./yr]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, flipud(redblue));
clim([-max(abs(trend(SWE_maps)),[],'all') +max(abs(trend(SWE_maps)),[],'all')])
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
exportgraphics(fi2,[dir_fig_trends '\T&C_SWE_spatial_trends.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

fi2 = figure('Renderer', 'painters', 'Position', [199.6667 435.6667 622.6666 270]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile
imagesc(x,flipud(y),flipud(nanmean(SSN_maps,3)),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Sublimation [mm w.e.]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap);
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
clim([0 max(nanmean(SSN_maps,3),[],'all')]);
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
nexttile
imagesc(x,flipud(y),flipud(trend(SSN_maps)),'AlphaData',~isnan(flipud(trend(SSN_maps)))); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Sublimation trend [mm w.e./yr]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, flipud(redblue));
clim([-max(abs(trend(SSN_maps)),[],'all') +max(abs(trend(SSN_maps)),[],'all')])
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
exportgraphics(fi2,[dir_fig_trends '\T&C_Sublimation_spatial_trends.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%%
fi2 = figure('Renderer', 'painters', 'Position', [199.6667 435.6667 622.6666 270]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile
imagesc(x,flipud(y),flipud(nanmean(Snow_freq_maps,3)),'AlphaData',flipud(mask_trend)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Snow cover frequency [-]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap);
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
clim([0 max(nanmean(Snow_freq_maps,3),[],'all','omitnan')]);
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
nexttile
imagesc(x,flipud(y),flipud(trend(Snow_freq_maps)),'AlphaData',~isnan(flipud(trend(Snow_freq_maps)))); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Trends in snow cover frequency [-/yr]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, flipud(redblue));
% clim([-max(abs(trend(Snow_freq_maps)),[],'all') +max(abs(trend(Snow_freq_maps)),[],'all')])
clim([-0.006 0.006]); axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
exportgraphics(fi2,[dir_fig_trends '\T&C_Snowfreq_spatial_trends.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Trends in snowfall 

cmap = cbrewer('seq','YlGnBu',50); cmap(cmap>1) = 1; cmap(cmap<0) = 0;
cmap_trend = cbrewer('seq', 'YlOrRd', 80); cmap_trend(cmap_trend<0)=0; cmap_trend(cmap_trend>1)=1;

fi2 = figure('Renderer', 'painters', 'Position', [199.6667 262.3333 638.6666 443.3334]);
tiledlayout(1,2,'TileSpacing','tight');
n1 = nexttile;
imagesc(x,flipud(y),flipud(nanmean(Psno_maps,3)),'AlphaData',flipud(mask_trend)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Snowfall [mm/yr]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap); axis image
xlim([min(demLons(~isnan(flipud(DTM)))) max(demLons(~isnan(flipud(DTM))))]); 
ylim([min(demLats(~isnan(flipud(DTM)))) max(demLats(~isnan(flipud(DTM))))]); 
clim([0 max(nanmean(Psno_maps,3),[],'all')]);
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))]); axis image
n2 = nexttile;
imagesc(x,flipud(y),flipud(trend(Psno_maps)).*10,'AlphaData',~isnan(flipud(trend(Psno_maps)))); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Trends in snowfall [mm/dec]'); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(n2, flipud(cmap_trend));
% clim([-max(abs(trend(Snow_freq_maps)),[],'all') +max(abs(trend(Snow_freq_maps)),[],'all')])
% clim([-0.006 0.006])
clim([-100 0]); axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
xlabel([glacier ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))])
exportgraphics(fi2,[dir_fig_trends '\T&C_Snowfall_spatial_trends.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
