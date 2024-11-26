% Plot spatial summary maps of runoff generation and contribution to runoff
% for a given period (script argument) and a given varaible.

dir_fig_spatial_flux = [dir_fig '\Spatial_fluxes'];
if ~exist(dir_fig_spatial_flux, 'dir'); mkdir(dir_fig_spatial_flux); end 


if ~exist('date_start_flux_p1','var') || ~exist('date_end_flux_p1','var')
    date_start_flux_p1 = date_m(1);
    date_end_flux_p1 = date_m(floor(length(date_m)/2));
end 

if ~exist('date_start_flux_p2','var') || ~exist('date_end_flux_p2','var')
    date_start_flux_p2 = date_m(1+floor(length(date_m)/2));
    date_end_flux_p2 = date_m(end);
end 

ind_start_p1 = find(date_m == date_start_flux_p1,1);
ind_end_p1 = find(date_m == date_end_flux_p1,1,'last');
nYears_p1 = yearfrac(date_start_flux_p1,date_end_flux_p1,0);

ind_start_p2 = find(date_m == date_start_flux_p2,1);
ind_end_p2 = find(date_m == date_end_flux_p2,1,'last');
nYears_p2 = yearfrac(date_start_flux_p2,date_end_flux_p2,0);

% Compute the average over the first period per variable

Pr_liq_map_p1 = nansum(PRAIN_map(:,:,ind_start_p1:ind_end_p1),3); Pr_liq_map_p1(MASK ~=1) = NaN;
Imelt_map_p1 = nansum(SMG_map(:,:,ind_start_p1:ind_end_p1),3); Imelt_map_p1(MASK~=1) = NaN;
Smelt_map_p1 = nansum(SMSm_map(:,:,ind_start_p1:ind_end_p1),3); Smelt_map_p1(MASK ~=1) = NaN;
SSN_map_p1 = nansum(SSN_map(:,:,ind_start_p1:ind_end_p1),3); SSN_map_p1(MASK ~=1) = NaN;
ET_map_p1 = nansum(ET_map(:,:,ind_start_p1:ind_end_p1),3); ET_map_p1(MASK ~=1) = NaN;

FLUXliq_map_p1 = cat(3,Pr_liq_map_p1, Imelt_map_p1, Smelt_map_p1);
FLUXvap_map_p1 = cat(3,SSN_map_p1, ET_map_p1);

Snowmelt_contrib_p1 = (nansum(Smelt_map_p1./nansum(FLUXliq_map_p1,3),3));
Snowmelt_contrib_p1(MASK ~=1) = NaN;

FLUXphase_map_p1 = cat(3,Imelt_map_p1, Smelt_map_p1, SSN_map_p1, ET_map_p1);
% 
% % Classify the main contributor of liquid water per pixel of the area
[~, FLUXphase_max_ID_p1] = max(FLUXphase_map_p1,[],3); FLUXphase_max_ID_p1(MASK ~=1) = NaN;
% [~, FLUXliq_max_ID] = max(FLUXliq_map,[],3);
% [~, FLUXvap_max_ID] = max(FLUXvap_map,[],3);
% [~, FLUXsnow_max_ID] = max(cat(3,Smelt_y_map,SSN_y_map),[],3);

% Compute the average over the first period per variable

Pr_liq_map_p2 = nansum(PRAIN_map(:,:,ind_start_p2:ind_end_p2),3); Pr_liq_map_p2(MASK ~=1) = NaN;
Imelt_map_p2 = nansum(SMG_map(:,:,ind_start_p2:ind_end_p2),3); Imelt_map_p2(MASK~=1) = NaN;
Smelt_map_p2 = nansum(SMSm_map(:,:,ind_start_p2:ind_end_p2),3); Smelt_map_p2(MASK ~=1) = NaN;
SSN_map_p2 = nansum(SSN_map(:,:,ind_start_p2:ind_end_p2),3); SSN_map_p2(MASK ~=1) = NaN;
ET_map_p2 = nansum(ET_map(:,:,ind_start_p2:ind_end_p2),3); ET_map_p2(MASK ~=1) = NaN;

FLUXliq_map_p2 = cat(3,Pr_liq_map_p2, Imelt_map_p2, Smelt_map_p2);
FLUXvap_map_p2 = cat(3,SSN_map_p2, ET_map_p2);

Snowmelt_contrib_p2 = (nansum(Smelt_map_p2./nansum(FLUXliq_map_p2,3),3));
Snowmelt_contrib_p2(MASK ~=1) = NaN;

FLUXphase_map_p2 = cat(3,Imelt_map_p2, Smelt_map_p2, SSN_map_p2, ET_map_p2);
% 
% % Classify the main contributor of liquid water per pixel of the area
[~, FLUXphase_max_ID_p2] = max(FLUXphase_map_p2,[],3); FLUXphase_max_ID_p2(MASK ~=1) = NaN;

% var_flux = {'Rain','Icemelt','Snowmelt','Sublimation','ET'};
% FLUX_map = cat(3,Pr_liq_map, Imelt_y_map, Smelt_y_map, SSN_y_map, ET_y_map);
% 
% % Classify the main contributor of liquid water per pixel of the area
% [~, FLUX_max_ID] = max(FLUX_map,[],3);
% [~, FLUXliq_max_ID] = max(FLUXliq_map,[],3);
% [~, FLUXvap_max_ID] = max(FLUXvap_map,[],3);
% [~, FLUXsnow_max_ID] = max(cat(3,Smelt_y_map,SSN_y_map),[],3);

cmap_ratio = cbrewer('seq','YlGnBu',50); cmap_ratio(cmap_ratio>1) = 1; cmap_ratio(cmap_ratio<0) = 0;
cmap_snowmelt = cbrewer('seq','YlOrBr',50); cmap_snowmelt(cmap_snowmelt>1) = 1; cmap_snowmelt(cmap_snowmelt<0) = 0;

%% Snowmelt figure
cmap_diff = redblue; cmap_diff = cmap_diff(1:255,:);

Snowmelt_contrib= (nansum((Smelt_map_p1+Smelt_map_p2)./nansum(FLUXliq_map_p1+FLUXliq_map_p2,3),3));
Snowmelt_contrib(MASK ~=1) = NaN;

Snowmelt_contrib_p1= (nansum((Smelt_map_p1)./nansum(FLUXliq_map_p1,3),3));
Snowmelt_contrib_p1(MASK ~=1) = NaN;

Snowmelt_contrib_p2= (nansum((Smelt_map_p2)./nansum(FLUXliq_map_p2,3),3));
Snowmelt_contrib_p2(MASK ~=1) = NaN;
%%
fi2 = figure('Renderer', 'painters', 'Position', [90.3333 163.6667 1.0153e+03 435.3333]);
tiledlayout(1,3,'TileSpacing','compact')
n1 = nexttile;
imageschs(DEM,flipud(Snowmelt_contrib),'caxis',[0 1],'colormap',cmap_ratio); hold on;
cb1 = n1.Colorbar; set(cb1,'Location','eastoutside')
ylabel(cb1,'Snowmelt contribution to runoff [-]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
xlim([min(x) max(x)]); ylim([min(y) max(y)])
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
% title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p1))])
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p2))])

n2 = nexttile;
imageschs(DEM,flipud(Snowmelt_contrib_p2 - Snowmelt_contrib_p1),'caxis',[-0.3 0.3],'colormap',flipud(cmap_diff)); hold on;
cb2 = n2.Colorbar; set(cb2,'Location','eastoutside')
ylabel(cb2,'\Delta Snowmelt contribution to runoff [-]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
xlim([min(x) max(x)]); ylim([min(y) max(y)])
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
% title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p1))])
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) ' vs ' ...
    num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1))])
exportgraphics(fi2,[dir_fig_spatial_flux '\Spatial_SnowmeltContributionRunoff_' num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1)) ...
    '_' num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Snowfall contribution to runoff figure

fi2 = figure('Renderer', 'painters', 'Position', [199.6667 21 551.3333 675.3334]);
tiledlayout(2,2,'TileSpacing','compact')
n1 = nexttile;
imageschs(DEM,flipud(Smelt_map_p1./nYears_p1),'colorbar',1,'colormap',cmap_snowmelt,'caxis',[0 1000]); hold on;
cb1 = n1.Colorbar; set(cb1,'Location','westoutside'); ylabel(cb1,'Snowmelt [mm/yr]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p1))])
n2 = nexttile;
imageschs(DEM,flipud(Smelt_map_p2./nYears_p2),'colorbar',1,'colormap',cmap_snowmelt,'caxis',[0 1000]); hold on;
cb2 = n2.Colorbar; set(cb2,'Location','eastoutside'); ylabel(cb2,'Snowmelt [mm/yr]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) ' - ' num2str(year(date_end_flux_p2))])
n3 = nexttile;
imageschs(DEM,flipud(Smelt_map_p2./nYears_p2 - Smelt_map_p1./nYears_p1),'colorbar',1,....
    'colormap',flipud(cmap_diff),'caxis',[-300 300]); hold on;
cb3 = n3.Colorbar; set(cb3,'Location','westoutside'); ylabel(cb3,'\Delta Snowmelt [mm/yr]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) ' vs ' ...
    num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1))])
n1 = nexttile;
imageschs(DEM,flipud(Snowmelt_contrib),'caxis',[0 1],'colormap',cmap_ratio); hold on;
cb1 = n1.Colorbar; set(cb1,'Location','eastoutside')
ylabel(cb1,'Snowmelt contribution to runoff [-]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
xlim([min(x) max(x)]); ylim([min(y) max(y)])
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
% title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p1))])
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p2))])
exportgraphics(fi2,[dir_fig_spatial_flux '\Spatial_Snowmelt_' num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1)) ...
    '_' num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Final vapor ratio figure

vap_liq_ratio_p1 = nansum(FLUXvap_map_p1,3)./nansum(FLUXliq_map_p1,3);
vap_liq_ratio_p2 = nansum(FLUXvap_map_p2,3)./nansum(FLUXliq_map_p2,3);

fi2 = figure('Renderer', 'painters', 'Position', [71 228.3333 644.6667 435.3333]);
tiledlayout(1,2,'TileSpacing','compact')
n1 = nexttile;
imageschs(DEM,flipud(vap_liq_ratio_p1),'colorbar',1,'colormap',cmap_snowmelt,'caxis',[0 1]); hold on;
cb1 = n1.Colorbar; set(cb1,'Location','westoutside'); ylabel(cb1,'Vapor/Liquid flux ratio [-]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p1))])
n2 = nexttile;
imageschs(DEM,flipud(vap_liq_ratio_p2),'colorbar',1,'colormap',cmap_snowmelt,'caxis',[0 1]); hold on;
cb2 = n2.Colorbar; set(cb2,'Location','eastoutside'); ylabel(cb2,'Vapor/Liquid flux ratio [-]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) ' - ' num2str(year(date_end_flux_p2))])
% n2 = nexttile;
% imageschs(DEM,flipud((vap_liq_ratio_p2-vap_liq_ratio_p1)),'colorbar',1,'colormap',flipud(cmap_diff),'caxis',[-0.5 0.5]); hold on;
% cb2 = n2.Colorbar; set(cb2,'Location','eastoutside'); ylabel(cb2,'Snowmelt [mm/yr]','FontSize',11)
% plot(catchShp_utm.X,catchShp_utm.Y,'r','LineWidth',0.7)
% set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
% for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
% axis image
% xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
% ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
% title([num2str(year(date_start_flux_p2)) ' - ' num2str(year(date_end_flux_p2))])
exportgraphics(fi2,[dir_fig_spatial_flux '\Spatial_Vapor_to_Liquid_FluxRatio_' num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1)) ...
    '_' num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) '.png'],'Resolution',300,'BackgroundColor','none')
% close(gcf)

%% Final ET figure

fi2 = figure('Renderer', 'painters', 'Position', [71 254.3333 822 409.3333]);
tiledlayout(1,3,'TileSpacing','compact')
n1 = nexttile;
imageschs(DEM,flipud(ET_map_p1./nYears_p1),'colorbar',1,'colormap',[1, 1, 1; cmap_snowmelt],'caxis',[0 600]); hold on;
cb1 = n1.Colorbar; set(cb1,'Location','westoutside'); ylabel(cb1,'Evapotranspiration [mm/yr]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p1))])
n2 = nexttile;
imageschs(DEM,flipud(ET_map_p2./nYears_p2),'colorbar',0,'colormap',[1, 1, 1; cmap_snowmelt],'caxis',[0 600]); hold on;
% cb2 = n2.Colorbar; set(cb2,'Location','eastoutside'); ylabel(cb2,'Evapotranspiration [mm/yr]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) ' - ' num2str(year(date_end_flux_p2))])
n3 = nexttile;
imageschs(DEM,flipud((ET_map_p2./nYears_p2-ET_map_p1./nYears_p1)),'colorbar',1,'colormap',cmap_diff,'caxis',[-100 100]); hold on;
cb3 = n3.Colorbar; set(cb3,'Location','eastoutside'); ylabel(cb3,'\DeltaEvapotranspiration [mm/yr]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) ' minus ' ...
    num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1))])
exportgraphics(fi2,[dir_fig_spatial_flux '\Spatial_ET_comp_' num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1)) ...
    '_' num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)


%% Final sublimation figure

sub_melt_ratio_p1 = (SSN_map_p1)./nYears_p1;%./(SSN_map_p1+Smelt_map_p1);
sub_melt_ratio_p2 = (SSN_map_p2)./nYears_p2;%./(SSN_map_p2+Smelt_map_p2);

fi2 = figure('Renderer', 'painters', 'Position', [71 183.6667 832.6667 479.9999]);
tiledlayout(1,3,'TileSpacing','compact')
n1 = nexttile;
imageschs(DEM,flipud(sub_melt_ratio_p1),'colorbar',1,'colormap',cmap_snowmelt,'caxis',[0 200]); hold on;
cb1 = n1.Colorbar; set(cb1,'Location','westoutside'); ylabel(cb1,'Sublimation [mm/yr]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p1))])
n2 = nexttile;
imageschs(DEM,flipud(sub_melt_ratio_p2),'colorbar',0,'colormap',cmap_snowmelt,'caxis',[0 200]); hold on;
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) ' - ' num2str(year(date_end_flux_p2))])
n2 = nexttile;
imageschs(DEM,flipud((sub_melt_ratio_p2-sub_melt_ratio_p1)),'colorbar',1,'colormap',flipud(cmap_diff),'caxis',[-50 50]); hold on;
cb2 = n2.Colorbar; set(cb2,'Location','eastoutside'); ylabel(cb2,'\DeltaSublimation [mm/yr]','FontSize',11)
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) ' vs ' ...
    num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1))])
exportgraphics(fi2,[dir_fig_spatial_flux '\Spatial_sublimation_' num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1)) ...
    '_' num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Figure on classified major phase change process

cmap_flux = [0.65 0.65 0.65; [0.45 0.85 0.9]; [240 100 10]./256; 0.0510 0.7686 0];
var_phase = {'Icemelt','Snowmelt','Sublimation','ET'};

fi4 = figure('Renderer', 'painters', 'Position', [71 222.3333 685.3333 441.3333]);
tiledlayout(1,2,'TileSpacing','compact')
n1 = nexttile;
imageschs(DEM,flipud(FLUXphase_max_ID_p1),'colorbar',1,'colormap',cmap_flux); hold on;
cb1 = n1.Colorbar; set(cb1,'Location','westoutside'); ylabel(cb1,'Process dominating phase change','FontSize',12,'FontWeight','bold')
set(cb1,'XTick',[1.375 2.1250 2.8750 3.6250]); set(cb1,'XTickLabel',var_phase); cb1.FontSize = 10; 
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]); clim([1 4])
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p1)) ' - ' num2str(year(date_end_flux_p1))])
n2 = nexttile;
imageschs(DEM,flipud(FLUXphase_max_ID_p2),'colorbar',1,'colormap',cmap_flux); hold on;
cb2 = n2.Colorbar; set(cb2,'Location','eastoutside');
set(cb2,'XTick',[1.375 2.1250 2.8750 3.6250]); set(cb2,'XTickLabel',var_phase); cb2.FontSize = 10; 
plot(catchShp_utm.X,catchShp_utm.Y,'k','LineWidth',0.7)
set(gca,'YTickLabels',[]);set(gca,'XTickLabels',[]);
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7 0.7]); end 
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-400 max(demLons(~isnan(flipud(DTM))))+400]); 
ylim([min(demLats(~isnan(flipud(DTM))))-400 max(demLats(~isnan(flipud(DTM))))+400]); 
title([num2str(year(date_start_flux_p2)) ' - ' num2str(year(date_end_flux_p2))])

exportgraphics(fi4,[dir_fig_spatial_flux '\Major_PhaseChangeProcess_' num2str(year(date_start_flux_p1)) '-' num2str(year(date_end_flux_p1)) ...
    '_' num2str(year(date_start_flux_p2)) '-' num2str(year(date_end_flux_p2)) '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
