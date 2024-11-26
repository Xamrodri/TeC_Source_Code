
cmap = cbrewer('seq','YlGnBu',50); cmap(cmap>1) = 1; cmap(cmap<0) = 0;
cmap_red = cbrewer('seq','YlOrBr',50); cmap_red(cmap_red>1) = 1; cmap_red(cmap_red<0) = 0;
cmap_red = [1,1,1; cmap_red];

Tot_precip = nansum(PRECIP_map,3);
Tot_precip(~MASK) = NaN;

Tot_snowfall = nansum(PSNOW_map,3);
Tot_snowfall(~MASK) = NaN;

Snowfall_ratio = Tot_snowfall./Tot_precip;

Tot_snowmelt = nansum(SMSm_map,3);
Tot_snowmelt(~MASK) = NaN;

Tot_sublim = nansum(SSN_map,3);
Tot_sublim(~MASK) = NaN;

Tot_icemelt = nansum(SMG_map,3);
Tot_icemelt(~MASK) = NaN;

Tot_ET = nansum(ET_map,3);
Tot_SSN = nansum(SSN_map,3);


if length(Years_no) > 1
    perYear = yearfrac(date_m(1), date_m(end));
else 
    perYear = 1;
end


if strcmp(glacier,'Langtang')  
   fi2 = figure('Renderer', 'painters', 'Position', [199.6667 109 758.0000 587.3334]);
else 
   fi2 = figure('Renderer', 'painters', 'Position', [199.6667 73.6667 476.0000 622.6667]);
end
ti = tiledlayout(2,2,'TileSpacing','tight');
nexttile
imagesc(x,flipud(y),flipud(Tot_precip./perYear),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title(['Total precipitation [mm/yr]']); clim([0 2000]); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap);
xlabel(['Mean : ' num2str(round(nanmean(Tot_precip./perYear,'all'))) ' mm (' upper(glacier(1:3)) ' : ' ...
    num2str(round(nanmean(Tot_precip(GLA_ID == gla_id)./perYear))) ' mm)'],'FontSize',9)
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 

nexttile
imagesc(x,flipud(y),flipud(Snowfall_ratio),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title(['Snowfall ratio [-]']); clim([0 1]); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap);
xlabel(['Mean : ' num2str(round(nanmean(Snowfall_ratio,'all'),2)) ' (' upper(glacier(1:3)) ' : ' ...
    num2str(round(nanmean(Snowfall_ratio(GLA_ID == gla_id)),2)) ' )'],'FontSize',9)
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 

nexttile
imagesc(x,flipud(y),flipud(Tot_snowfall./perYear),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title(['Snowfall [mm/yr]']); clim([0 2000]); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap);
xlabel(['Mean : ' num2str(round(nanmean(Tot_snowfall./perYear,'all'))) ' mm (' upper(glacier(1:3)) ' : ' ...
    num2str(round(nanmean(Tot_snowfall(GLA_ID == gla_id)./perYear))) ' mm)'],'FontSize',9)
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
nexttile
imagesc(x,flipud(y),flipud(Tot_snowmelt./perYear),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Snowmelt [mm/yr]'); clim([0 1500]); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); colormap(gca, cmap_red);
xlabel(['Mean : ' num2str(round(nanmean(Tot_snowmelt./perYear,'all'))) ' mm (' upper(glacier(1:3)) ' : ' ...
    num2str(round(nanmean(Tot_snowmelt(GLA_ID == gla_id)./perYear))) ' mm)'],'FontSize',9)
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
title(ti,[PP_dataset ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))],'FontSize',11);
exportgraphics(fi2,[dir_fig '\TC_snowfall_melt_sublim_yearly_maps.png'],'Resolution',300,'BackgroundColor','none')
% close(gcf)

%% Sublimation and ET

if strcmp(glacier,'Langtang')  
   fi2 = figure('Renderer', 'painters', 'Position', [199.6667 109 758.0000 587.3334]);
else 
   fi2 = figure('Renderer', 'painters', 'Position',[199.6667 177.6667 714.6666 518.6667]);
end
ti = tiledlayout(1,2,'TileSpacing','compact');
nexttile
imagesc(x,flipud(y),flipud(Tot_ET./perYear),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('ET [mm]'); clim([0 600]); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);colormap(gca, cmap_red);
xlabel(['Mean : ' num2str(round(nanmean(Tot_ET./perYear,'all'))) ' mm (' upper(glacier(1:3)) ' : ' ...
    num2str(round(nanmean(Tot_ET(GLA_ID == gla_id)./perYear))) ' mm)'],'FontSize',9)
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 

nexttile
imagesc(x,flipud(y),flipud(Tot_sublim./perYear),'AlphaData',flipud(MASK)); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
set(gca,'YDir','normal');
title('Snow surface sublimation [mm]'); clim([0 300]); set(gca,'Color',[0.8 0.8 0.8])
colorbar; set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);colormap(gca, cmap_red);
xlabel(['Mean : ' num2str(round(nanmean(Tot_sublim./perYear,'all'))) ' mm (' upper(glacier(1:3)) ' : ' ...
    num2str(round(nanmean(Tot_sublim(GLA_ID == gla_id)./perYear))) ' mm)'],'FontSize',9)
axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 
title(ti,[PP_dataset ', ' num2str(Years_no(1)) '-' num2str(Years_no(end))],'FontSize',11);
exportgraphics(fi2,[dir_fig '\TC_annual_ET_sublimation_maps.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)