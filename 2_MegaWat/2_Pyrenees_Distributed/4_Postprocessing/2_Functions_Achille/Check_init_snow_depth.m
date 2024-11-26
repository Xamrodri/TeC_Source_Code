%Needs:
% -Datestamp of the first timestep (startDate(1))
% -path to Landsat scenes (path_L8S2)
% -daily snow depth from T&C output (snow_depth)
% -DTM

if ~exist('startDate','var')
    startDate = Date(1);
end 

% Load closest in time Sentinel-Landsat snow cover

[d,ix] = min(abs(datenum(year(startDate(1)),month(startDate(1)),day(startDate(1)))-datenum(date_l8s2)));

l8s2_rgb_load=  imread([path_L8S2 '\TIF\' fn_maps(ix).name]);
l8s2_r_ini = l8s2_rgb_load(:,:,1);
l8s2_g_ini = l8s2_rgb_load(:,:,2);
l8s2_b_ini = l8s2_rgb_load(:,:,3);

l8s2_sno_load=  imread([path_L8S2 '\TIF\' fn_snowmaps(ix).name]);
l8s2_sno_ini = l8s2_sno_load; 

iniSND_show = flipud(snow_depth(:,:,1));
iniSND_show(flipud(DTM) <= 0) = NaN;
c4 = colormap(getPyPlot_cMap('GnBu')); c4 = [1 1 1; c4];

l8s2_rgb = cat(3,l8s2_r_ini, l8s2_g_ini, l8s2_b_ini);
l8s2_rgb = l8s2_rgb + 0.2;
l8s2_rgb(l8s2_rgb>1) = 1;

% Figure

fi3 = figure('Renderer', 'painters', 'Position', [360 263.6667 568.3333 434.3333]) ;
tiledlayout(1,2,'TileSpacing','compact');
nexttile
imagesc(x_l8s2(1,:),y_l8s2(:,1),l8s2_rgb,'AlphaData',l8s2_r_ini ~= -0.9999); clim([0 1]); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.9 0.9 0.9])
xlabel([datestr(date_l8s2(ix),'dd-mmm-yYYy')])
plot(catchShp_utm.X,catchShp_utm.Y,'Color',[0 0 0 0.9])
title([extractBefore(fn_snowmaps(ix).name,'_Snow') ' - RGB'])
set(gca,'YDir','normal'); axis image
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 

nexttile
imageschs(DEM,iniSND_show,'colormap',c4,'colorbar',1); hold on;
cb =  colorbar;  ylabel(cb,'Snow depth [m]','FontSize',12);
set(gca,'XTickLabels',[],'YTickLabels',[]); set(gca,'Color',[0.9 0.9 0.9])
plot(catchShp_utm.X,catchShp_utm.Y,'Color',[0 0 0 0.9])
title('Model initial conditions'); axis image
set(gca,'YDir','normal'); clim([0 round(prctile(snow_depth(:,:,1),95,'all'),1)]); xlabel(datestr(startDate,'dd-mmm-yyyy'))
xlim([min(demLons(~isnan(flipud(DTM))))-100 max(demLons(~isnan(flipud(DTM))))+100]); 
ylim([min(demLats(~isnan(flipud(DTM))))-100 max(demLats(~isnan(flipud(DTM))))+100]); 

exportgraphics(fi3,[dir_fig '\Initial_snow_depth.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)