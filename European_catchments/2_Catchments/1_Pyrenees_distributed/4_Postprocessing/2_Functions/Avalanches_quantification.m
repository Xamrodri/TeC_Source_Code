if ~exist('gla_id_ava','var')
    gla_id_ava = gla_id; % In case no glaciers IDs other than the main glaciers were provided
end 

%%% Total over the whole period

Aval_period = nansum(AVA_map,3);
Aval_Main = Aval_period;
Aval_Main(GLA_ID ~= gla_id) = NaN;
aval_Main = nanmean(Aval_Main,'all')./(size(AVA_map,3)/12); % -0.56m w.e./yr for 2015-2020 ERA5-Land

if nansum(Aval_period(:)) ~= 0  % No need for avalanche plots if there are no avalanched

    fi3 = figure('Renderer', 'painters', 'Position',[360 208.3333 808.3333 489.6667]);
    tiledlayout(1,2,"TileSpacing","compact")
    nexttile
    imagesc(demLons(1,:), demLats(:,1),flipud(Aval_period./(size(AVA_map,3)/12)), 'AlphaData', flipud(MASK))
%     clim([-round(abs(prctile(Aval_Main(:)./(size(AVA_map,3)/12),2)),-1) abs(round(prctile(Aval_Main(:)./(size(AVA_map,3)/12),2),-1))])
    clim([-2000 2000])
    hold on; colormap(flipud(redblue))
    plot(catchShp_utm.X, catchShp_utm.Y,'k')
    cb1 = colorbar; ylabel(cb1,'Snow avalanched [mm/yr]','FontSize',12)
    for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
    set(gca,'YDir','normal'); set(gca,'Color',[0.8 0.8 0.8])
    set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); axis image; 
    xlabel([glacier ' avalanches supply: ' num2str(round(aval_Main)) ' mm w.e./yr'],'FontSize',10)
    title(['T&C: ' datestr(date_m(end),'dd-mmm-yyyy') ' - ' datestr(date_m(1),'dd-mmm-yyyy')]);
    n2 = nexttile;
    imagesc(demLons(1,:), demLats(:,1),flipud(snow_depth(:,:,end)), 'AlphaData', flipud(MASK))
    hold on; colormap(n2, cbrewer('seq','BuPu',10))
    plot(catchShp_utm.X, catchShp_utm.Y,'k'); ; clim([0 10])
    cb1 = colorbar; ylabel(cb1,'Final snow depth [m w.e.]','FontSize',12)
    for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
    title(datestr(date_m(end)))
    set(gca,'YDir','normal'); set(gca,'Color',[0.8 0.8 0.8])
    set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); axis image

    exportgraphics(fi3,[dir_fig '\Avalanching\' glacier '_Avalanches_'  datestr(date_m(end),'yyyy-mm-dd') '-' datestr(date_m(1),'yyyy-mm-dd')...
    '.png'],'Resolution',300,'BackgroundColor','none')

    AVAL_obj = DEM;
    AVAL_obj.Z = flipud(nansum(AVA_map,3));

    GRIDobj2geotiff(AVAL_obj,[dir_fig '\Avalanching\' glacier '_Avalanches_'  datestr(date_m(end),'yyyy-mm-dd') '-' datestr(date_m(1),'yyyy-mm-dd') '.tif'])

end 

%% For avalanche years

Hydro_years = (year(date_m(find(month(date_m) == 10 & day(date_m) == 1,1))):...
    1:year(date_m(find(month(date_m) == 9,1,'last'))))';

Hydro_year_label = strcat(string(num2str(Hydro_years(1:end-1))), "-", ...
   string(num2str(Hydro_years(2:end))));

if length(Hydro_years) > 1

for yy = 1:length(Hydro_years)-1

    ind_start = find(year(date_m) == Hydro_years(yy) & month(date_m) == 10 & day(date_m) == 1,1);
    ind_end =   find(year(date_m) == Hydro_years(yy+1) & month(date_m) == 9 & day(date_m) == 1,1);
    ind_period = ind_start:ind_end;

    Ava_hy_map = nansum(AVA_map(:,:,ind_start:ind_end),3);
    Ava_hy_map(GLA_ID ~= gla_id_ava) = NaN; 
    Ava_hy(yy) = nanmean(Ava_hy_map,'all');

    Ava_hy_map_all = nansum(AVA_map(:,:,ind_start:ind_end),3);
    Ava_hy_map_all(MASK ~= 1) = NaN; 

if nansum(Ava_hy(yy)) ~= 0  % No need for avalanche plots if there are no avalanched

    fi3 = figure('Renderer', 'painters', 'Position', [360 273 474.3333 425]);
    imagesc(demLons(1,:), demLats(:,1),flipud(Ava_hy_map_all), 'AlphaData', flipud(MASK))
    clim([-3000 3000])
    hold on; colormap(flipud(redblue))
    plot(catchShp_utm.X, catchShp_utm.Y,'k')
    cb1 = colorbar; ylabel(cb1,'Snow avalanched [mm/yr]','FontSize',12)
    for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
    set(gca,'YDir','normal'); set(gca,'Color',[0.8 0.8 0.8])
    set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
    xlabel([glacier ' avalanches supply: ' num2str(round(Ava_hy(yy))) ' mm w.e./yr'],'FontSize',10)
    title(['T&C, hydro-year: ' num2str(Hydro_year_label(yy)) ]); axis image
    exportgraphics(fi3,[dir_fig '\Avalanching\' glacier '_Avalanches_hydroYear_' num2str(Hydro_year_label(yy)) ...
    '.png'],'Resolution',300,'BackgroundColor','none')
    close(gcf)
end 

end 
end 

clear gla_id_ava
