% Load MODIS maps and corresponding DEM

load([root '\Remote_sensing\MODIS\Processed_data\' glacier '\' glacier '_modis_40.mat'])
foSC = [root '\Remote_sensing\MODIS\Data\' glacier ' snow cover'];
fns = dir([foSC '\MOD10A1.006_NDSI_Snow_Cover_doy*.tif']);
imInfo = geotiffinfo([foSC '\' fns(1).name]);
imRM = imInfo.RefMatrix;
imSize = [imInfo.Height,imInfo.Width];
[imLons,imLats] = pixcenters(imRM,imSize);
[imLons,imLats] = meshgrid(imLons,imLats);
[imX,imY] = ll2utm(imLats,imLons, utm_zone);

fi2 = figure('Renderer', 'painters', 'Position', [185.6667 290.3333 750.6666 406.6667]);

for dm = 1:length(date_m)
    id_modis_i = find(year(TRA_modis.date) == year(date_m(dm)) & month(TRA_modis.date) == month(date_m(dm)) & TRA_modis.fData > 0.95,1);
    if ~isempty(id_modis_i)
        id_modis(dm) = id_modis_i; %% indices of MODIS images for which there is a corresponding T&C date
    else 
        id_modis(dm) = NaN;
    end
end 
 id_modis(isnan(id_modis)) = [];

for ii = 1:length(id_modis)

id_tc = find(datetime(year(date_h),month(date_h),day(date_h)) == TRA_modis.date(id_modis(ii)));

tiledlayout(1,2,'TileSpacing','compact');

s3 = nexttile;
imagesc(x,y,snow_depth(:,:,id_tc),'AlphaData',DTM > 0); %clim([0 1]); hold on;
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]); set(gca,'Color',[0.9 0.9 0.9])
title(['T&C with ' forcing_nm ' : Snow depth'])
set(gca,'YDir','normal'); set(gca,'YAxisLocation','Left')
colormap(s3,[1,1,1; c4])
cb = colorbar('Location','westoutside'); 
ylabel(cb,'Snow depth [m]','FontSize',11)
text(4.4771*10^5, 3.08909*10^6,datestr(Date_d(id_tc),'dd-mmm-yYYy'),'FontSize',11,'FontWeight','bold')
clim([0 1.5])
xlabel(['fsc = ' num2str(round(scas(id_tc),2)) ', SLA = ' num2str(round(SPAVG_dm.SLE(id_tc),-1)) ' m'] )

nexttile
imagesc(imX(1,:),imY(:,1),TRA_modis.ims(:,:,id_modis(ii)),'AlphaData',TRA_modis.mask & TRA_modis.ims_data(:,:,id_modis(ii)))
set(gca,'YDir','normal');  set(gca,'Color',[0.9 0.9 0.9])
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
text(4.4771*10^5, 3.08909*10^6,datestr(TRA_modis.date(id_modis(ii)),'dd-mmm-yYYy'),'FontSize',11,'FontWeight','bold')
xlabel(['fsc = ' num2str(round(MOD_sca_40.fSnow(id_modis(ii)),2)) ', SLA = ' num2str(round(MOD_sca_40.se(id_modis(ii)),-1)) ' m'] )

exportgraphics(fi2,[dir_fig '\Snow_cover\' glacier '_SnowCover_MODIS_validation_' datestr(TRA_modis.date(id_modis(ii)),'yYYy-mm-dd') ... 
    '.png'],'Resolution',300,'BackgroundColor','none')
end