% if contains(sim_nm, '2015_2020') || contains(sim_nm, '2015_2022') % Only those two cases so far 
% 
% path_SMB = [root '\Study_sites\' glacier '\Spatial_data\GMB\Miles2021\Trambau_2015-2019_01_mar_2022'];
% 
% SMB_miles = GRIDobj([path_SMB '\15.03448_SMB.tif']);
% DH_miles = GRIDobj([path_SMB '\15.03448_DH.tif']);
% DEM_miles = GRIDobj([path_SMB '\15.03448_DEM.tif']);
% DEM_miles_orig = DEM_miles;
% DEM_miles.Z(isnan(DH_miles.Z)) = NaN;
% 
% demSize = size(DEM_miles.Z);
% demInfo = geotiffinfo([path_SMB '\15.03448_DEM.tif']);
% demRM = demInfo.RefMatrix;
% utm_zone = demInfo.Zone;
% [SMB_x,SMB_y] = pixcenters(demRM,demSize,'makegrid');
% 
% dEL=50; % width of elevation bins
% ELs = nanmin(DEM_miles.Z,[],'all'):dEL:nanmax(DEM_miles.Z,[],'all');
%     
%     for iel = 1:numel(ELs)
%         cur=(DEM_miles.Z<(ELs(iel)+dEL/2))&(DEM_miles.Z>=(ELs(iel)-dEL/2)); %current section of DEM
%         pSMB_miles(iel)=nanmean(SMB_miles.Z(cur)); % mean elevation change within that elevation band
%         pSMB_std_miles(iel)=nanstd(SMB_miles.Z(cur)); % mean elevation change within that elevation band
%         pDH_miles(iel)=nanmean(DH_miles.Z(cur)); % mean elevation change within that elevation band
%     end 
% 
% end 

%% Check remotely-sensed SMB vs DH data 

% if contains(sim_nm, '2015_2020') || contains(sim_nm, '2015_2022') % Only those two cases so far 
% 
% fi2 = figure('Renderer', 'painters', 'Position', [119 217.6667 1.0227e+03 479.3333]);
% tiledlayout(1,3,'TileSpacing','compact');
% nexttile
% imagesc(SMB_miles.Z,'AlphaData',~isnan(SMB_miles.Z)); hold on;
% [C, h] = contour(DEM_miles_orig.Z,[2000:300:6400],'Color',[0.2 0.2 0.2],'EdgeAlpha',0.5);
% clabel(C,h,[2000:300:6400],'Color',[0.2 0.2 0.2],'LabelSpacing',450,'FontWeight','bold','FontSize',8)
% colormap(flipud(redblue)); 
% set(gca,'Color',[0.8 0.8 0.8]); clim([-3 3])
% cb = colorbar('Location','Westoutside'); ylabel(cb,'SMB [m w.e.yr^{-1}]','FontSize',11)
% set(gca,'XTickLabel',[],'YTickLabel',[])
% 
% nexttile
% plot(pSMB_miles,ELs,'LineWidth',0.9); grid on; hold on;
% plot(pDH_miles,ELs,'LineWidth',0.9); grid on; hold on;
% plot(pDH_hug,ELs_hug,'LineWidth',0.9); grid on; hold on;
% xline(0,'--k','LineWidth',0.7)
% legend('SMB [m w.e. yr^{-1}]','DH [m.yr^{-1}]','DH Hugonnet','Location','NorthWest')
% xlim([-3 1])
% title({[glacier ' - 2015-2019 SMB/DH (Miles et al. 2021)']},{' '})
% 
% nexttile
% imagesc(DH_miles.Z,'AlphaData',~isnan(DH_miles.Z)); hold on;
% [C, h] = contour(DEM_miles_orig.Z,[2000:300:6400],'Color',[0.2 0.2 0.2],'EdgeAlpha',0.5);
% clabel(C,h,[2000:300:6400],'Color',[0.2 0.2 0.2],'LabelSpacing',450,'FontWeight','bold','FontSize',8)
% colormap(flipud(redblue)); set(gca,'Color',[0.8 0.8 0.8]); clim([-3 3])
% cb = colorbar; ylabel(cb,'DH [m.yr^{-1}]','FontSize',11)
% set(gca,'XTickLabel',[],'YTickLabel',[])
% 
% exportgraphics(fi2,[dir_fig '\SMB_Miles_2021.png'],'Resolution',300,'BackgroundColor','none')
% end 


%% GMB comparison with Miles SMB
% if year(Date(1)) <= 2015 && year(Date(end)) >= 2019
% 
% fi2 = figure('Renderer', 'painters', 'Position', [119 217.6667 1.0227e+03 479.3333]);
% tiledlayout(1,3,'TileSpacing','compact');
% nexttile
% imagesc(x,y,flipud(GMB_Hug_Trambau_tc.*0.001),'AlphaData',~isnan(flipud(GMB_Hug_Trambau_tc))); hold on;
% [C, h] = contour(x,y,DEM.Z,[2000:300:6400],'Color',[0.2 0.2 0.2],'EdgeAlpha',0.5);
% clabel(C,h,[2000:300:6400],'Color',[0.2 0.2 0.2],'LabelSpacing',450,'FontWeight','bold','FontSize',8)
% colormap(flipud(redblue)); xlim([x(1) x(end)]); ylim([y(1) y(end)])
% set(gca,'Color',[0.8 0.8 0.8]); clim([-3 3])
% cb = colorbar('Location','Westoutside'); ylabel(cb,'SMB [m w.e.yr^{-1}]','FontSize',11)
% set(gca,'XTickLabel',[],'YTickLabel',[])
% xlabel({'Modelled',['Mean SMB = ' num2str(round(nanmean(GMB_Hug_Trambau_tc.*0.001,'all'),2))]})
% 
% nexttile
% imagesc(SMB_x(1,:), flipud(SMB_y(:,1)),SMB_miles.Z,'AlphaData',~isnan(SMB_miles.Z)); hold on;
% [C, h] = contour(SMB_x(1,:),flipud(SMB_y(:,1)),DEM_miles_orig.Z,[2000:300:6400],'Color',[0.2 0.2 0.2],'EdgeAlpha',0.5);
% clabel(C,h,[2000:300:6400],'Color',[0.2 0.2 0.2],'LabelSpacing',450,'FontWeight','bold','FontSize',8)
% colormap(flipud(redblue)); xlim([x(1) x(end)]); ylim([y(1) y(end)])
% set(gca,'Color',[0.8 0.8 0.8]); clim([-3 3])
% set(gca,'XTickLabel',[],'YTickLabel',[])
% title({[glacier ' - 2015-2019 SMB']},{' '})
% xlabel({'Miles et al. 2021',['Mean SMB = ' num2str(round(nanmean(SMB_miles.Z,'all'),2))]})
% 
% nexttile
% shadedErrorBar(double(ELs),double(pSMB_miles),double(pSMB_std_miles),'lineProps',{'k','LineWidth',1.1}); grid on; hold on;
% shadedErrorBar(double(ELs_tc),double(pSMB_tc).*0.001,double(pSMB_std_tc).*0.001,'lineProps',{'r','LineWidth',1.1})
% yline(0,'--k','LineWidth',0.7); ylabel(['SMB ' char(177) ' 1std [m w.e. yr^{-1}]']);
% xlabel('Elevation [m a.s.l.]')
% legend('Miles et al. 2021','Modelled','Location','NorthWest')
%  set(gca,'XAxisLocation','top')
% ylim([-7 3]); view([90, -90])
% exportgraphics(fi2,[dir_fig '\GMB\' glacier '_2015-2019_SMB_validation.png'],'Resolution',300,'BackgroundColor','none')
% 
% end
