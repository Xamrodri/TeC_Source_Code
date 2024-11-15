% Create heatmaps, from the monthly variables calculated in Seasonal_runoff

% Where to store them:
dir_fig_heatmap = [dir_fig '\Heatmaps'];
if ~exist(dir_fig_heatmap, 'dir'); mkdir(dir_fig_heatmap); end 

%% CAtchment air temperature and precipitation anomaly

ISO_0_ym(isnan(ISO_0_ym)) = min(DTM(:));

cmap_precip = cbrewer('seq', 'BuPu', 80); cmap_precip(cmap_precip<0)=0; cmap_precip(cmap_precip>1)=1;
cmap_precip_anom = cbrewer('div', 'BrBG', 80); cmap_precip_anom(cmap_precip_anom<0)=0; cmap_precip_anom(cmap_precip_anom>1)=1;

fi5 = figure('Renderer', 'painters', 'Position',[289 235.6667 727.3333 423.3333]);
tiledlayout(1,2,"TileSpacing","compact")
nexttile
Ta_yearly_anom = nanmean(Ta_yearly_ym,2) - nanmean(nanmean(Ta_yearly_ym,2));
imagesc([Ta_yearly_ym - nanmean(Ta_yearly_ym,1) Ta_yearly_anom]); clim([-4 4]); colormap(redblue);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;
cb1 = colorbar; ylabel(cb1,'Air temp anomaly [°C]','FontSize',11);  xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1)
n2 = nexttile;
Pr_yearly_anom = (nansum(Pr_ym,2) - nanmean(nansum(Pr_ym,2)))./nanmean(nansum(Pr_ym,2));
Pr_spring_anom = (nansum(Pr_ym(:,5:9),2) - nanmean(nansum(Pr_ym(:,5:9),2)))./nanmean(nansum(Pr_ym(:,5:9),2));
imagesc([100*(Pr_ym - nanmean(Pr_ym,1))./nanmean(Pr_ym,1) 100*Pr_yearly_anom]); colormap(n2,colorbrewer.div.BrBG{1, 11}./255);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;
cb2 = colorbar; ylabel(cb2,'Precipitation anomaly [%]','FontSize',11); xline(size(Pr_ym,2)+0.5,'k','LineWidth',1.1)
clim([-200 200]);
exportgraphics(fi5,[dir_fig_heatmap '\Catchment_monthly_Ta_Precip_anomaly_HydroYears.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

Pr_liq_ym = Pr_ym - Pr_sno_ym;

%% Glacier monthly snowfall anomaly

fi5 = figure('Renderer', 'painters', 'Position',[289 235.6667 727.3333 423.3333]);
tiledlayout(1,2,"TileSpacing","compact")
n1 = nexttile;
imagesc(Pr_sno_ym_gla); colormap(n1,cmap_precip);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); 
cb2 = colorbar; ylabel(cb2,'Snowfall [mm w.e.]','FontSize',11); 
clim([0 350]);
n2 = nexttile;
Pr_gla_sno_anom_y = (nansum(Pr_sno_ym_gla,2) - nanmean(nansum(Pr_sno_ym_gla,2)))./nanmean(nansum(Pr_sno_ym_gla,2));
imagesc([100*(Pr_sno_ym_gla - nanmean(Pr_sno_ym_gla,1))./nanmean(Pr_sno_ym_gla,1) 100*Pr_gla_sno_anom_y]); colormap(n2,cmap_precip_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'Annual']);hold on;
cb2 = colorbar; ylabel(cb2,'Snowfall anomaly [%]','FontSize',11); xline(size(Pr_sno_ym_gla,2)+0.5,'k','LineWidth',1.1)
clim([-150 150]);
exportgraphics(fi5,[dir_fig_heatmap '\Glacier_Monthly_Precip_anomaly_HydroYears.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Catchment monthly snowfall anomaly

fi5 = figure('Renderer', 'painters', 'Position',[289 235.6667 727.3333 423.3333]);
tiledlayout(1,2,"TileSpacing","compact")
n1 = nexttile;
imagesc(Pr_sno_ym); colormap(n1,cmap_precip);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); 
cb2 = colorbar; ylabel(cb2,'Snowfall [mm w.e.]','FontSize',11); 
% clim([-150 150]);
n2 = nexttile;
Pr_sno_anom_y = (nansum(Pr_sno_ym,2) - nanmean(nansum(Pr_sno_ym,2)))./nanmean(nansum(Pr_sno_ym,2));
imagesc([100*(Pr_sno_ym - nanmean(Pr_sno_ym,1))./nanmean(Pr_sno_ym,1) 100*Pr_sno_anom_y]); colormap(n2,cmap_precip_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); cb2 = colorbar; ylabel(cb2,'Snowfall anomaly [%]','FontSize',11); xline(size(Pr_sno_ym,2)+0.5,'k','LineWidth',1.1)
clim([-150 150]);
exportgraphics(fi5,[dir_fig_heatmap '\Catchment_Monthly_Precip_anomaly_HydroYears.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Monthly glacier mass balance

cmap_smb = cbrewer('div','PuOr',100); cmap_smb(cmap_smb<0)=0; cmap_smb(cmap_smb>1)=1;

fi5 = figure('Renderer', 'painters', 'Position',[289 235.6667 727.3333 423.3333]);
tiledlayout(1,2,"TileSpacing","compact")
n1 = nexttile;
imagesc([SMB_ym.*0.001, nansum(SMB_ym,2).*0.001] ); clim([-1 1]); colormap(n1,flipud(redblue));
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); cb2 = colorbar; ylabel(cb2,'Snowfall anomaly [%]','FontSize',11); xline(size(Pr_sno_ym,2)+0.5,'k','LineWidth',1.1)
cb1 = colorbar; ylabel(cb1,'Glacier SMB [m w.e.]','FontSize',11); xline(size(SMB_ym,2)+0.5,'k','LineWidth',1.1)
n2 = nexttile;
imagesc([SMB_ym.*0.001 - nanmean(SMB_ym.*0.001,1), 0.001*nansum(SMB_ym - nanmean(SMB_ym,1),2)]); 
clim([-0.5 0.5]); colormap(n2,cmap_smb);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); cb2 = colorbar; ylabel(cb2,'Snowfall anomaly [%]','FontSize',11); xline(size(Pr_sno_ym,2)+0.5,'k','LineWidth',1.1)
cb1 = colorbar; ylabel(cb1,'Glacier SMB anomaly [m w.e.]','FontSize',11); xline(size(SMB_ym,2)+0.5,'k','LineWidth',1.1)
exportgraphics(fi5,[dir_fig_heatmap '\Glacier_SMB_monthly_anomaly_HydroYears.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

% Numbers for the paper

find_Jul22 = find(year(date_ym)==2022 & month(date_ym)==7);

disp(['Kyzylsu mass balance July 2022: ' num2str(SMB_ym(find_Jul22)*0.001,2) ' m w.e.']);

disp(['Kyzylsu mass balance for July (2000-2023): ' num2str(nanmean(SMB_ym(month(date_ym)==7))*0.001,2) ' ' char(177) ' ' ...
    num2str(nanstd(SMB_ym(month(date_ym)==7))*0.001,2) ' m w.e.']);


%% Glacier monthly icemelt, snowfall and SMB anomaly

fi5 = figure('Renderer', 'painters', 'Position',[19.6667 279.6667 1.0007e+03 329.3333]);
tiledlayout(1,3,"TileSpacing","compact")
n1 = nexttile;
Pr_gla_sno_anom_y = nansum(Pr_sno_ym_gla,2) - nanmean(nansum(Pr_sno_ym_gla,2));
imagesc([0.001*(Pr_sno_ym_gla-nanmean(Pr_sno_ym_gla,1)) 0.001*Pr_gla_sno_anom_y]); colormap(n1,cmap_precip_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'Annual']);hold on;
cb2 = colorbar; ylabel(cb2,'On-Glacier Snowfall anomaly [m w.e.]','FontSize',11); xline(size(Pr_sno_ym_gla,2)+0.5,'k','LineWidth',1.1)
clim([-0.5 0.5]);
n2 = nexttile;
Imelt_gla_anom_y = nansum(Imelt_ym_gla,2) - nanmean(nansum(Imelt_ym_gla,2));
imagesc([0.001*(Imelt_ym_gla-nanmean(Imelt_ym_gla,1)) 0.001*Imelt_gla_anom_y]); colormap(n2,cmap_precip_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'Annual']);hold on;
cb2 = colorbar; ylabel(cb2,'Icemelt anomaly [m w.e.]','FontSize',11); xline(size(Imelt_ym_gla,2)+0.5,'k','LineWidth',1.1)
clim([-0.5 0.5]);
n3 = nexttile;
imagesc([SMB_ym.*0.001 - nanmean(SMB_ym.*0.001,1), 0.001*nansum(SMB_ym - nanmean(SMB_ym,1),2)]); 
clim([-0.5 0.5]); colormap(n3,cmap_smb);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); cb2 = colorbar; ylabel(cb2,'Snowfall anomaly [%]','FontSize',11); xline(size(Pr_sno_ym,2)+0.5,'k','LineWidth',1.1)
cb1 = colorbar; ylabel(cb1,'Glacier SMB anomaly [m w.e.]','FontSize',11); xline(size(SMB_ym,2)+0.5,'k','LineWidth',1.1)
exportgraphics(fi5,[dir_fig_heatmap '\Glacier_SMB_Snowfall_Icemelt_monthly_anomaly_HydroYears.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Heat map of MODIS snow cover and snow cover anomalies

cmap_fsc = cbrewer('seq', 'Greys', 80); cmap_fsc(cmap_fsc<0)=0; cmap_fsc(cmap_fsc>1)=1;
cmap_fsc_anom = cbrewer('div', 'BrBG', 80); cmap_fsc_anom(cmap_fsc_anom<0)=0; cmap_fsc_anom(cmap_fsc_anom>1)=1;

fi5 = figure('Renderer', 'painters', 'Position',[289 235.6667 727.3333 423.3333]);
tiledlayout(1,2,"TileSpacing","compact")
nexttile
imagesc([fsc_MODIS_ym nanmean(fsc_MODIS_ym,2)],'AlphaData',~isnan([fsc_MODIS_ym nanmean(fsc_MODIS_ym,2)])); clim([0 1]); colormap(flipud(cmap_fsc));
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on; 
set(gca,'Color',[1.0000    0.8118    0.9608])
cb1 = colorbar; ylabel(cb1,'Snow cover fraction [-]','FontSize',11);  xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1)
n2 = nexttile;
MODIS_yearly_anom = nanmean(fsc_MODIS_ym,2) - nanmean(nanmean(fsc_MODIS_ym,2));
imagesc([fsc_MODIS_ym - nanmean(fsc_MODIS_ym,1) MODIS_yearly_anom],'AlphaData',~isnan([fsc_MODIS_ym MODIS_yearly_anom])); clim([-0.10 0.10]); colormap(n2,cmap_fsc_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;  set(gca,'Color',[1.0000    0.8118    0.9608])
cb1 = colorbar; ylabel(cb1,'Snow cover fraction anomaly [-]','FontSize',11);  xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1)
exportgraphics(fi5,[dir_fig_heatmap '\MODIS_snow_cover_anomaly_HydroYears.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Heat map of MODIS snow cover and snow cover anomalies

cmap_fsc = cbrewer('seq', 'Greys', 80); cmap_fsc(cmap_fsc<0)=0; cmap_fsc(cmap_fsc>1)=1;
cmap_fsc_anom = cbrewer('div', 'BrBG', 80); cmap_fsc_anom(cmap_fsc_anom<0)=0; cmap_fsc_anom(cmap_fsc_anom>1)=1;

fi5 = figure('Renderer', 'painters', 'Position',[103 69 846 641.3333]);
tiledlayout(2,2,"TileSpacing","compact")
nexttile
imagesc([fsc_MODIS_ym nanmean(fsc_MODIS_ym,2)],'AlphaData',~isnan([fsc_MODIS_ym nanmean(fsc_MODIS_ym,2)])); clim([0 1]); colormap(flipud(cmap_fsc));
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on; 
set(gca,'Color',[1.0000    0.8118    0.9608]);
% xlabel('MODIS','FontSize',11,'FontWeight','bold'); 
set(gca,'XAxisLocation','top'); title('MODIS')
%cb1 = colorbar; ylabel(cb1,'Snow cover fraction [-]','FontSize',11);  
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1)

n3 = nexttile;
imagesc([fsc_tc_ym nanmean(fsc_tc_ym,2)],'AlphaData',~isnan([fsc_tc_ym nanmean(fsc_tc_ym,2)])); clim([0 1]); colormap(n3,flipud(cmap_fsc));
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on; 
set(gca,'Color',[1.0000    0.8118    0.9608])
cb1 = colorbar; ylabel(cb1,'Snow cover fraction [-]','FontSize',11);  
xlabel('T&C','FontSize',11,'FontWeight','bold'); set(gca,'XAxisLocation','top')
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1)

n2 = nexttile;
MODIS_yearly_anom = nanmean(fsc_MODIS_ym,2) - nanmean(nanmean(fsc_MODIS_ym,2));
imagesc([fsc_MODIS_ym - nanmean(fsc_MODIS_ym,1) MODIS_yearly_anom],'AlphaData',~isnan([fsc_MODIS_ym MODIS_yearly_anom])); clim([-0.20 0.20]); colormap(n2,cmap_fsc_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;  set(gca,'Color',[1.0000    0.8118    0.9608])
%cb1 = colorbar; ylabel(cb1,'Snow cover fraction anomaly [-]','FontSize',11);  
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1)

n4 = nexttile;
TCfsc_yearly_anom = nanmean(fsc_tc_ym,2) - nanmean(nanmean(fsc_tc_ym,2));
imagesc([fsc_tc_ym - nanmean(fsc_tc_ym,1) TCfsc_yearly_anom],'AlphaData',~isnan([fsc_tc_ym TCfsc_yearly_anom])); clim([-0.2 0.2]); colormap(n4,cmap_fsc_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;  set(gca,'Color',[1.0000    0.8118    0.9608])
cb1 = colorbar; ylabel(cb1,'Snow cover fraction anomaly [-]','FontSize',11);  
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1)

exportgraphics(fi5,[dir_fig_heatmap '\TC_MODIS_snow_cover_anomaly_HydroYears.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Heat-map of 0°C isotherm

cmap_iso = cbrewer('div', 'Spectral', 80); cmap_iso(cmap_iso<0)=0; cmap_iso(cmap_iso>1)=1;
cmap_iso_anom = cbrewer('div', 'RdBu', 80); cmap_iso_anom(cmap_iso_anom<0)=0; cmap_iso_anom(cmap_iso_anom>1)=1;


fi5 = figure('Renderer', 'painters', 'Position',[289 235.6667 727.3333 423.3333]);
tiledlayout(1,2,"TileSpacing","compact")
n1 = nexttile;
imagesc(ISO_0_ym); colormap(n1,flipud(cmap_iso));
set(gca,'YTick',1:2:length(ISO_0_ym)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); 
cb2 = colorbar; ylabel(cb2,'Isotherm 0°C [m a.s.l.]','FontSize',12); 
% clim([-150 150]);
n2 = nexttile;
imagesc(ISO_0_ym - nanmean(ISO_0_ym,1)); colormap(n2,flipud(cmap_iso_anom));
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); cb2 = colorbar; 
ylabel(cb2,'Isotherm 0°C anomaly [m]','FontSize',12); xline(size(Pr_sno_ym,2)+0.5,'k','LineWidth',1.1)
clim([-500 500]);
exportgraphics(fi5,[dir_fig_heatmap '\Catchment_Isotherm0_anomaly_HydroYears.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Kyzylsu paper, figure on snowfall heat map + snow cover anomalies

fi5 = figure('Renderer', 'painters', 'Position',[42.3333 195.6667 1.1574e+03 456.6667]);
tiledlayout(1,4,"TileSpacing","compact")
n1 = nexttile;
imagesc(Pr_sno_ym); colormap(n1,cmap_precip);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:12); set(gca,'XTickLabel',hydro_month_labels(1:2:12)); 
cb2 = colorbar('Location','northoutside'); ylabel(cb2,'Snowfall [mm w.e.]','FontSize',11); 
yline(find(Years_no_seas==2018)-0.5,'r','LineWidth',1.2); ylabel('Hydrological years','FontSize',11); clim([0 250])

n2 = nexttile;
Pr_sno_anom_y = (nansum(Pr_sno_ym,2) - nanmean(nansum(Pr_sno_ym,2)))./nanmean(nansum(Pr_sno_ym,2));
imagesc([100*(Pr_sno_ym - nanmean(Pr_sno_ym,1))./nanmean(Pr_sno_ym,1) 100*Pr_sno_anom_y]); colormap(n2,cmap_precip_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']); 
set(gca,'YTickLabels',[])
cb2 = colorbar('Location','northoutside'); ylabel(cb2,'Snowfall anomaly [%]','FontSize',11); xline(size(Pr_sno_ym,2)+0.5,'k','LineWidth',1.1)
clim([-150 150]); yline(find(Years_no_seas==2018)-0.5,'r','LineWidth',1.2)

n3 = nexttile;
TCfsc_yearly_anom = nanmean(fsc_tc_ym,2) - nanmean(nanmean(fsc_tc_ym,2));
imagesc([fsc_tc_ym - nanmean(fsc_tc_ym,1) TCfsc_yearly_anom],'AlphaData',~isnan([fsc_tc_ym TCfsc_yearly_anom])); 
clim([-0.3 0.3]); colormap(n3,cmap_fsc_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;  set(gca,'Color',[1.0000    0.8118    0.9608])
cb1 = colorbar('Location','northoutside'); ylabel(cb1,{'Snow persistence anomaly [-]','Simulated'},'FontSize',11);  
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1); set(gca,'YTickLabels',[])
yline(find(Years_no_seas==2018)-0.5,'r','LineWidth',1.2)

n4 = nexttile;
MODIS_yearly_anom = nanmean(fsc_MODIS_ym,2) - nanmean(nanmean(fsc_MODIS_ym,2));
imagesc([fsc_MODIS_ym - nanmean(fsc_MODIS_ym,1) MODIS_yearly_anom],'AlphaData',~isnan([fsc_MODIS_ym MODIS_yearly_anom])); 
clim([-0.3 0.3]); colormap(n4,cmap_fsc_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;  
set(gca,'Color',[1.0000    0.8118    0.9608]); set(gca,'YAxisLocation','right')
cb3 = colorbar('Location','northoutside'); ylabel(cb3,{'Snow persistence anomaly [-]','Observed (MODIS)'},'FontSize',11);  
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1); yline(find(Years_no_seas==2018)-0.5,'r','LineWidth',1.2)

exportgraphics(fi5,[dir_fig_heatmap '\Catchment_Monthly_Snowfall_Snowcover_anomaly_HydroYears_brown.png'],'Resolution',300,'BackgroundColor','none')
% close(gcf)
%% MODIS vs T&C monthly fsc anomalies scatter

fsc_tc_anom_ym = fsc_tc_ym - nanmean(fsc_tc_ym,1);
fsc_MODIS_anom_ym = fsc_MODIS_ym - nanmean(fsc_MODIS_ym,1);

fsc_MODIS_anom_ym = fsc_MODIS_anom_ym(:,5:12);
fsc_tc_anom_ym = fsc_tc_anom_ym(:,5:12);

mean_fsc_anom_comp = nanmean(fsc_tc_anom_ym(:)-fsc_MODIS_anom_ym(:));
rmse_fsc_anom_comp = rmse(fsc_tc_anom_ym(:),fsc_MODIS_anom_ym(:),"omitnan");
r_fsc_anom = fitlm(fsc_tc_anom_ym(:),fsc_MODIS_anom_ym(:)).Rsquared.Adjusted;

fi5 = figure('Renderer', 'painters', 'Position',[289 291 392.6667 368]);
scatter(fsc_MODIS_anom_ym(:),fsc_tc_anom_ym(:),20,'filled')
xlim([-0.3 0.3]); ylim([-0.3 0.3]); grid on; hold on;
plot([-0.5 0.5],[-0.5 0.5],'--k');
annotation('textbox',[0.18 0.77 0.1 0.1],'string',{['ME: ' num2str(round(mean_fsc_anom_comp,3))], ['RMSE: ' num2str(round(rmse_fsc_anom_comp,3)) ], ...
   ['R^2 =' num2str(round(r_fsc_anom,3))]},'FontSize',8,'BackgroundColor',[1 1 1])
axis square
ylabel('Modelled monthly FSC anomaly [-]')
xlabel('MODIS monthly FSC anomaly [-]')
title({'2000-2023 Feb-Sep snow cover anomaly', ' '})
exportgraphics(fi5,[dir_fig_heatmap '\Snowcover_anomaly_MODIS_comp_scatter.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%%
fi5 = figure('Renderer', 'painters', 'Position',[42.3333 276.3333 996.6667 376.0001]);
tiledlayout(1,3,"TileSpacing","compact")

n3 = nexttile;
TCfsc_yearly_anom = nanmean(fsc_tc_ym,2) - nanmean(nanmean(fsc_tc_ym,2));
imagesc([fsc_tc_ym - nanmean(fsc_tc_ym,1) TCfsc_yearly_anom],'AlphaData',~isnan([fsc_tc_ym TCfsc_yearly_anom])); clim([-0.2 0.2]); colormap(n3,cmap_fsc_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;  set(gca,'Color',[1.0000    0.8118    0.9608])
cb1 = colorbar('Location','northoutside'); ylabel(cb1,'Modelled snow persistence anomaly [-]','FontSize',11);  
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1); %set(gca,'YTickLabels',[])
yline(19.5,'r','LineWidth',1.2)

n4 = nexttile;
MODIS_yearly_anom = nanmean(fsc_MODIS_ym,2) - nanmean(nanmean(fsc_MODIS_ym,2));
imagesc([fsc_MODIS_ym - nanmean(fsc_MODIS_ym,1) MODIS_yearly_anom],'AlphaData',~isnan([fsc_MODIS_ym MODIS_yearly_anom])); 
clim([-0.20 0.20]); colormap(n4,cmap_fsc_anom);
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;  
set(gca,'Color',[1.0000    0.8118    0.9608]); set(gca,'YAxisLocation','right')
cb3 = colorbar('Location','northoutside'); ylabel(cb3,'MODIS snow persistence anomaly [-]','FontSize',11);  set(gca,'YTickLabels',[])
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1); yline(19.5,'r','LineWidth',1.2)

n5 = nexttile;
imagesc((fsc_MODIS_ym - nanmean(fsc_MODIS_ym,1)) - (fsc_tc_ym - nanmean(fsc_tc_ym,1)),'AlphaData',~isnan(fsc_MODIS_ym)); 
clim([-0.20 0.20]); colormap(n5,flipud(redblue));
set(gca,'YTick',1:2:length(Years_no_seas)); set(gca,'YTickLabel',Hydro_year_label(1:2:end)); 
set(gca,'XTick',1:2:13); set(gca,'XTickLabel',[hydro_month_labels(1:2:12) 'H-year']);hold on;  
set(gca,'Color',[1.0000    0.8118    0.9608]); set(gca,'YAxisLocation','right')
cb5 = colorbar('Location','northoutside'); ylabel(cb5,'Modeled - MODIS persistence anomaly [-]','FontSize',11);  
xline(size(Ta_yearly_ym,2)+0.5,'k','LineWidth',1.1); yline(19.5,'r','LineWidth',1.2)
exportgraphics(fi5,[dir_fig_heatmap '\Snowcover_anomaly_MODIS_comp_heatmap.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)