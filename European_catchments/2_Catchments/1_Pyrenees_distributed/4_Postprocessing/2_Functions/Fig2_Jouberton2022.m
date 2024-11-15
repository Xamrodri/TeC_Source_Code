%% Figure 2. equivalent of Jouberton et al. 2022

% needs:
% - Monthly or daily maps of: air temperature, sold and total
% precipitation, avalanches, snowmelt, icemelt, evaporation (Ice and snow)
% - Monthly or daily catchment discharge
% - Corresponding time vector (date_f2)

if ~exist('date_f2','var')
    date_f2 = date_m;
end 

if ~exist('hydro_year_f2_i','var')
    hydro_year_f2_i = 1;
end 

if ~exist('show_ratio','var')
    show_ratio = 1;
end 


ind_monsoon = ismember(month(date_f2), [6 7 8 9]); clear cmap
Q_mmh_to_m3s_catch = (nCatchPix.*(DEM.cellsize).^2/1000)/3600;
Years_no_f2 = unique(year(date_f2));

period_start = min(find(month(date_f2) == 10 & day(date_f2) == 1,1), find(month(date_f2) == 11 & day(date_f2) == 1,1));

period_end = find(month(date_f2) == 9,1,'last');

Hydro_years_f2 = unique(year(date_f2(period_start:period_end)));
Hydro_years_f2 = Hydro_years_f2(1:end-1);

if hydro_year_f2_i == 1; Years_no_f2 = Hydro_years_f2; end 

[Ta_yearly_y, Pr_y, Pr_sno_y, SMB_y, Imelt_y, Smelt_y, Q_y, ESN_y,ET_y] = deal(NaN(numel(Years_no_f2),1));

if ~isempty(Years_no_f2)

for yy = 1:length(Years_no_f2)

    if hydro_year_f2_i == 1
        ind_start = find(year(date_m) == Years_no_f2(yy) & month(date_m) == 10 & day(date_m) == 1,1);
        ind_end =   find(year(date_m) == Years_no_f2(yy)+1 & month(date_m) == 9 & day(date_m) == 1,1);
        if yy == 1; ind_start = find(year(date_m) == Years_no_f2(yy) & month(date_m) == 11 & day(date_m) == 1,1); end
    else 
        ind_start = find(year(date_m) == Years_no_f2(yy),1,'first');
        ind_end =   find(year(date_m) == Years_no_f2(yy),1,'last');
    end 

    Ta_yearly_y(yy) = nanmean(TA_map(:,:,ind_start:ind_end),'all');

    Pr_y_map = nansum(PRECIP_map(:,:,ind_start:ind_end),3); Pr_y_map(MASK ~=1) = NaN;
    Pr_sno_y_map = nansum(PSNOW_map(:,:,ind_start:ind_end),3); Pr_sno_y_map(MASK ~=1) = NaN;

    Pr_y(yy) = nanmean(Pr_y_map,'all'); 
    Pr_sno_y(yy) = nanmean(Pr_sno_y_map,'all'); 

    SMB_mapi = nansum(PSNOW_map(:,:,ind_start:ind_end),3) + nansum(AVA_map(:,:,ind_start:ind_end),3) ...
        - nansum(SMSm_map(:,:,ind_start:ind_end),3) - nansum(ESN_map(:,:,ind_start:ind_end),3) ...
        - nansum(SMG_map(:,:,ind_start:ind_end),3) - nansum(EICE_map(:,:,ind_start:ind_end),3);
    SMB_mapi(GLA_ID ~= gla_id) = NaN; % GMB at main glacier
    SMB_y(yy) = nanmean(SMB_mapi,'all');

    Imelt_fig_3 = nansum(SMG_map(:,:,ind_start:ind_end),3);
    Imelt_fig_3(GLA_ID ~= gla_id) = NaN; % GMB at main glacier
    Imelt_fig_3_y(yy) = nanmean(Imelt_fig_3,'all');   

    Imelt_y_map = nansum(SMG_map(:,:,ind_start:ind_end),3); Imelt_y_map(MASK ~=1) = NaN;
    Imelt_y(yy) = nanmean(Imelt_y_map,'all');

    Smelt_y_map = nansum(SMSm_map(:,:,ind_start:ind_end),3); Smelt_y_map(MASK ~=1) = NaN;
    Smelt_y(yy) = nanmean(Smelt_y_map,'all');

    ESN_y_map = nansum(ESN_map(:,:,ind_start:ind_end),3); ESN_y_map(MASK ~=1) = NaN;
    ESN_y(yy) = nanmean(ESN_y_map,'all');

    ET_y_map = nansum(ET_map(:,:,ind_start:ind_end),3); ET_y_map(MASK ~=1) = NaN;
    ET_y(yy) = nanmean(ET_y_map,'all');

    if exist('SPAVG_dm','var')
     ind_m_q = year(SPAVG_dm.Date) == Years_no_f2(yy);
     Q_y(yy) = nanmean(SPAVG_dm.Q_channel_tg(ind_m_q).*Q_mmh_to_m3s_catch);
    end 
end 

Pr_liq_y = Pr_y - Pr_sno_y;
Runoff_y = Pr_liq_y + Imelt_y + Smelt_y;

Imelt_ctb_y = Imelt_y./Runoff_y;
Smelt_ctb_y = Smelt_y./Runoff_y;
Rain_ctb_y = Pr_liq_y./Runoff_y;

%Prepare the colorbar for the monsoon temperature anomaly

Ta_devia = (Ta_yearly_y - nanmean(Ta_yearly_y))';

%Associate a color to each air temperature anonamly
cmab_rb = redblue(50);
t_anom = linspace(-1.5,1.5,50);
Ta_devia_col = NaN(numel(Ta_devia),3);

for t_i = 1:numel(Ta_devia)
   col_id = find(abs(Ta_devia(t_i) - t_anom) == min(abs(Ta_devia(t_i) - t_anom)),1);
   Ta_devia_col(t_i,:) = cmab_rb(col_id,:);
end 

r_snow_monsoon = Pr_sno_y./ Pr_y;
Mov_snowr_monsoon = movmean(r_snow_monsoon,10);

if exist('Q_y','var'); M_Q_tot = movmean(Q_y,10); end
M = movmean(SMB_y,10);

%% Figure 2

Hug_period_i = '2000-01-01_2020-01-01'; % 2000-2020
GMB_HUG_MAIN_mean = Hug_gmb.(['RGI_' num2str(gla_id)]).dmdtda(categorical(cellstr(Hug_period_i)) == Hug_gmb.(['RGI_' num2str(gla_id)]).period);
GMB_HUG_MAIN_err = Hug_gmb.(['RGI_' num2str(gla_id)]).err_dmdtda(categorical(cellstr(Hug_period_i)) == Hug_gmb.(['RGI_' num2str(gla_id)]).period);
GMB_HUG_MAIN_tc = nanmean(SMB_y(Years_no_f2>1998 & Years_no_f2 < 2019));
GMB_KH9_SRTM_tc = nanmean(SMB_y(Years_no_f2>1973 & Years_no_f2 < 2000));

FS = 11;
x_ticks = Years_no_f2(mod(Years_no_f2,5) == 0);


fi10 = figure('Renderer', 'painters', 'Position', [345.6667 46.3333 542 580.6667]);
tiledlayout(3,1,"TileSpacing",'compact')

nexttile
b = bar(Years_no_f2, Pr_y,'HandleVisibility','off','FaceColor','flat');
ylabel('Precipitation (mm)','FontSize',FS); ylim([0 1.2.*max(Pr_y)]); hold on; grid on;
yticks([0, 1000, 2000]); ylim([0 3000])
if show_ratio == 1
yyaxis right
plot(Years_no_f2, r_snow_monsoon,'--k','LineWidth',0.6','DisplayName','Snowfall ratio'); hold on;
plot(Years_no_f2, Mov_snowr_monsoon,'--k','LineWidth',1.1,'HandleVisibility','off'); b.FaceColor = 'flat';
h2 = plot(Years_no_f2, Mov_snowr_monsoon,'--k','LineWidth',1.1,'HandleVisibility','off');
h2 = plot(Years_no_f2, Mov_snowr_monsoon,'-k','LineWidth',1.1,'DisplayName','Snowfall ratio - 10yr average');
lg1 = legend('Box','off'); lg1.NumColumns = 2;
ylim([0 1]); yticks([0.1 0.4 0.7 1]); ylabel('Jun-Sep snowfall ratio (-)','FontSize',FS); 
end
b.CData = Ta_devia_col; b.EdgeColor = [0.5 0.5 0.5]; %clim([-1.5 1.5])
ax = gca; ax.YColor = [0 0 0]; 
xlim([Years_no_f2(1)-0.5 Years_no_f2(end)+0.5]); xticks(x_ticks)%xlabel('Year')
hc = colorbar('Location','northoutside');%
set(hc,'Position',[0.3 0.65 0.4 0.017]); 
annotation('textbox', [0.23 0.68 0 0], 'string', '-1.5°C','FontSize',8)
annotation('textbox', [0.70 0.68 0.3 0], 'string', '+1.5°C (Ta anomaly)','FontSize',8,'EdgeColor','none')
colormap(redblue(50))
set(gca, 'XAxisLocation', 'top')

nexttile
bar(Years_no_f2, SMB_y.*0.001,'FaceColor',[0.4 0.4 0.4],'DisplayName','Individual years'); hold on; grid on;
plot(Years_no_f2(1:min(5,length(Years_no_f2)/2)), M(1:min(5,length(Years_no_f2)/2)).*0.001,'--r','Marker','none','LineWidth',2,'HandleVisibility','off')
plot(Years_no_f2(end-min(4,length(Years_no_f2)/2):end), M(end-min(4,length(Years_no_f2)/2):end).*0.001,'--r','Marker','none','LineWidth',2,'HandleVisibility','off')
plot(Years_no_f2(min(5,length(Years_no_f2)/2):end-min(4,length(Years_no_f2)/2)), M(min(5,length(Years_no_f2)/2):end-min(4,length(Years_no_f2)/2)).*0.001,'-r','Marker','none','LineWidth',1.6,'DisplayName','10yr-average')
if strcmp(glacier,'Parlung4')
    shadedErrorBar(1974:2000, (1974:2000).*0 - 0.18, (1974:2000).*0+0.14,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(5,:)./255 0.45],'LineWidth',2,'DisplayName','KH-9 - SRTM','LineStyle','-'}); hold on
    plot(1974:2000, (1974:2000).*0 + GMB_KH9_SRTM_tc.*0.001,'Color',colorbrewer.qual.Dark2{1, 8}(5,:)./255,'LineWidth',2,'DisplayName','KH-9 - SRTM (modelled)','LineStyle','--'); hold on
elseif strcmp(glacier,'Kyzylsu')
    shadedErrorBar(1973:2000, (1973:2000).*0 +0.05, (1973:2000).*0+0.2,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(5,:)./255 0.45],'LineWidth',2,'DisplayName','KH-9 - SRTM','LineStyle','-'}); hold on
    plot(1973:2000, (1973:2000).*0 + GMB_KH9_SRTM_tc.*0.001,'Color',colorbrewer.qual.Dark2{1, 8}(5,:)./255,'LineWidth',2,'DisplayName','KH-9 - SRTM (modelled)','LineStyle','--'); hold on
end

shadedErrorBar(2000:2020, (2000:2020).*0+GMB_HUG_MAIN_mean, (2000:2020).*0+GMB_HUG_MAIN_err,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(3,:)./255 0.45],'LineWidth',2,'DisplayName','Hugonnet 2021','LineStyle','-'}); hold on
plot(2000:2020, (2000:2020).*0+GMB_HUG_MAIN_tc.*0.001,'Color',[colorbrewer.qual.Dark2{1, 8}(3,:)./255],'LineWidth',2,'DisplayName','Hugonnet 2021 (modelled)','LineStyle','--'); hold on

xlim([Years_no_f2(1)-0.5 Years_no_f2(end)+0.5]); xticks(x_ticks)%xlabel('Year')
set(gca,'XTickLabel',[])
ylabel('Glacier \Deltamass (m w.e.)','FontSize',FS); ylim([-2 1])
%set(c,'Position',[xo (yo+1*(space+height)+0.02) width height]); %set(gca,'XTickLabel',[])
lg1 = legend('Box','off');%,'Orientation','Horizontal')
lg1.NumColumns = 2; lg1.Location = 'SouthWest';

nexttile
b1 = bar(Years_no_f2, [Imelt_y Smelt_y Pr_liq_y],'stacked');
xlim([Years_no_f2(1)-0.5 Years_no_f2(end)+0.5]); xticks(x_ticks)%xlabel('Year')
b1(1, 2).FaceColor = [0.65 0.95 1];
b1(1, 1).FaceColor = [0.85 0.85 0.85];
b1(1, 3).FaceColor = [0 0.6 1]; grid on;
lg2 = legend('Icemelt','Snowmelt','Rain','Discharge','Discharge 10yr-avg','Orientation','horizontal');
lg2.ItemTokenSize = [16,10];
ylabel('Runoff [mm w.e./yr]','FontSize',FS)
exportgraphics(fi10,[dir_fig '\T&C_Figure2_Jouberton2022_' num2str(Years_no_f2(1)) '_' ...
   num2str(Years_no_f2(end)) '.png'],'Resolution',300,'BackgroundColor','none')

%% Breakpoint analysis

ind_change_Prsno = find(ischange(Pr_sno_y,'mean','MaxNumChanges',1));
ind_change_Ta = find(ischange(Ta_yearly_y,'mean','MaxNumChanges',1));
ind_change_Pr = find(ischange(Pr_y,'mean','MaxNumChanges',1));
ind_change_SMB = find(ischange(SMB_y,'mean','MaxNumChanges',1));

if ~isempty(ind_change_Prsno + ind_change_Ta + ind_change_Pr + ind_change_SMB)

fi10 = figure('Renderer', 'painters', 'Position', [345.6667 91 610.6666 536]);
tiledlayout(4,1,"TileSpacing","compact")
nexttile
plot(Years_no_f2,Pr_sno_y,'-sqk','MarkerFaceColor','k','MarkerSize',5); grid on; hold on;
plot([Years_no_f2(1) Years_no_f2(ind_change_Prsno)], [mean(Pr_sno_y(1:ind_change_Prsno)) mean(Pr_sno_y(1:ind_change_Prsno))],'--k')
xline(Years_no_f2(ind_change_Prsno),'r','LineWidth',1); ylim([0 1200]);
plot([Years_no_f2(ind_change_Prsno) Years_no_f2(end)], [mean(Pr_sno_y(ind_change_Prsno:end)) mean(Pr_sno_y(ind_change_Prsno:end))],'--k')
lg1 = legend('Hydro-year','Mean-bef','Breakpoint','Mean-Aft','Location','SouthWest');
lg1.NumColumns = 2;
xlim([1999 2023]); ylabel('Snowfall [mm]','FontSize',10); xticks(Years_no_f2(1):2:Years_no_f2(end))
nexttile
plot(Years_no_f2,Pr_y,'-sqk','MarkerFaceColor','k','MarkerSize',5); grid on; hold on;
plot([Years_no_f2(1) Years_no_f2(ind_change_Pr)], [mean(Pr_y(1:ind_change_Pr)) mean(Pr_y(1:ind_change_Pr))],'--k')
plot([Years_no_f2(ind_change_Pr) Years_no_f2(end)], [mean(Pr_y(ind_change_Pr:end)) mean(Pr_y(ind_change_Pr:end))],'--k')
xline(Years_no_f2(ind_change_Pr),'r','LineWidth',1); ylim([0 1800]);
xlim([1999 2023]); ylabel('Precipitation [mm]','FontSize',10); xticks(Years_no_f2(1):2:Years_no_f2(end))
nexttile
plot(Years_no_f2,Ta_yearly_y,'-sqk','MarkerFaceColor','k','MarkerSize',5); grid on; hold on;
%plot([Years_no_f2(1) Years_no_f2(ind_change_Ta)], [mean(Ta_yearly_y(1:ind_change_Ta)) mean(Ta_yearly_y(1:ind_change_Ta))],'--k')
%plot([Years_no_f2(ind_change_Ta) Years_no_f2(end)], [mean(Ta_yearly_y(ind_change_Ta:end)) mean(Ta_yearly_y(ind_change_Ta:end))],'--k')
plot([Years_no_f2(1) Years_no_f2(end-4)], [mean(Ta_yearly_y(1:end-4)) mean(Ta_yearly_y(1:end-4))],'--k')
plot([Years_no_f2(end-4) Years_no_f2(end)], [mean(Ta_yearly_y(end-4:end)) mean(Ta_yearly_y(end-4:end))],'--k')
xline(Years_no_f2(ind_change_Ta),'r','LineWidth',1); %ylim([0 1800]);
xlim([1999 2023]); ylabel('Ta [°C]','FontSize',10); xticks(Years_no_f2(1):2:Years_no_f2(end))
nexttile
plot(Years_no_f2,SMB_y,'-sqk','MarkerFaceColor','k','MarkerSize',5); grid on; hold on;
plot([Years_no_f2(1) Years_no_f2(ind_change_SMB)], [mean(SMB_y(1:ind_change_SMB)) mean(SMB_y(1:ind_change_SMB))],'--k')
plot([Years_no_f2(ind_change_SMB) Years_no_f2(end)], [mean(SMB_y(ind_change_SMB:end)) mean(SMB_y(ind_change_SMB:end))],'--k')
xline(Years_no_f2(ind_change_SMB),'r','LineWidth',1); %ylim([0 1800]);
xlim([1999 2023]); ylabel('SMB [m w.e.]','FontSize',10); xticks(Years_no_f2(1):2:Years_no_f2(end))
exportgraphics(fi10,[dir_fig '\Breakpoint_analysis_' num2str(Years_no_f2(1)) '_' ...
   num2str(Years_no_f2(end)) '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
end 
end 