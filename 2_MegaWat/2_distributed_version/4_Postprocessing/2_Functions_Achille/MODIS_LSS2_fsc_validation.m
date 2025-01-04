%% MODIS and L8/S2 fractional snow cover figure

%Date_d : timestamp of snow cover maps
% date_tg : timestamp of precipitation, temperature

if exist('SPAVG_dm','var') % Fully distributed simulations
    date_tg = SPAVG_dm.Date; % Corresponds to the retime to daily of catchment output table
    date_d = SPAVG_dm.Date; % Corresponds to the retime to daily of catchment output table
    Ta_tg = SPAVG_dm.Ta_tg; % Mean daily air temperature
    Pr_tg = SPAVG_ds.Pr_tg; % Sum of daily precipitation
    Pr_sno_tg = SPAVG_ds.Pr_sno_tg; % Sum of daily solid precipitation
    Date_d = date_h-hours(12); % Corresponds to daily map timestamp
elseif exist('spatialize','var') && spatialize == 1 % cluster simulation
    date_tg = DATES.TA_map; % Corresponds to the retime to daily of catchment output table
    date_d = DATES.TA_map; % Corresponds to the retime to daily of catchment output table
    Ta_tg = squeeze(nanmean(TA_map,[1 2])); % Mean daily air temperature
    Pr_tg = squeeze(nanmean(PRECIP_map,[1 2])); % Sum of daily precipitation
    Pr_sno_tg = squeeze(nanmean(PSNOW_map,[1 2])); % Sum of daily solid precipitation
end 

if exist('date_h')
    Date_d = date_h-hours(12); % Corresponds to daily map timestamp
end 


datenum_modis = datenum(modis_fsc_40.dateTime);
datenum_modis(isnan(Modis_fsc_mv)) = NaN;

min_precip = 0;
max_precip = round(1.5*max(Pr_tg),-1);

%%% Compute performance metrics

% First make a synchronize table of fsc and SLA

LS_fsc_comp = table2timetable(LS_sno(:,[1 3])); LS_fsc_comp.Properties.VariableNames = {'fsc_l8s2'};
MODIS_fsc_comp = timetable(modis_fsc_40.dateTime, fsc_modis_filt_40,'VariableNames',{'fsc_MODIS'});
TC_fsc_comp = timetable(Date_d,scas,'VariableNames',{'fsc_TC'});
fsc_comp = synchronize(TC_fsc_comp, MODIS_fsc_comp, LS_fsc_comp); clear LS_fsc_comp MODIS_fsc_comp TC_fsc_comp

LS_sla_comp = table2timetable(LS_sno(:,[1 2])); LS_sla_comp.Properties.VariableNames = {'sla_l8s2'};
LS_sla_comp.sla_l8s2(LS_sno.fsnow_hybrid<0.05) = NaN;
MODIS_sla_comp = timetable(modis_fsc_40.dateTime, se_modis_filt_40,'VariableNames',{'sla_MODIS'});
TC_sla_comp = timetable(Date_d,se,'VariableNames',{'sla_TC'});
sla_comp = synchronize(TC_sla_comp, MODIS_sla_comp, LS_sla_comp); clear LS_sla_comp MODIS_sla_comp TC_sla_comp

% Compute performance metrics

ME_fsc_MODIS = nanmean(fsc_comp.fsc_TC - fsc_comp.fsc_MODIS);
ME_fsc_l8s2 = nanmean(fsc_comp.fsc_TC - fsc_comp.fsc_l8s2);

RMSE_fsc_MODIS = rmse(fsc_comp.fsc_TC, fsc_comp.fsc_MODIS,"omitnan");
RMSE_fsc_l8s2 = rmse(fsc_comp.fsc_TC, fsc_comp.fsc_l8s2,"omitnan");

r_fsc_MODIS = fitlm(fsc_comp.fsc_TC,fsc_comp.fsc_MODIS).Rsquared.Adjusted;
r_fsc_l8s2 = fitlm(fsc_comp.fsc_TC,fsc_comp.fsc_l8s2).Rsquared.Adjusted;

disp(['Snow cover performance: Mean-Error (MODIS, L8S2) = ' num2str(ME_fsc_MODIS,2) ', ' num2str(ME_fsc_l8s2,2) ' [-]'])
disp(['Snow cover performance: R^2 (MODIS, L8S2) = ' num2str(r_fsc_MODIS,2) ', ' num2str(r_fsc_l8s2,2) ' [-]'])
disp(['Snow cover performance: RMSE (MODIS, L8S2) = ' num2str(RMSE_fsc_MODIS,2) ', ' num2str(RMSE_fsc_l8s2,2) ' [-]'])

% Snow line elevation metrics
ME_sla_MODIS = nanmean(sla_comp.sla_TC - sla_comp.sla_MODIS);
ME_sla_l8s2 = nanmean(sla_comp.sla_TC - sla_comp.sla_l8s2);

RMSE_sla_MODIS = rmse(sla_comp.sla_TC, sla_comp.sla_MODIS,"omitnan");
RMSE_sla_l8s2 = rmse(sla_comp.sla_TC, sla_comp.sla_l8s2,"omitnan");

r_sla_MODIS = fitlm(sla_comp.sla_TC,sla_comp.sla_MODIS).Rsquared.Adjusted;
r_sla_l8s2 = fitlm(sla_comp.sla_TC,sla_comp.sla_l8s2).Rsquared.Adjusted;

disp(['Snow line altitude performance: Mean-Error (MODIS, L8S2) = ' num2str(ME_sla_MODIS,2) ', ' num2str(ME_sla_l8s2,2) ' [m]'])
disp(['Snow line altitude performance: R^2 (MODIS, L8S2) = ' num2str(r_sla_MODIS,2) ', ' num2str(r_sla_l8s2,2) ' [m]'])
disp(['Snow line altitude: RMSE (MODIS, L8S2) = ' num2str(RMSE_sla_MODIS,2) ', ' num2str(RMSE_sla_l8s2,2) ' [m]'])

%%% If simulation period is shorter than two years %%%%%%%%%%%%
if nYears <2.2

fi2 = figure('Renderer', 'painters', 'Position', [89 419.6667 1108 277.3333]) ;
colororder([0 0 0; 0 0.6 1])
yyaxis left;
shadedErrorBar(datenum_modis,Modis_fsc_mv.*100, [(Modis_fsc_mv_25-Modis_fsc_mv).*100 (-Modis_fsc_mv_45+Modis_fsc_mv).*100] ,'LineProps',{'k','LineWidth',0.6}); hold on; grid on;
plot(datenum(Date_d),movmean(scas.*100,10),'-r'); hold on;
scatter(datenum(LS_sno.dateTime), LS_sno.fsnow_hybrid.*100,30,'+k','LineWidth',1.3)
xticks(linspace(datenum(Date_d(1)),datenum(Date_d(end)),10));  % date_end NHM = 14-Sep-2019 12:00:00
datetick('x','mmm-yyyy','keepticks')
ylabel('Snow covered area [%]','FontSize',10)
xlim([datenum(Date_d(1)) datenum(Date_d(end))]);
grid on; ylim([0 105])

s1 = surface([datenum(date_tg)';datenum(date_tg)'],[102*ones(size(date_tg))'; 105*ones(size(date_tg))'],...
    [zeros(size(date_tg))';zeros(size(date_tg))'], [Ta_tg';Ta_tg'],...
        'facecol','no',...
        'edgecol','interp',...
        'linew',2);
cm2 = colormap(redblue); clim([-20 10]); cb2 = colorbar('eastoutside');
ylabel(cb2,'Daily Ta [°C]','FontSize',10)

% 
yyaxis right
b1 = bar(datenum(date_d),Pr_sno_tg,1,'FaceColor',[0 0.6 1],'EdgeColor','none','LineStyle','-');
b2 = bar(datenum(date_d),Pr_tg,1,'FaceColor','none','EdgeColor',[0.7 0.7 0.7],'LineStyle','-');
ylim([min_precip max_precip])
lg1 = legend('MODIS','T&C','L8/S2','Location','northwest'); lg1.NumColumns = 2;
ylabel('Snowfall [mm]','FontSize',10)
title({[glacier '-catchment ' num2str(year(Date_d(1))) '-' num2str(year(Date_d(end))) ' - ' forcing_nm], ' '})
exportgraphics(fi2,[dir_fig '\Snow_cover\MODIS_L8S2_vs_TC_fsc.png'],'Resolution',300,'BackgroundColor','none')

%%% If simulation period is shorter longer than 2 years but shorter than 20 years %%%%%%%%%%%%

elseif (length(Years_no) < 20) && (nYears >= 2.2)

fi2 = figure('Renderer', 'painters', 'Position', [11 50.3333 1.2053e+03 656.0000]) ;
tiledlayout(ceil(length(Years_no)./2),2,'TileSpacing','compact');

for yy = 1:length(Years_no)
nexttile
id_start = find(year(Date)==Years_no(yy),1);
id_end = find(year(Date)==Years_no(yy),1,'last');
colororder([0 0 0; 0 0.6 1])
scatter(datenum_modis,fsc_modis_filt_40*100,7,'sq','MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
plot([datenum_modis'; datenum_modis'], [(fsc_modis_filt_40'-(fsc_modis_filt_40'-fsc_modis_filt_25')); (fsc_modis_filt_40'-(fsc_modis_filt_40'-fsc_modis_filt_45'))].*100,'k','HandleVisibility','off','Color',[0.3 0.3 0.3],'LineStyle','-');%,'sr','MarkerFaceColor','r','MarkerSize',2); hold on; grid on;
plot(datenum(Date_d),scas.*100,'-r','LineWidth',1.1);
scatter(datenum(LS_sno.dateTime), LS_sno.fsnow_hybrid.*100,30,'+k','LineWidth',1.5,'MarkerEdgeColor',[62 150 81]./255)
xticks(linspace(datenum(Date(id_start)),datenum(Date(id_end)),7));
datetick('x','mmm','keepticks')
if mod(yy,2)==1; ylabel('Snow cover [%]'); end 
xlim([datenum(datetime(Years_no(yy),1,1)) datenum(datetime(Years_no(yy),12,31))]);
grid on; ylim([0 105])
text(datenum(datetime(Years_no(yy),1,15)), 35, num2str(Years_no(yy)),'FontSize',10,'FontWeight','bold')

end
lg1 = legend('MODIS','T&C','Landsat/S2'); 
lg1.NumColumns = 3; 
lg1.Position = [0.1668    0.94    0.3017    0.0254];

exportgraphics(fi2,[dir_fig '\Snow_cover\MODIS_L8S2_vs_TC_SLA_fsc.png'],'Resolution',300,'BackgroundColor','none')

else
%%% If simulation period is than 20 years %%%%%%%%%%%%

Years_kept = Years_no(Years_no > 1999); 

fi2 = figure('Renderer', 'painters', 'Position', [11 50.3333 761.3333 656.0000]) ;
tiledlayout(ceil(length(Years_kept)./4),2,'TileSpacing','compact');

for yy = 1:ceil(length(Years_kept)/2)
nexttile
id_start = find(year(Date)==Years_kept(yy),1);
id_end = find(year(Date)==Years_kept(yy),1,'last');
scatter(datenum_modis,fsc_modis_filt_40*100,7,'sq','MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
plot([datenum_modis'; datenum_modis'], [(fsc_modis_filt_40'-(fsc_modis_filt_40'-fsc_modis_filt_25')); (fsc_modis_filt_40'-(fsc_modis_filt_40'-fsc_modis_filt_45'))].*100,'k','HandleVisibility','off','Color',[0.3 0.3 0.3],'LineStyle','-');%,'sr','MarkerFaceColor','r','MarkerSize',2); hold on; grid on;
plot(datenum(Date_d),scas.*100,'-r','LineWidth',1.1);
scatter(datenum(LS_sno.dateTime), LS_sno.fsnow_hybrid.*100,30,'+k','LineWidth',1.5,'MarkerEdgeColor',[62 150 81]./255)
xticks(linspace(datenum(Date(id_start)),datenum(Date(id_end)),7));
datetick('x','mmm','keepticks')
if mod(yy,2)==1; ylabel('Snow cover [%]'); end 
xlim([datenum(datetime(Years_kept(yy),1,1)) datenum(datetime(Years_kept(yy),12,31))]);
grid on; ylim([0 105])
text(datenum(datetime(Years_kept(yy),1,15)), 35, num2str(Years_kept(yy)),'FontSize',10,'FontWeight','bold')

end
lg1 = legend('MODIS','T&C','Landsat/S2'); 
lg1.NumColumns = 3; 
lg1.Position = [0.1668    0.94    0.3017    0.0254];

exportgraphics(fi2,[dir_fig '\Snow_cover\MODIS_L8S2_vs_TC_SLA_SI_fsc_' num2str(Years_kept(1)) '-' num2str(Years_kept(ceil(length(Years_kept)/2))) ...
    '.png'],'Resolution',300,'BackgroundColor','none')


fi2 = figure('Renderer', 'painters', 'Position', [11 50.3333 761.3333 656.0000]) ;
tiledlayout(ceil(length(Years_kept)./4),2,'TileSpacing','compact');

for yy = (1+ceil(length(Years_kept)/2)):length(Years_kept)
nexttile
id_start = find(year(Date)==Years_kept(yy),1);
id_end = find(year(Date)==Years_kept(yy),1,'last');
scatter(datenum_modis,fsc_modis_filt_40*100,7,'sq','MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
plot([datenum_modis'; datenum_modis'], [(fsc_modis_filt_40'-(fsc_modis_filt_40'-fsc_modis_filt_25')); (fsc_modis_filt_40'-(fsc_modis_filt_40'-fsc_modis_filt_45'))].*100,'k','HandleVisibility','off','Color',[0.3 0.3 0.3],'LineStyle','-');%,'sr','MarkerFaceColor','r','MarkerSize',2); hold on; grid on;
plot(datenum(Date_d),scas.*100,'-r','LineWidth',1.1);
scatter(datenum(LS_sno.dateTime), LS_sno.fsnow_hybrid.*100,30,'+k','LineWidth',1.5,'MarkerEdgeColor',[62 150 81]./255)
xticks(linspace(datenum(Date(id_start)),datenum(Date(id_end)),7));
datetick('x','mmm','keepticks')
if mod(yy,2)==1; ylabel('Snow cover [%]'); end 
xlim([datenum(datetime(Years_kept(yy),1,1)) datenum(datetime(Years_kept(yy),12,31))]);
grid on; ylim([0 105])
text(datenum(datetime(Years_kept(yy),1,15)), 35, num2str(Years_kept(yy)),'FontSize',10,'FontWeight','bold')
end
lg1 = legend('MODIS','T&C','Landsat/S2'); 
lg1.NumColumns = 3; 
lg1.Position = [0.1668    0.94    0.3017    0.0254];

exportgraphics(fi2,[dir_fig '\Snow_cover\MODIS_L8S2_vs_TC_SI_fsc_' num2str(Years_kept((1+ceil(length(Years_kept)/2)))) '-' num2str(Years_kept(end)) ...
    '.png'],'Resolution',300,'BackgroundColor','none')
end 

%% MODIS and L8/S2 snowline elevation figure

min_SLA = min([min(se), min(se_modis_filt_25), min(LS_sno.se_hybrid)]) ;
max_SLA = max([max(se), max(se_modis_filt_45), max(LS_sno.se_hybrid)]) ;

if nYears <2 

fi2 = figure('Renderer', 'painters', 'Position', [89 419.6667 1108 277.3333]) ;
nexttile
shadedErrorBar(datenum_modis,Modis_se_mv, [(Modis_se_mv_25-Modis_se_mv) (-Modis_se_mv_45+Modis_se_mv)] ,'LineProps',{'k','LineWidth',0.6}); hold on; grid on;
plot(datenum(Date_d),se,'r');
scatter(datenum(LS_sno.dateTime(LS_sno.fsnow_hybrid>0.05)), LS_sno.se_hybrid(LS_sno.fsnow_hybrid>0.05),30,'+k','LineWidth',1.3)
xticks(linspace(datenum(Date_d(1)),datenum(Date_d(end)),10));
datetick('x','mmm-yyyy','keepticks'); grid on; 
ylabel('Snow line elevation [m a.s.l.]')
xlim([datenum(Date_d(1)) datenum(Date_d(end))]);
title([glacier ' catchment ' num2str(year(Date_d(1))) '-' num2str(year(Date_d(end))) ' - ' forcing_nm])
exportgraphics(fi2,[dir_fig '\MODIS_L8S2_vs_TC_SLA.png'],'Resolution',300,'BackgroundColor','none')

elseif length(Years_no) < 20

fi2 = figure('Renderer', 'painters', 'Position', [11 50.3333 980 656]) ;
tiledlayout(ceil(length(Years_no)./2),2,'TileSpacing','compact');

for yy = 1:length(Years_no)
nexttile
id_start = find(year(Date)==Years_no(yy),1);
id_end = find(year(Date)==Years_no(yy),1,'last');
% shadedErrorBar(datenum_modis,Modis_se_mv, [(Modis_se_mv_25-Modis_se_mv) (-Modis_se_mv_45+Modis_se_mv)] ,'LineProps',{'k','LineWidth',0.6}); hold on; grid on;
scatter(datenum(modis_fsc_40.dateTime),se_modis_filt_40,10,'sq','MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
plot([datenum(modis_fsc_40.dateTime)'; datenum(modis_fsc_40.dateTime)'], [(se_modis_filt_40'-(se_modis_filt_40'-se_modis_filt_25')); (se_modis_filt_40'-(se_modis_filt_40'-se_modis_filt_45'))],'k','HandleVisibility','off','Color',[0.3 0.3 0.3]);%,'sr','MarkerFaceColor','r','MarkerSize',2); hold on; grid on;
plot(datenum(Date_d),se,'-r','LineWidth',1.1);
scatter(datenum(LS_sno.dateTime(LS_sno.fsnow_hybrid>0.05)), LS_sno.se_hybrid(LS_sno.fsnow_hybrid>0.05),30,'+','MarkerEdgeColor',[62 150 81]./255,'LineWidth',1.3)
xticks(linspace(datenum(Date(id_start)),datenum(Date(id_end)),7));
datetick('x','mmm','keepticks'); ylim([min_SLA max_SLA])
if ismember(yy, 1:2:length(Years_no)); ylabel('Elevation [m a.s.l.]'); end
xlim([datenum(datetime(Years_no(yy),1,1)) datenum(datetime(Years_no(yy),12,31))]);
grid on; 
text(datenum(datetime(Years_no(yy),1,23)), 5700, num2str(Years_no(yy)),'FontSize',10,'FontWeight','bold')
if yy == 1
lg1 = legend('MODIS','T&C','L8/S2'); 
lg1.NumColumns = 3; 
end 
if ismember(yy, [1 2])
    title({[glacier ' catchment ' num2str(year(Date_d(1))) '-' num2str(year(Date_d(end))) ' - ' forcing_nm], ' '})
end 

end
exportgraphics(fi2,[dir_fig '\Snow_cover\MODIS_L8S2_vs_TC_SLA.png'],'Resolution',300,'BackgroundColor','none')

else

Years_kept = Years_no(Years_no > 1999); 

fi2 = figure('Renderer', 'painters', 'Position', [11 50.3333 761.3333 656.0000]) ;
tiledlayout(ceil(length(Years_kept)./4),2,'TileSpacing','compact');

for yy = (1+ceil(length(Years_kept)/2)):length(Years_kept)
nexttile
id_start = find(year(Date)==Years_kept(yy),1);
id_end = find(year(Date)==Years_kept(yy),1,'last');
% shadedErrorBar(datenum_modis,Modis_se_mv, [(Modis_se_mv_25-Modis_se_mv) (-Modis_se_mv_45+Modis_se_mv)] ,'LineProps',{'k','LineWidth',0.6}); hold on; grid on;
scatter(datenum(modis_fsc_40.dateTime),se_modis_filt_40,10,'sq','MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
plot([datenum(modis_fsc_40.dateTime)'; datenum(modis_fsc_40.dateTime)'], [(se_modis_filt_40'-(se_modis_filt_40'-se_modis_filt_25')); (se_modis_filt_40'-(se_modis_filt_40'-se_modis_filt_45'))],'k','HandleVisibility','off','Color',[0.3 0.3 0.3]);%,'sr','MarkerFaceColor','r','MarkerSize',2); hold on; grid on;
plot(datenum(Date_d),se,'-r','LineWidth',1.1);
scatter(datenum(LS_sno.dateTime(LS_sno.fsnow_hybrid>0.05)), LS_sno.se_hybrid(LS_sno.fsnow_hybrid>0.05),30,'+','MarkerEdgeColor',[62 150 81]./255,'LineWidth',1.3)
xticks(linspace(datenum(Date(id_start)),datenum(Date(id_end)),7));
datetick('x','mmm','keepticks'); ylim([min_SLA max_SLA])
if mod(yy,4)==1; ylabel('Elevation [m a.s.l.]'); end
xlim([datenum(datetime(Years_kept(yy),1,1)) datenum(datetime(Years_kept(yy),12,31))]);
grid on; yticks([3000 4000 5000])
text(datenum(datetime(Years_kept(yy),1,23)), 5400, num2str(Years_kept(yy)),'FontSize',10,'FontWeight','bold')
end 
lg1 = legend('MODIS','T&C','Landsat/S2'); 
lg1.NumColumns = 3; 
lg1.Position = [0.1668    0.94    0.3017    0.0254];

exportgraphics(fi2,[dir_fig '\Snow_cover\MODIS_L8S2_vs_TC_SLA_' num2str(Years_kept((1+ceil(length(Years_kept)/2)))) '-' num2str(Years_kept(end)) ...
    '.png'],'Resolution',300,'BackgroundColor','none')

fi3 = figure('Renderer', 'painters', 'Position', [11 50.3333 761.3333 656.0000]) ;
tiledlayout(ceil(length(Years_kept)./4),2,'TileSpacing','compact');

for yy = 1:ceil(length(Years_kept)/2)
nexttile
id_start = find(year(Date)==Years_kept(yy),1);
id_end = find(year(Date)==Years_kept(yy),1,'last');
% shadedErrorBar(datenum_modis,Modis_se_mv, [(Modis_se_mv_25-Modis_se_mv) (-Modis_se_mv_45+Modis_se_mv)] ,'LineProps',{'k','LineWidth',0.6}); hold on; grid on;
scatter(datenum(modis_fsc_40.dateTime),se_modis_filt_40,10,'sq','MarkerFaceColor','k','MarkerEdgeColor','k'); hold on
plot([datenum(modis_fsc_40.dateTime)'; datenum(modis_fsc_40.dateTime)'], [(se_modis_filt_40'-(se_modis_filt_40'-se_modis_filt_25')); (se_modis_filt_40'-(se_modis_filt_40'-se_modis_filt_45'))],'k','HandleVisibility','off','Color',[0.3 0.3 0.3]);%,'sr','MarkerFaceColor','r','MarkerSize',2); hold on; grid on;
plot(datenum(Date_d),se,'-r','LineWidth',1.1);
scatter(datenum(LS_sno.dateTime(LS_sno.fsnow_hybrid>0.05)), LS_sno.se_hybrid(LS_sno.fsnow_hybrid>0.05),30,'+','MarkerEdgeColor',[62 150 81]./255,'LineWidth',1.3)
xticks(linspace(datenum(Date(id_start)),datenum(Date(id_end)),7));
datetick('x','mmm','keepticks'); ylim([min_SLA max_SLA])
if mod(yy,4)==1; ylabel('Elevation [m a.s.l.]'); end
xlim([datenum(datetime(Years_kept(yy),1,1)) datenum(datetime(Years_kept(yy),12,31))]);
grid on; yticks([3000 4000 5000])
text(datenum(datetime(Years_kept(yy),1,23)), 5400, num2str(Years_kept(yy)),'FontSize',10,'FontWeight','bold')
end 
lg1 = legend('MODIS','T&C','Landsat/S2'); 
lg1.NumColumns = 3; 
lg1.Position = [0.1668    0.94    0.3017    0.0254];

exportgraphics(fi3,[dir_fig '\Snow_cover\MODIS_L8S2_vs_TC_SLA_' num2str(Years_kept(1)) '-' num2str(Years_kept(ceil(length(Years_kept)/2))) ...
    '.png'],'Resolution',300,'BackgroundColor','none')

end

