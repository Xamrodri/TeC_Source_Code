%  Needs:
% - the name of the outlet pixel (outlet_nm)
% - monthly or daily maps of liquid precip, snowmelt and icemelt

Q_mmh_to_m3s_catch = (nCatchPix.*(DEM.cellsize).^2/1000)/3600;
Q_mmh_to_m3s = ((DEM.cellsize).^2/1000)/3600;

 % Load discharge output
if labelled_output == 0
    outlet_table = readtable([dir_tcout '\POSTPROCESSED\TABLES\TrackedPixel___Streamgauge_lake___NaNm']);
else   
    outlet_table = readtable([dir_tcout '\OUTPUT_' glacier '_PIXEL_' outlet_nm]);
end 

outlet_h = table2timetable(outlet_table); %outlet_h.NA = [];
outlet_m = retime(outlet_h,'monthly',@nanmean);
outlet_y = retime(outlet_h,'yearly',@nanmean);

Q_mm_to_m3s = 0.001*(sim_res*sim_res)/(3600);
Q_mm_to_m3s_catch = 0.001*nCatchPix*(sim_res*sim_res)/(3600);

if strcmp(glacier,'Kyzylsu')
        outlet_table_gla = readtable([dir_tcout '\OUTPUT_' glacier '_PIXEL_Gauge1']);
        outlet_gla_h = table2timetable(outlet_table_gla); %outlet_h.NA = [];
        outlet_gla_m = retime(outlet_gla_h,'monthly',@nanmean);
        outlet_gla_y = retime(outlet_gla_h,'yearly',@nanmean);
end 
 
 % Compute monthly average of icemelt, snowmelt and precipitation
 
Sms_tt = timetable(date_m, SMSm_mean,'VariableNames',{'SMS'});
Sms_tt = retime(Sms_tt,'monthly',@nansum);
Smg_tt = timetable(date_m, SMG_mean,'VariableNames',{'SMG'}); Smg_tt = retime(Smg_tt,'monthly',@nansum);
Precip_tt = timetable(date_m,Pliq_mean, Psnow_mean,'VariableNames',{'Pliq','Psno'});
Precip_tt = retime(Precip_tt,'monthly',@nansum);

Hydro_m = synchronize(Sms_tt, Smg_tt, Precip_tt);
Hydro_y = retime(Hydro_m,'yearly',@nansum);
%  Hydrograph figure

if length(Years_no) < 2
    
fi4 = figure('Renderer', 'painters', 'Position', [86.3333 145 1126 265.3333]);
a1 = area(Hydro_m.date_m, [Hydro_m.Pliq  Hydro_m.SMS Hydro_m.SMG]); hold on; 
a1(3).FaceColor = [0.85 0.85 0.85]; a1(2).FaceColor = [0.65 0.95 1]; a1(1).FaceColor = [0 0.6 1];
a1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none';
ylabel('Runoff generation (mm/month)','FontSize',11); ylim([0 500])
yyaxis right
plot(outlet_m.Date(2:end) + days(30), outlet_m.QpointC(2:end)*Q_mm_to_m3s,'Color',[0 0 0],'LineWidth',1.3); hold on; grid on;
ax = gca; ax.YColor = [0 0 0]; ylim([0 round(1.2*max(outlet_m.QpointC(2:end)*Q_mm_to_m3s),-1)]);
ylabel('Discharge [m^3/s]','FontSize',11)
lg2 = legend('Rain','Snowmelt','Icemelt','Outlet discharge');
lg2.NumColumns = 2; lg2.Location = 'NorthEast';
title([glacier ' catchment ' num2str(year(date_h(1))) '-' num2str(year(date_h(end)))])
exportgraphics(fi4,[dir_fig '\T&C_hydrograph.png'],'Resolution',300,'BackgroundColor','none')

else 

fi4 = figure('Renderer', 'painters', 'Position', [84.3333 107 969.3333 473.3333]);
tiledlayout(2,1,'TileSpacing','compact')
nexttile
a1 = area(Hydro_m.date_m, [Hydro_m.Pliq  Hydro_m.SMS Hydro_m.SMG]); hold on; 
a1(3).FaceColor = [0.85 0.85 0.85]; a1(2).FaceColor = [0.65 0.95 1]; a1(1).FaceColor = [0 0.6 1];
a1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none';
ylabel('Runoff generation (mm/month)','FontSize',11); ylim([0 500])
yyaxis right
plot(outlet_m.Date(2:end) + days(30), outlet_m.QpointC(2:end)*Q_mm_to_m3s,'Color',[0 0 0],'LineWidth',1.3); hold on; grid on;
% plot(outlet_gla_m.Date(2:end) + days(30), outlet_gla_m.QpointC(2:end)*Q_mm_to_m3s,'Color',[1 0 0],'LineWidth',1.3); hold on; grid on;
ax = gca; ax.YColor = [0 0 0]; ylim([0 round(1.2*max(outlet_m.QpointC(2:end)*Q_mm_to_m3s),-1)]);
ylabel('Discharge [m^3/s]','FontSize',11)
lg2 = legend('Rain','Snowmelt','Icemelt','Outlet discharge');
lg2.NumColumns = 2; lg2.Location = 'NorthEast';
title([glacier ' catchment ' num2str(year(date_h(1))) '-' num2str(year(date_h(end)))])
xlim([datetime(Years_no(1),9,30) datetime(Years_no(ceil(length(Years_no)/2)),9,30)])

nexttile
a1 = area(Hydro_m.date_m, [Hydro_m.Pliq  Hydro_m.SMS Hydro_m.SMG]); hold on; 
a1(3).FaceColor = [0.85 0.85 0.85]; a1(2).FaceColor = [0.65 0.95 1]; a1(1).FaceColor = [0 0.6 1];
a1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none';
ylabel('Runoff generation (mm/month)','FontSize',11); ylim([0 500])
yyaxis right
plot(outlet_m.Date(2:end) + days(30), outlet_m.QpointC(2:end)*Q_mm_to_m3s,'Color',[0 0 0],'LineWidth',1.3); hold on; grid on;
% plot(outlet_gla_m.Date(2:end) + days(30), outlet_gla_m.QpointC(2:end)*Q_mm_to_m3s,'Color',[1 0 0],'LineWidth',1.3); hold on; grid on;
ax = gca; ax.YColor = [0 0 0]; ylim([0 round(1.2*max(outlet_m.QpointC(2:end)*Q_mm_to_m3s),-1)]);
ylabel('Discharge [m^3/s]','FontSize',11)
xlim([datetime(Years_no(ceil(length(Years_no)/2)),9,30) datetime(Years_no(end),9,30)])
exportgraphics(fi4,[dir_fig '\Hyddrograph_outlet_runoff.png'],'Resolution',300,'BackgroundColor','none')    

end
%% Investigate relationship between catchment outlet discharge and runoff component magnitude

r_sms_Q = fitlm(Hydro_y.SMS(2:end),outlet_y.QpointC(2:end)*Q_mm_to_m3s).Rsquared.Adjusted;
r_smg_Q = fitlm(Hydro_y.SMG(2:end),outlet_y.QpointC(2:end)*Q_mm_to_m3s).Rsquared.Adjusted;
r_pliq_Q = fitlm(Hydro_y.Pliq(2:end),outlet_y.QpointC(2:end)*Q_mm_to_m3s).Rsquared.Adjusted;
r_psno_Q = fitlm(Hydro_y.Psno(2:end),outlet_y.QpointC(2:end)*Q_mm_to_m3s).Rsquared.Adjusted;
r_sm_Q = fitlm(Hydro_y.SMS(2:end)+Hydro_y.SMG(2:end),outlet_y.QpointC(2:end)*Q_mm_to_m3s).Rsquared.Adjusted;

%% Compute mean monthly runoff for the two periods, at Kyzylsu outlet and catchment outlet

for yy = 1:length(Years_no_seas)-1
    for mm = 1:12

      if hydro_month(mm) > 9
        ind_start = find(year(date_m) == Years_no_seas(yy) & month(date_m) == hydro_month(mm),1,'first');
        ind_end =   find(year(date_m) == Years_no_seas(yy) & month(date_m) == hydro_month(mm),1,'last');

        date_ym(yy,mm) = datetime(Years_no_seas(yy), hydro_month(mm),15);
      else
        ind_start = find(year(date_m) == (Years_no_seas(yy)+1) & month(date_m) == hydro_month(mm),1,'first');
        ind_end =   find(year(date_m) == (Years_no_seas(yy)+1) & month(date_m) == hydro_month(mm),1,'last');

        date_ym(yy,mm) = datetime(Years_no_seas(yy)+1, hydro_month(mm),15);
      end 

    Q_yearly_ym(yy,mm) = nanmean(outlet_m.QpointC(ind_start:ind_end)*Q_mm_to_m3s);
    Q_yearly_gla_ym(yy,mm) = nanmean(outlet_gla_m.QpointC(ind_start:ind_end)*Q_mm_to_m3s);
    end
end 

if ~exist('date_start_runoff_p1','var') || ~exist('date_end_runoff_p1','var')
    date_start_runoff_p1 = date_m(1);
    date_end_runoff_p1 = date_m(226);
end 

if ~exist('date_start_runoff_p2','var') || ~exist('date_end_runoff_p2','var')
    date_start_runoff_p2 = date_m(227);
    date_end_runoff_p2 = date_m(end);
end 

ind_p1 = (date_ym > date_start_runoff_p1) & (date_ym < date_end_runoff_p1);
ind_p2 = (date_ym > date_start_runoff_p2) & (date_ym < date_end_runoff_p2);

Q_ym_p1 = Q_yearly_ym; Q_ym_p1(~ind_p1) = NaN;
Q_ym_p2 = Q_yearly_ym; Q_ym_p2(~ind_p2) = NaN;

Q_gla_ym_p1 = Q_yearly_gla_ym; Q_gla_ym_p1(~ind_p1) = NaN;
Q_gla_ym_p2 = Q_yearly_gla_ym; Q_gla_ym_p2(~ind_p2) = NaN;

%% Discharge difference between the two periods

color_p2 = (colorbrewer.qual.Set1{1, 9}(1,:))./255;
color_p1 = (colorbrewer.qual.Set1{1, 9}(2,:))./255;

fi2 = figure('Renderer', 'painters', 'Position',[131.6667 304.3333 798 333.3333]);
tiledlayout(1,2,'TileSpacing','compact')
nexttile
plot(1:12, nanmean(Q_ym_p1,1),'Color',color_p1,'LineWidth',1.2); grid on; hold on;
plot(1:12, nanmean(Q_ym_p2,1),'Color',color_p2,'LineWidth',1.2); ylim([0 13])
set(gca,'XTickLabels',hydro_month_labels(1:2:12));
ylabel('Discharge [m^3/s]','FontSize',11)
title('Kyzylsu - Cathment outlet')
xticks([1:2:12]); xlim([1 12]);
nexttile
plot(1:12, nanmean(Q_gla_ym_p1,1),'Color',color_p1,'LineWidth',1.2); grid on; hold on;
plot(1:12, nanmean(Q_gla_ym_p2,1),'Color',color_p2,'LineWidth',1.2); ylim([0 13]) 
set(gca,'XTickLabels',hydro_month_labels(1:2:12));
ylabel('Discharge [m^3/s]','FontSize',11)
title('Kyzylsu - Proglacial stream')
xticks([1:2:12]); xlim([1 12]);
exportgraphics(fi2,[dir_fig '\Discharge_2periods_outlet_progla.png'],'Resolution',300,'BackgroundColor','none')    

%% Number for the paper:

Q_outlet_ano_perc = 100*(nanmean(Q_ym_p2,'all')-nanmean(Q_ym_p1,'all'))./nanmean(Q_ym_p1,'all');
Q_gla_ano_perc = 100*(nanmean(Q_gla_ym_p2,'all')-nanmean(Q_gla_ym_p1,'all'))./nanmean(Q_gla_ym_p1,'all');

disp(['Catchment decrease between 2000-2018 and 2018-2023: ' num2str(Q_outlet_ano_perc,2) ' %'])
disp(['Proglacial decrease between 2000-2018 and 2018-2023: ' num2str(Q_gla_ano_perc,2) ' %'])