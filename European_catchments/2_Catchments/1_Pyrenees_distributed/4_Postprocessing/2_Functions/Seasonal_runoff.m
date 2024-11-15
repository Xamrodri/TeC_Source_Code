%%% Calculate monthly mean runoff fluxes for analysis at the catcment scale

% Where to store runoff figures:
dir_fig_runoff = [dir_fig '\Runoff'];
if ~exist(dir_fig_runoff, 'dir'); mkdir(dir_fig_runoff); end 

if ~exist('gla_id_smb','var')
    gla_id_smb = gla_id; % In case no glaciers IDs other than the main glaciers were provided
end 

if ~exist('date_seas','var')
    date_seas = date_m;
end 

% Compute hydrological years starts and ends
period_start = min(find(month(date_seas) == 10 & day(date_seas) == 1,1), find(month(date_seas) == 11 & day(date_seas) == 1,1));
period_end = find(month(date_seas) == 9,1,'last');

Years_no_seas = unique(year(date_seas(period_start:period_end)));

Hydro_year_label = strcat(string(num2str(Years_no_seas(1:end-1)-100*floor(Years_no_seas(1:end-1)/100),'%02d ')), '/', ...
    string(num2str(Years_no_seas(2:end)-100*floor(Years_no_seas(2:end)/100),'%02d ')));

Month_seas = month(date_seas(period_start:period_end));

hydro_month = [10,11,12,1,2,3,4,5,6,7,8,9];
hydro_month_labels = {'Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'};

% Settings elevation bands analysis
dEL=100; % width of elevation bins
ELs = nanmin(DTM,[],'all'):dEL:nanmax(DTM,[],'all');

[Ta_yearly_ym, Pr_ym, Pr_sno_ym, SMB_ym, Imelt_ym, Smelt_ym, Q_ym, ESN_ym,SSN_ym, ET_ym, fsc_MODIS_ym, fsc_MODIS_y_std, ...
    fsc_tc_ym, fsc_tc_y_std, ISO_0_ym, Pr_sno_ym_gla, Imelt_ym_gla] = deal(NaN(numel(Years_no_seas)-1, 12));

date_ym = NaT(numel(Years_no_seas)-1, 12);

[Pr_liq_ym_el, Imelt_ym_el, Smelt_ym_el, SSN_ym_el, ET_ym_el, Pr_ym_el] = deal(NaN(numel(Years_no_seas)-1, 12, numel(ELs)));


%% Calculations

month_m = date_m;%-calmonths(1);

for yy = 1:length(Years_no_seas)-1
    for mm = 1:12
      if hydro_month(mm) > 9
        ind_start = find(year(month_m) == Years_no_seas(yy) & month(month_m) == hydro_month(mm),1,'first');
        ind_end =   find(year(month_m) == Years_no_seas(yy) & month(month_m) == hydro_month(mm),1,'last');

        ind_start_modis = find(year(date_modis) == Years_no_seas(yy) & month(date_modis) == hydro_month(mm),1,'first');
        ind_end_modis =   find(year(date_modis) == Years_no_seas(yy) & month(date_modis) == hydro_month(mm),1,'last');

        ind_start_d = find(year(Date_d) == Years_no_seas(yy) & month(Date_d) == hydro_month(mm),1,'first');
        ind_end_d =   find(year(Date_d) == Years_no_seas(yy) & month(Date_d) == hydro_month(mm),1,'last');

        date_ym(yy,mm) = datetime(Years_no_seas(yy), hydro_month(mm),15);
      else
        ind_start = find(year(month_m) == (Years_no_seas(yy)+1) & month(month_m) == hydro_month(mm),1,'first');
        ind_end =   find(year(month_m) == (Years_no_seas(yy)+1) & month(month_m) == hydro_month(mm),1,'last');

        ind_start_modis = find(year(date_modis) == (Years_no_seas(yy)+1) & month(date_modis) == hydro_month(mm),1,'first');
        ind_end_modis =   find(year(date_modis) == (Years_no_seas(yy)+1) & month(date_modis) == hydro_month(mm),1,'last');

        ind_start_d = find(year(Date_d) == (Years_no_seas(yy)+1) & month(Date_d) == hydro_month(mm),1,'first');
        ind_end_d =   find(year(Date_d) == (Years_no_seas(yy)+1) & month(Date_d) == hydro_month(mm),1,'last');

        date_ym(yy,mm) = datetime(Years_no_seas(yy)+1, hydro_month(mm),15);
      end 

    Ta_y_map = nanmean(TA_map(:,:,ind_start:ind_end),3); Ta_y_map(MASK ~=1) = NaN;
    Ta_yearly_ym(yy,mm) = nanmean(TA_map(:,:,ind_start:ind_end),'all');

    ISO_0_ym(yy,mm) = nanmean(DTM(Ta_y_map>-0.5 & Ta_y_map<0.5));

    Pr_y_map = nansum(PRECIP_map(:,:,ind_start:ind_end),3); Pr_y_map(MASK ~=1) = NaN;
    Pr_sno_y_map = nansum(PSNOW_map(:,:,ind_start:ind_end),3); Pr_sno_y_map(MASK ~=1) = NaN;
    Pr_liq_y_map = Pr_y_map - Pr_sno_y_map;

    Pr_sno_y_map_gla = nansum(PSNOW_map(:,:,ind_start:ind_end),3); Pr_sno_y_map_gla(GLA_ID ~= gla_id_smb) = NaN;
   
    Pr_ym(yy,mm) = nanmean(Pr_y_map,'all'); 
    Pr_sno_ym(yy,mm) = nanmean(Pr_sno_y_map,'all'); 
    Pr_sno_ym_gla(yy,mm) = nanmean(Pr_sno_y_map_gla,'all'); 

    SMB_mapi = nansum(PSNOW_map(:,:,ind_start:ind_end),3) + nansum(AVA_map(:,:,ind_start:ind_end),3) ...
        - nansum(SMSm_map(:,:,ind_start:ind_end),3) - nansum(ESN_map(:,:,ind_start:ind_end),3) ...
        - nansum(SMG_map(:,:,ind_start:ind_end),3) - nansum(EICE_map(:,:,ind_start:ind_end),3);
    SMB_mapi(GLA_ID ~= gla_id_smb) = NaN; % GMB at main glacier
    SMB_ym(yy,mm) = nanmean(SMB_mapi,'all');

    Imelt_y_map = nansum(SMG_map(:,:,ind_start:ind_end),3); Imelt_y_map(MASK~=1) = NaN;
    Imelt_ym(yy,mm) = nanmean(Imelt_y_map,'all');

    Imelt_y_map_gla = nansum(SMG_map(:,:,ind_start:ind_end),3); Imelt_y_map_gla(GLA_ID ~= gla_id_smb) = NaN;
    Imelt_ym_gla(yy,mm) = nanmean(Imelt_y_map_gla,'all');

    Smelt_y_map = nansum(SMSm_map(:,:,ind_start:ind_end),3); Smelt_y_map(MASK ~=1) = NaN;
    Smelt_ym(yy,mm) = nanmean(Smelt_y_map,'all');

    ESN_y_map = nansum(ESN_map(:,:,ind_start:ind_end),3); ESN_y_map(MASK ~=1) = NaN;
    ESN_ym(yy,mm) = nanmean(ESN_y_map,'all');

    SSN_y_map = nansum(SSN_map(:,:,ind_start:ind_end),3); SSN_y_map(MASK ~=1) = NaN;
    SSN_ym(yy,mm) = nanmean(SSN_y_map,'all');     

    ET_y_map = nansum(ET_map(:,:,ind_start:ind_end),3); ET_y_map(MASK ~=1) = NaN;
    ET_ym(yy,mm) = nanmean(ET_y_map,'all');

    ind_ym = (year(date_modis) == Years_no(yy)) & (month(date_modis) == mm);

    fsc_MODIS_ym(yy,mm) = nanmean(fsc_modis_filt_40(ind_start_modis:ind_end_modis));
    fsc_MODIS_y_std(yy,mm) = nanstd(fsc_modis_filt_40(ind_start_modis:ind_end_modis));

    ind_ym = (year(Date_d) == Years_no(yy)) & (month(Date_d) == mm);

    fsc_tc_ym(yy,mm) = nanmean(scas(ind_start_d:ind_end_d));
    fsc_tc_y_std(yy,mm) = nanstd(scas(ind_start_d:ind_end_d));

        for iel = 1:numel(ELs)
              cur=(DTM<(ELs(iel)+dEL/2))&(DTM>=(ELs(iel)-dEL/2)); %current section of DEM
              Pr_liq_ym_el(yy,mm,iel) = nanmean(Pr_liq_y_map(cur));
              Pr_ym_el(yy,mm,iel) = nanmean(Pr_y_map(cur));
              Imelt_ym_el(yy,mm,iel) = nanmean(Imelt_y_map(cur));
              Smelt_ym_el(yy,mm,iel) = nanmean(Smelt_y_map(cur));
              SSN_ym_el(yy,mm,iel) = nanmean(SSN_y_map(cur));
              ET_ym_el(yy,mm,iel) = nanmean(ET_y_map(cur));     
        end 
    end
end 

Pr_liq_ym = Pr_ym - Pr_sno_ym;
TCfsc_yearly_anom = nanmean(fsc_tc_ym,2) - nanmean(nanmean(fsc_tc_ym,2));
Ta_yearly_anom = nanmean(Ta_yearly_ym,2) - nanmean(nanmean(Ta_yearly_ym,2));
Pr_spring_anom = (nansum(Pr_ym(:,5:9),2) - nanmean(nansum(Pr_ym(:,5:9),2)))./nanmean(nansum(Pr_ym(:,5:9),2));

for iel = 1:numel(ELs)
    cur=(DTM<(ELs(iel)+dEL/2))&(DTM>=(ELs(iel)-dEL/2)); %current section of DEM
    Hypso_el(iel) = nansum(cur==1,'all').*sim_res.*sim_res;
    Hypso_gla_el(iel) = nansum(cur==1 & (GLA_ID > 0),'all').*sim_res.*sim_res;
end 

%% Bar precipitation

fi5 = figure('Renderer', 'painters', 'Position',[289 406.3333 350 252.6667]);
b1 = bar(1:12, [nanmean(Pr_sno_ym,1)'  nanmean(Pr_liq_ym,1)'],'stacked'); hold on; 
grid on; ylim([0  1.3.*max([nanmean(Pr_sno_ym,1)'+nanmean(Pr_liq_ym,1)'])])
b1(1,1).FaceColor = [200 100 16]/256;
b1(1,2).FaceColor = [0 0 1];
ylabel('Precipitation [mm]'); legend('Solid','Liquid')
exportgraphics(fi5,[dir_fig '\T&C_precipitation_seasonal_composition.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Seasonal contribution for all different years - STACKED

date_ym_vec = reshape(date_ym',numel(date_ym),1);
Pr_liq_ym_vec = reshape(Pr_liq_ym',numel(Pr_liq_ym),1);
Imelt_ym_vec = reshape(Imelt_ym',numel(Imelt_ym),1);
Smelt_ym_vec = reshape(Smelt_ym',numel(Smelt_ym),1);
SSN_ym_vec = reshape(SSN_ym',numel(SSN_ym),1);
ET_ym_vec = reshape(ET_ym',numel(ET_ym),1);

Net_runoff_ym = Pr_liq_ym + Imelt_ym + Smelt_ym - SSN_ym - ET_ym;
Net_runoff_y = nansum(Net_runoff_ym,2);

if numel(Years_no_seas) < 10

fi5 = figure('Renderer', 'painters', 'Position',[289 127 760 532]);
a1 = area(date_ym_vec, [Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec]); hold on; 
a2 = area(date_ym_vec, [-SSN_ym_vec -ET_ym_vec]);
a1(2).FaceColor = [0.85 0.85 0.85]; a1(3).FaceColor = [0.65 0.95 1]; a1(1).FaceColor = [0 0.6 1];
a1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none'; 
a2(1).FaceColor = [240 100 10]./256; a2(1).FaceAlpha = 0.7;
a2(1).EdgeColor = 'none'; a2(2).EdgeColor = 'none';
plot(date_ym_vec, nansum([Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec -ET_ym_vec -SSN_ym_vec],2),'--k')
ylim([-100 430]); %text(8,400,strcat(Hydro_year_label(chosen_year_ids(ii)),',',' ',string(cluster_label(ii)')),'FontWeight','bold','FontSize',9)
yline(0,'k','HandleVisibility','off','LineWidth',0.8)
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
lg2=legend('Rain','Icemelt','Snowmelt','Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg2,8,"points");  

else

fi2 = figure('Renderer', 'painters', 'Position',[39.6667 157 1.0320e+03 505.3334]);
tiledlayout(2,1,'TileSpacing','compact')
nexttile
a1 = area(date_ym_vec, [Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec]); hold on; 
a2 = area(date_ym_vec, [-SSN_ym_vec -ET_ym_vec]);
a1(2).FaceColor = [0.85 0.85 0.85]; a1(3).FaceColor = [0.65 0.95 1]; a1(1).FaceColor = [0 0.6 1];
a1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none'; 
a2(1).FaceColor = [240 100 10]./256; a2(1).FaceAlpha = 0.7;
a2(1).EdgeColor = 'none'; a2(2).EdgeColor = 'none';
plot(date_ym_vec, nansum([Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec -ET_ym_vec -SSN_ym_vec],2),'--k')
xticks(datetime(year(date_ym_vec(1)),1,1):years:datetime(year(date_ym_vec(end)),1,1))
datetick('x','yyyy','keepticks')
ylim([-100 550]); 
yline(0,'k','HandleVisibility','off','LineWidth',0.8)
xlim([datetime(Years_no_seas(1),7,1) datetime(Years_no_seas(floor(length(Years_no_seas)/2)),7,31)])
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
yyaxis right
scatter(datetime(Years_no_seas(2:end),7,1), Net_runoff_y,30,'+','LineWidth',1.6); ylim([0 1400])
ylabel('Annual runoff [mm]','FontSize',12)
nexttile
a1 = area(date_ym_vec, [Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec]); hold on; 
a2 = area(date_ym_vec, [-SSN_ym_vec -ET_ym_vec]);
a1(2).FaceColor = [0.85 0.85 0.85]; a1(3).FaceColor = [0.65 0.95 1]; a1(1).FaceColor = [0 0.6 1];
a1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none'; 
a2(1).FaceColor = [240 100 10]./256; a2(1).FaceAlpha = 0.7;
a2(1).EdgeColor = 'none'; a2(2).EdgeColor = 'none';
plot(date_ym_vec, nansum([Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec -ET_ym_vec -SSN_ym_vec],2),'--k')
ylim([-100 550]); 
yline(0,'k','HandleVisibility','off','LineWidth',0.8)
xticks(datetime(year(date_ym_vec(1)),1,1):years:datetime(year(date_ym_vec(end)),1,1))
datetick('x','yyyy','keepticks')
xlim([datetime(Years_no_seas(floor(length(Years_no_seas)/2)),7,31) datetime(Years_no_seas(end),7,31)])
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
yyaxis right
scatter(datetime(Years_no_seas(2:end),7,1), Net_runoff_y,30,'+','LineWidth',1.6); ylim([0 1400])
ylabel('Annual runoff [mm]','FontSize',12)
lg2=legend('Rain','Icemelt','Snowmelt','Sublimation','Evapotranspiration','Net','Annual Sum','location','northeast'); fontsize(lg2,8,"points");  
lg2.NumColumns = 4;

exportgraphics(fi2,[dir_fig_runoff '\Year_all_monthly_runoff_springP.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
end 

%% Seasonal contribution for all different years - LINES

if numel(Years_no_seas) < 10

fi5 = figure('Renderer', 'painters', 'Position',[289 127 760 532]);
pl1 = plot(date_ym_vec, [Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec],'LineWidth',2); hold on; grid on;
pl2 = plot(date_ym_vec, [-SSN_ym_vec -ET_ym_vec],'LineWidth',2);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 0.9];
% pl1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none'; 
pl2(1).Color = [240 100 10 180]./256; 
% pl2(1).EdgeColor = 'none'; a2(2).EdgeColor = 'none';
plot(date_ym_vec, nansum([Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec -ET_ym_vec -SSN_ym_vec],2),'--k')
ylim([-100 430]); %text(8,400,strcat(Hydro_year_label(chosen_year_ids(ii)),',',' ',string(cluster_label(ii)')),'FontWeight','bold','FontSize',9)
xticks(date_ym_vec(1):years:date_ym_vec(end))
datetick('x','yyyy','keepticks')
yline(0,'k','HandleVisibility','off','LineWidth',0.8)
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
lg2=legend('Rain','Icemelt','Snowmelt','Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg2,8,"points");  

else

fi2 = figure('Renderer', 'painters', 'Position',[39.6667 157 1.0320e+03 505.3334]);
tiledlayout(2,1,'TileSpacing','compact')
nexttile
pl1 = plot(date_ym_vec, [Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec],'LineWidth',2); hold on; grid on;
pl2 = plot(date_ym_vec, [-SSN_ym_vec -ET_ym_vec],'LineWidth',1.3);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
% pl1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none'; 
pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
% pl2(1).EdgeColor = 'none'; a2(2).EdgeColor = 'none';
plot(date_ym_vec, nansum([Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec -ET_ym_vec -SSN_ym_vec],2),'--k','LineWidth',0.7)
ylim([-100 400]);
xticks(datetime(year(date_ym_vec(1)),1,1):years:datetime(year(date_ym_vec(end)),1,1))
datetick('x','yyyy','keepticks')
yline(0,'k','HandleVisibility','off','LineWidth',0.8)
xlim([datetime(Years_no_seas(1),7,1) datetime(Years_no_seas(floor(length(Years_no_seas)/2)),7,31)])
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
yyaxis right
scatter(datetime(Years_no_seas(2:end),7,1), Net_runoff_y,30,'+','LineWidth',1.6); ylim([0 1400])
ylabel('Annual runoff [mm]','FontSize',12)
nexttile
pl1 = plot(date_ym_vec, [Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec],'LineWidth',2); hold on; grid on;
pl2 = plot(date_ym_vec, [-SSN_ym_vec -ET_ym_vec],'LineWidth',1.3);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
% pl1(1).EdgeColor = 'none'; a1(2).EdgeColor = 'none'; a1(3).EdgeColor = 'none'; 
pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
% pl2(1).EdgeColor = 'none'; a2(2).EdgeColor = 'none';
plot(date_ym_vec, nansum([Pr_liq_ym_vec  Imelt_ym_vec Smelt_ym_vec -ET_ym_vec -SSN_ym_vec],2),'--k','LineWidth',0.7)
ylim([-100 400]); 
yline(0,'k','HandleVisibility','off','LineWidth',0.8)
xticks(datetime(year(date_ym_vec(1)),1,1):years:datetime(year(date_ym_vec(end)),1,1))
datetick('x','yyyy','keepticks')
xlim([datetime(Years_no_seas(floor(length(Years_no_seas)/2)),7,31) datetime(Years_no_seas(end),7,31)])
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
yyaxis right
scatter(datetime(Years_no_seas(2:end),7,1), Net_runoff_y,30,'+','LineWidth',1.6); ylim([0 1400])
ylabel('Annual runoff [mm]','FontSize',12)
lg2=legend('Rain','Icemelt','Snowmelt','Sublimation','Evapotranspiration','Net','Annual Sum','location','northeast'); fontsize(lg2,8,"points");  
lg2.NumColumns = 4;

exportgraphics(fi2,[dir_fig_runoff '\Year_all_monthly_runoff_lines.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
end 

%% Monthly contribution difference between two different period

if ~exist('date_start_runoff_p1','var') || ~exist('date_end_runoff_p1','var')
    date_start_runoff_p1 = date_m(1);
    date_end_runoff_p1 = date_m(floor(length(date_m)/2));
end 

if ~exist('date_start_runoff_p2','var') || ~exist('date_end_runoff_p2','var')
    date_start_runoff_p2 = date_m(1+floor(length(date_m)/2));
    date_end_runoff_p2 = date_m(end);
end 

ind_p1 = (date_ym >= date_start_runoff_p1) & (date_ym <= date_end_runoff_p1);
ind_p2 = (date_ym >= date_start_runoff_p2) & (date_ym <= date_end_runoff_p2);

Pr_liq_ym_p1 = Pr_liq_ym; Pr_liq_ym_p1(~ind_p1) = NaN;
Pr_ym_p1 = Pr_ym; Pr_ym_p1(~ind_p1) = NaN;
Imelt_ym_p1 = Imelt_ym; Imelt_ym_p1(~ind_p1) = NaN;
Smelt_ym_p1 = Smelt_ym; Smelt_ym_p1(~ind_p1) = NaN;
SSN_ym_p1 = SSN_ym; SSN_ym_p1(~ind_p1) = NaN;
ET_ym_p1 = ET_ym; ET_ym_p1(~ind_p1) = NaN;

Pr_liq_ym_p2 = Pr_liq_ym; Pr_liq_ym_p2(~ind_p2) = NaN;
Pr_ym_p2 = Pr_ym; Pr_ym_p2(~ind_p2) = NaN;
Imelt_ym_p2 = Imelt_ym; Imelt_ym_p2(~ind_p2) = NaN;
Smelt_ym_p2 = Smelt_ym; Smelt_ym_p2(~ind_p2) = NaN;
SSN_ym_p2 = SSN_ym; SSN_ym_p2(~ind_p2) = NaN;
ET_ym_p2 = ET_ym; ET_ym_p2(~ind_p2) = NaN;

% Compute changes in runoff contribution between the two periods

Icemelt_contrib_p1 = nanmean(nansum(Imelt_ym_p1,2)./(nansum(Imelt_ym_p1,2)+nansum(Smelt_ym_p1,2) + nansum(Pr_liq_ym_p1,2)));
Icemelt_contrib_p2 = nanmean(nansum(Imelt_ym_p2,2)./(nansum(Imelt_ym_p2,2)+nansum(Smelt_ym_p2,2) + nansum(Pr_liq_ym_p2,2)));

Snowmelt_contrib_p1 = nanmean(nansum(Smelt_ym_p1,2)./(nansum(Imelt_ym_p1,2)+nansum(Smelt_ym_p1,2) + nansum(Pr_liq_ym_p1,2)));
Snowmelt_contrib_p2 = nanmean(nansum(Smelt_ym_p2,2)./(nansum(Imelt_ym_p2,2)+nansum(Smelt_ym_p2,2) + nansum(Pr_liq_ym_p2,2)));

Rain_contrib_p1 = nanmean(nansum(Pr_liq_ym_p1,2)./(nansum(Imelt_ym_p1,2)+nansum(Smelt_ym_p1,2) + nansum(Pr_liq_ym_p1,2)));
Rain_contrib_p2 = nanmean(nansum(Pr_liq_ym_p2,2)./(nansum(Imelt_ym_p2,2)+nansum(Smelt_ym_p2,2) + nansum(Pr_liq_ym_p2,2)));

Runoff_deficit = nansum(Imelt_ym_p2 + Smelt_ym_p2 + Pr_liq_ym_p2 - SSN_ym_p2 - ET_ym_p2 ,'all')./nYears_p2 - nansum(Imelt_ym_p1 + Smelt_ym_p1 + Pr_liq_ym_p1 - SSN_ym_p1 - ET_ym_p1 ,'all')./nYears_p1;
Runoff_woIcemelt = nansum(Smelt_ym_p2 + Pr_liq_ym_p2 - SSN_ym_p2 - ET_ym_p2 ,'all')./nYears_p2 - nansum(Smelt_ym_p1 + Pr_liq_ym_p1 - SSN_ym_p1 - ET_ym_p1 ,'all')./nYears_p1;
Runoff_woET = nansum(Imelt_ym_p2 + Smelt_ym_p2 + Pr_liq_ym_p2 - SSN_ym_p2 ,'all')./nYears_p2 - nansum(Imelt_ym_p1 + Smelt_ym_p1 + Pr_liq_ym_p1 - SSN_ym_p1 ,'all')./nYears_p1;
Icemelt_change = nansum(Imelt_ym_p2,'all')./nYears_p2 - nansum(Imelt_ym_p1,'all')./nYears_p1;
ET_change = nansum(ET_ym_p2,'all')./nYears_p2 - nansum(ET_ym_p1,'all')./nYears_p1;
Smelt_change = nansum(Smelt_ym_p2,'all')./nYears_p2 - nansum(Smelt_ym_p1,'all')./nYears_p1;
Rain_change = nansum(Pr_liq_ym_p2,'all')./nYears_p2 - nansum(Pr_liq_ym_p1,'all')./nYears_p1;
Precip_change = nansum(Pr_ym_p2,'all')./nYears_p2 - nansum(Pr_ym_p1,'all')./nYears_p1;
Precip_change_perc = 100*(nansum(Pr_ym_p2,'all')./nYears_p2 - nansum(Pr_ym_p1,'all')./nYears_p1)/(nansum(Pr_ym_p1,'all')./nYears_p1);
Snowfall_change = nansum(Pr_ym_p2-Pr_liq_ym_p2,'all')./nYears_p2 - nansum(Pr_ym_p1-Pr_liq_ym_p1,'all')./nYears_p1;
Snowfall_change_perc = 100*(nansum(Pr_ym_p2-Pr_liq_ym_p2,'all')./nYears_p2 - nansum(Pr_ym_p1-Pr_liq_ym_p1,'all')./nYears_p1)./(nansum(Pr_ym_p1-Pr_liq_ym_p1,'all')./nYears_p1);
Sublimation_change =  nansum(SSN_ym_p2,'all')./nYears_p2 - nansum(SSN_ym_p1,'all')./nYears_p1;

Compensation_icemelt_deficit = Icemelt_change./abs(Runoff_woIcemelt);
Compensation_ET_deficit = (Runoff_woET-Runoff_deficit)./abs(Runoff_woET);

% Compute averages over the whole period

Precip_period_average = nansum(Pr_ym,'all')./nYears;
Snowfall_period_average = nansum(Pr_sno_ym,'all')./nYears;
SnowfallFraction_period_average = Snowfall_period_average./Precip_period_average;

SnowfallFraction_period_1 = nansum(Pr_ym_p1-Pr_liq_ym_p1,'all')./nansum(Pr_ym_p1,'all');
SnowfallFraction_period_2 = nansum(Pr_ym_p2-Pr_liq_ym_p2,'all')./nansum(Pr_ym_p2,'all');

SnowfallFraction_Oct_Dec_period_1 = nansum(Pr_ym_p1(:,1:3)-Pr_liq_ym_p1(:,1:3),'all')./nansum(Pr_ym_p1(:,1:3),'all');
SnowfallFraction_Oct_Dec_period_2 = nansum(Pr_ym_p2(:,1:3)-Pr_liq_ym_p2(:,1:3),'all')./nansum(Pr_ym_p2(:,1:3),'all');

SnowfallFraction_Jan_Mar_period_1 = nansum(Pr_ym_p1(:,4:6)-Pr_liq_ym_p1(:,4:6),'all')./nansum(Pr_ym_p1(:,4:6),'all');
SnowfallFraction_Jan_Mar_period_2 = nansum(Pr_ym_p2(:,4:6)-Pr_liq_ym_p2(:,4:6),'all')./nansum(Pr_ym_p2(:,4:6),'all');

SnowfallFraction_Apr_Jun_period_1 = nansum(Pr_ym_p1(:,7:9)-Pr_liq_ym_p1(:,7:9),'all')./nansum(Pr_ym_p1(:,7:9),'all');
SnowfallFraction_Apr_Jun_period_2 = nansum(Pr_ym_p2(:,7:9)-Pr_liq_ym_p2(:,7:9),'all')./nansum(Pr_ym_p2(:,7:9),'all');

SnowfallFraction_Jul_Sep_period_1 = nansum(Pr_ym_p1(:,10:12)-Pr_liq_ym_p1(:,10:12),'all')./nansum(Pr_ym_p1(:,10:12),'all');
SnowfallFraction_Jul_Sep_period_2 = nansum(Pr_ym_p2(:,10:12)-Pr_liq_ym_p2(:,10:12),'all')./nansum(Pr_ym_p2(:,10:12),'all');

disp(['Changes in catchment ET between the two periods: ' num2str(ET_change,1) ' mm w.e., (' num2str(100*ET_change./(nansum(ET_ym_p1,'all')./nYears_p1),1) ...
    ' %)'])

disp(['Changes in catchment rain between the two periods: ' num2str(Rain_change,3) ' mm w.e., (' num2str(100*Rain_change./(nansum(Pr_liq_ym_p1,'all')./nYears_p1),2) ...
    ' %)'])


%% Final figure, hydrograph monthly fluxes between two periods
fi2 = figure('Renderer', 'painters', 'Position',[131.6667 304.3333 798 333.3333]);
tiledlayout(1,2,'TileSpacing','tight')
nexttile
pl1 = plot(1:12, [nanmean(Pr_liq_ym_p1,1)'  nanmean(Imelt_ym_p1,1)' nanmean(Smelt_ym_p1,1)'],'LineWidth',2.5); hold on; 
pl2 = plot(1:12, [-nanmean(SSN_ym_p1,1)' -nanmean(ET_ym_p1,1)'],'LineWidth',2.5);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 180]./256; 
plot(1:12, nanmean(nanmean(Pr_liq_ym_p1,1)' + nanmean(Imelt_ym_p1,1)'+ nanmean(Smelt_ym_p1,1)' ...
    -nanmean(SSN_ym_p1,1)' -nanmean(ET_ym_p1,1)',2),'--k','LineWidth',1.2)
xticks([1:2:12]); xlim([1 12]); ylim([-100 300]); text(8,270,[num2str(year(date_start_runoff_p1)) ' - ' num2str(year(date_end_runoff_p1))],'FontWeight','bold','FontSize',10)
set(gca,'XTickLabels',hydro_month_labels(1:2:12)); yline(0,'k','HandleVisibility','off','LineWidth',0.8)
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
lg2=legend(['Rain (' num2str(Rain_contrib_p1*100,2) '%)'],['Icemelt (' num2str(Icemelt_contrib_p1*100,2) '%)'],...
    ['Snowmelt (' num2str(Snowmelt_contrib_p1*100,2) '%)'],'Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg2,8,"points");  
nexttile
pl1 = plot(1:12, [nanmean(Pr_liq_ym_p2,1)'  nanmean(Imelt_ym_p2,1)' nanmean(Smelt_ym_p2,1)'],'LineWidth',2.5); hold on; 
pl2 = plot(1:12, [-nanmean(SSN_ym_p2,1)' -nanmean(ET_ym_p2,1)'],'LineWidth',2.5);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 180]./256; 
plot(1:12, nanmean(nanmean(Pr_liq_ym_p2,1)' + nanmean(Imelt_ym_p2,1)'+ nanmean(Smelt_ym_p2,1)' ...
    -nanmean(SSN_ym_p2,1)' -nanmean(ET_ym_p2,1)',2),'--k','LineWidth',1.2)
xticks([1:2:12]); xlim([1 12]); ylim([-100 300]); text(8,270,[num2str(year(date_start_runoff_p2)) ' - ' num2str(year(date_end_runoff_p2))],'FontWeight','bold','FontSize',10)
set(gca,'XTickLabels',hydro_month_labels(1:2:12)); yline(0,'k','HandleVisibility','off','LineWidth',0.8)
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
set(gca,'YAxisLocation','right')
lg2=legend(['Rain (' num2str(Rain_contrib_p2*100,2) '%)'],['Icemelt (' num2str(Icemelt_contrib_p2*100,2) '%)'],...
    ['Snowmelt (' num2str(Snowmelt_contrib_p2*100,2) '%)'],'Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg2,8,"points");  
exportgraphics(fi2,[dir_fig_runoff '\Seasonal_runoff_2periods_lines.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Compute mean elevation fluxes for two different periods

Year_start = 2000;
Year_middle = 2018;
Year_end = 2023;

ind_p1 = Years_no_seas(1:end-1) < Year_middle;
ind_p2 = Years_no_seas(1:end-1) >= Year_middle;

Pr_liq_y_el_p1 = squeeze(nansum(Pr_liq_ym_el,2)); Pr_liq_y_el_p1 = nanmean(Pr_liq_y_el_p1(ind_p1,:),1)';
Imelt_y_el_p1 = squeeze(nansum(Imelt_ym_el,2)); Imelt_y_el_p1 = nanmean(Imelt_y_el_p1(ind_p1,:),1)';
Smelt_y_el_p1 = squeeze(nansum(Smelt_ym_el,2)); Smelt_y_el_p1 = nanmean(Smelt_y_el_p1(ind_p1,:),1)';
SSN_y_el_p1 = squeeze(nansum(SSN_ym_el,2)); SSN_y_el_p1 = nanmean(SSN_y_el_p1(ind_p1,:),1)';
ET_y_el_p1 = squeeze(nansum(ET_ym_el,2)); ET_y_el_p1 = nanmean(ET_y_el_p1(ind_p1,:),1)';

Pr_liq_y_el_p2 = squeeze(nansum(Pr_liq_ym_el,2)); Pr_liq_y_el_p2 = nanmean(Pr_liq_y_el_p2(ind_p2,:),1)';
Imelt_y_el_p2 = squeeze(nansum(Imelt_ym_el,2)); Imelt_y_el_p2 = nanmean(Imelt_y_el_p2(ind_p2,:),1)';
Smelt_y_el_p2 = squeeze(nansum(Smelt_ym_el,2)); Smelt_y_el_p2 = nanmean(Smelt_y_el_p2(ind_p2,:),1)';
SSN_y_el_p2 = squeeze(nansum(SSN_ym_el,2)); SSN_y_el_p2 = nanmean(SSN_y_el_p2(ind_p2,:),1)';
ET_y_el_p2 = squeeze(nansum(ET_ym_el,2)); ET_y_el_p2 = nanmean(ET_y_el_p2(ind_p2,:),1)';

%% Prepare catchment hypsometry figure

%%% 1 Fir  (evergr.)
%%% 2 Larch (decid.)
%%% 3 Grass C3
%%% 4 Shrub (decid.)
%%% 5 Broadleaf evergreen
%%% 6 Broadleaf deciduous
%%% 7 Rock
%%% 8 Ice
%%% 9 Debris-coverd ice

lc_labels = {'Fir','Larch','Grass','Shrub','Broadlf evergreen','Broadleaf deci.','Rock','Clean ice',sprintf('Debris-covered\\newline          ice')};

colors = [[255, 90, 1]./255; [0.9294    0.6941    0.1255] ; [90, 182, 20]./255; ...
          [20, 107, 24]./255;  [0.3 0.7 0.6]  ; [0.3 0.8 0.7]; ...
          [96 95 95]./255 ; [121, 208, 235]./255  ;  [193, 182, 122]./255];

VEG_show = VEG_CODE; VEG_show(isnan(DTM)) = NaN;
VEG_show(VEG_show == 8 & DEB_MAP > 0) = max(VEG_CODE(:)) + 1;
VEG_map = VEG_show;

% area_tot = nansum(~isnan(DTM),'all').*100*100*10^-6;  % 161.9 km^2

lc_el = [];

for lci = 1:length(colors)
for iel = 1:numel(ELs)
    cur=(DTM<(ELs(iel)+dEL/2))&(DTM>=(ELs(iel)-dEL/2)); %current section of DEM
    lc_el(iel,lci) = nansum(VEG_map(cur) == lci).*100*100*10^-6;
end
end 

for ii = 1:length(lc_labels)
    if nansum(lc_el(:,ii)) == 0
        lc_labels{ii} = '';
    end 
end 

%% Final figure 6, without hypsometry

% fi2 = figure('Renderer', 'painters', 'Position',[50 50 1.0474e+03 576.6667]);
% tiledlayout(2,5,'TileSpacing','compact')
% 
% nexttile(1,[1 2])
% pl1 = plot(1:12, [nanmean(Pr_liq_ym_p1,1)'  nanmean(Imelt_ym_p1,1)' nanmean(Smelt_ym_p1,1)'],'LineWidth',2.5); hold on; 
% pl2 = plot(1:12, [-nanmean(SSN_ym_p1,1)' -nanmean(ET_ym_p1,1)'],'LineWidth',2.5);
% pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
% pl2(1).Color = [240 100 10 180]./256; 
% plot(1:12, nanmean(nanmean(Pr_liq_ym_p1,1)' + nanmean(Imelt_ym_p1,1)'+ nanmean(Smelt_ym_p1,1)' ...
%     -nanmean(SSN_ym_p1,1)' -nanmean(ET_ym_p1,1)',2),'--k','LineWidth',1.2)
% xticks([1:1:12]); xlim([1 12]); ylim([-100 300]); text(8,270,[num2str(year(date_start_runoff_p1)) ' - ' num2str(year(date_end_runoff_p1))],'FontWeight','bold','FontSize',10)
% set(gca,'XTickLabels',hydro_month_labels(1:1:12)); yline(0,'k','HandleVisibility','off','LineWidth',0.8)
% ylabel('Runoff contribution (mm)','FontSize',12); grid on;
% lg2=legend(['Rain (' num2str(Rain_contrib_p1*100,2) '%)'],['Icemelt (' num2str(Icemelt_contrib_p1*100,2) '%)'],...
%     ['Snowmelt (' num2str(Snowmelt_contrib_p1*100,2) '%)'],'Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg2,8,"points");  
% 
% nexttile(6,[1 2])
% pl1 = plot(1:12, [nanmean(Pr_liq_ym_p2,1)'  nanmean(Imelt_ym_p2,1)' nanmean(Smelt_ym_p2,1)'],'LineWidth',2.5); hold on; 
% pl2 = plot(1:12, [-nanmean(SSN_ym_p2,1)' -nanmean(ET_ym_p2,1)'],'LineWidth',2.5);
% pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
% pl2(1).Color = [240 100 10 180]./256; 
% plot(1:12, nanmean(nanmean(Pr_liq_ym_p2,1)' + nanmean(Imelt_ym_p2,1)'+ nanmean(Smelt_ym_p2,1)' ...
%     -nanmean(SSN_ym_p2,1)' -nanmean(ET_ym_p2,1)',2),'--k','LineWidth',1.2)
% xticks([1:1:12]); xlim([1 12]); ylim([-100 300]); text(8,270,[num2str(year(date_start_runoff_p2)) ' - ' num2str(year(date_end_runoff_p2))],'FontWeight','bold','FontSize',10)
% set(gca,'XTickLabels',hydro_month_labels(1:1:12)); yline(0,'k','HandleVisibility','off','LineWidth',0.8)
% ylabel('Runoff contribution (mm)','FontSize',12); grid on;
% set(gca,'YAxisLocation','left')
% lg2=legend(['Rain (' num2str(Rain_contrib_p2*100,2) '%)'],['Icemelt (' num2str(Icemelt_contrib_p2*100,2) '%)'],...
%     ['Snowmelt (' num2str(Snowmelt_contrib_p2*100,2) '%)'],'Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg2,8,"points");  
% 
% nexttile(3,[2 1])
% pl1 = plot(ELs, [Pr_liq_y_el_p1  Imelt_y_el_p1 Smelt_y_el_p1].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
% pl2 = plot(ELs, [-SSN_y_el_p1 -ET_y_el_p1].*c_fact'.*10^-6,'LineWidth',2.5);
% pl3 = plot(ELs, (Pr_liq_y_el_p1+Imelt_y_el_p1+Smelt_y_el_p1-SSN_y_el_p1-ET_y_el_p1).*c_fact'.*10^-6,'--k','LineWidth',1.1);
% pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
% pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
% yline(0,'Color',[0.2 0.2 0.2])
% text(5500,2,[num2str(year(date_start_runoff_p1)) ' - ' num2str(year(date_end_runoff_p1))],'FontWeight','bold','FontSize',10)
% area(ELs, -12+ 60*(Hypso_el./nansum(Hypso_el)), 'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeColor','none','BaseValue',-12)
% area(ELs, -12+ 60*(Hypso_gla_el./nansum(Hypso_el)), 'FaceColor',[121, 208, 235]./255,'FaceAlpha',0.5,'EdgeColor','none','BaseValue',-12)
% view(90,-90);grid on; ylim([-12 17])
% ylabel('m^3 w.e. (*10^6)'); 
% set(gca,'XTickLabel',[]); xlim([min(ELs) max(ELs)])
% 
% nexttile(4,[2 1])
% pl1 = plot(ELs, [Pr_liq_y_el_p2  Imelt_y_el_p2 Smelt_y_el_p2].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
% pl2 = plot(ELs, [-SSN_y_el_p2 -ET_y_el_p2].*c_fact'.*10^-6,'LineWidth',2.5);
% pl3 = plot(ELs, (Pr_liq_y_el_p2+Imelt_y_el_p2+Smelt_y_el_p2-SSN_y_el_p2-ET_y_el_p2).*c_fact'.*10^-6,'--k','LineWidth',1.1);
% pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
% pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
% yline(0,'Color',[0.2 0.2 0.2])
% text(5500,2,[num2str(year(date_start_runoff_p2)) ' - ' num2str(year(date_end_runoff_p2))],'FontWeight','bold','FontSize',10)
% view(90,-90);grid on; ylim([-12 17]); xlim([min(ELs) max(ELs)])
% %xlabel('Elevation [m a.s.l.]','FontSize',12); 
% ylabel('m^3 w.e. (*10^6)'); 
% set(gca,'XAxisLocation','top'); set(gca,'XTickLabel',[]);
% 
% nexttile(5,[2 1])
% pl1 = plot(ELs, [Pr_liq_y_el_p2-Pr_liq_y_el_p1  Imelt_y_el_p2-Imelt_y_el_p1 Smelt_y_el_p2-Smelt_y_el_p1].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
% pl2 = plot(ELs, [-SSN_y_el_p2+SSN_y_el_p1 -ET_y_el_p2+ET_y_el_p2].*c_fact'.*10^-6,'LineWidth',2.5);
% pl3 = plot(ELs, (Pr_liq_y_el_p2+Imelt_y_el_p2+Smelt_y_el_p2-SSN_y_el_p2-ET_y_el_p2 - (Pr_liq_y_el_p1+Imelt_y_el_p1+Smelt_y_el_p1-SSN_y_el_p1-ET_y_el_p1)).*c_fact'.*10^-6,'--k','LineWidth',1.1);
% pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
% pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
% yline(0,'Color',[0.2 0.2 0.2])
% text(5500,1,'\DeltaPeriods','FontWeight','bold','FontSize',10)
% view(90,-90);grid on; ylim([-10 10]); xlim([min(ELs) max(ELs)])
% xlabel('Elevation [m a.s.l.]','FontSize',12); ylabel('m^3 w.e. (*10^6)'); 
% set(gca,'XAxisLocation','top'); 
% 
% exportgraphics(fi2,[dir_fig_runoff '\Elevation_fluxes_elev_2periods_lines.png'],'Resolution',300,'BackgroundColor','none')
% % close(gcf)

%% Final figure 6, with hypsometry

if strcmp(glacier,'Langtang')
    ylim_vol = [-20 80];%[-12 17];
    ylim_vol_diff = [-40 10];
    ymax_area = 40; %20;
    ylim_runMo = [-100 400]; %[-100 300];
    yaxis_text_elev = 6500;%5500; 
elseif strcmp(glacier,'Kyzylsu')
    ylim_vol = [-12 17];
    ylim_vol_diff = [-10 10];
    ymax_area = 20;
    ylim_runMo = [-100 300];
    yaxis_text_elev = 5500;
elseif strcmp(glacier,'Parlung4')
    ylim_vol = [-12 17];
    ylim_vol_diff = [-10 10];
    ymax_area = 20;
    ylim_runMo = [-100 600];
    yaxis_text_elev = 5800;
end 

fi2 = figure('Renderer', 'painters', 'Position',[50 50 885.3667 784.6667]);
tiledlayout(3,4,'TileSpacing','tight')

nexttile(1,[2 1])
pl1 = plot(ELs, [Pr_liq_y_el_p1  Imelt_y_el_p1 Smelt_y_el_p1].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
pl2 = plot(ELs, [-SSN_y_el_p1 -ET_y_el_p1].*c_fact'.*10^-6,'LineWidth',2.5);
pl3 = plot(ELs, (Pr_liq_y_el_p1+Imelt_y_el_p1+Smelt_y_el_p1-SSN_y_el_p1-ET_y_el_p1).*c_fact'.*10^-6,'--k','LineWidth',1.1);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
yline(0,'Color',[0.2 0.2 0.2])
text(yaxis_text_elev,2,[num2str(year(date_start_runoff_p1)) ' - ' num2str(year(date_end_runoff_p1))],'FontWeight','bold','FontSize',10)
view(90,-90);grid on; ylim(ylim_vol)
ylabel('m^3 w.e. (*10^6)'); 
xlim([min(ELs) max(ELs)])
xlabel('Elevation [m a.s.l.]','FontSize',12); 
lg4 = legend('Rain ','Icemelt ','Snowmelt', 'Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg4,9,"points");  
lg4.NumColumns = 6; lg4.Position = [0.095    0.936    0.6860    0.0246];

nexttile(2,[2 1])
pl1 = plot(ELs, [Pr_liq_y_el_p2  Imelt_y_el_p2 Smelt_y_el_p2].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
pl2 = plot(ELs, [-SSN_y_el_p2 -ET_y_el_p2].*c_fact'.*10^-6,'LineWidth',2.5);
pl3 = plot(ELs, (Pr_liq_y_el_p2+Imelt_y_el_p2+Smelt_y_el_p2-SSN_y_el_p2-ET_y_el_p2).*c_fact'.*10^-6,'--k','LineWidth',1.1);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
yline(0,'Color',[0.2 0.2 0.2])
text(yaxis_text_elev,2,[num2str(year(date_start_runoff_p2)) ' - ' num2str(year(date_end_runoff_p2))],'FontWeight','bold','FontSize',10)
view(90,-90);grid on; ylim(ylim_vol); xlim([min(ELs) max(ELs)])
%xlabel('Elevation [m a.s.l.]','FontSize',12); 
ylabel('m^3 w.e. (*10^6)'); 
set(gca,'XAxisLocation','top'); set(gca,'XTickLabel',[]);

nexttile(3,[2 1])
pl1 = plot(ELs, [Pr_liq_y_el_p2-Pr_liq_y_el_p1  Imelt_y_el_p2-Imelt_y_el_p1 Smelt_y_el_p2-Smelt_y_el_p1].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
pl2 = plot(ELs, [-SSN_y_el_p1+SSN_y_el_p2 -ET_y_el_p1+ET_y_el_p2].*c_fact'.*10^-6,'LineWidth',1.5);
pl3 = plot(ELs, (Pr_liq_y_el_p2+Imelt_y_el_p2+Smelt_y_el_p2-SSN_y_el_p2-ET_y_el_p2 - (Pr_liq_y_el_p1+Imelt_y_el_p1+Smelt_y_el_p1-SSN_y_el_p1-ET_y_el_p1)).*c_fact'.*10^-6,'--k','LineWidth',1.1);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
yline(0,'Color',[0.2 0.2 0.2])
text(yaxis_text_elev,1,'\DeltaPeriods','FontWeight','bold','FontSize',10)
view(90,-90);grid on; ylim(ylim_vol_diff); xlim([min(ELs) max(ELs)])
ylabel('m^3 w.e. (*10^6)'); set(gca,'XTickLabel',[]);
set(gca,'XAxisLocation','top'); 

nexttile(4, [2 1])
b1 = area(ELs,lc_el,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
for i = 1:length(b1)
    set(b1(i), 'FaceColor', colors(i,:)) 
end
ylabel('Area [km^2]'); grid on;
xlim([min(ELs) max(ELs)]); view(-90,90); ylim([0 ymax_area]); set(gca,'XAxisLocation','top')
lg3 = legend(lc_labels,'FontSize',8.5,'Box','off','Location','NorthEast');
lg3.Position = [0.728    0.780    0.1502    0.1468];

xlabel('Elevation [m a.s.l.]','FontSize',12); 

nexttile(9,[1 2])
pl1 = plot(1:12, [nanmean(Pr_liq_ym_p1,1)'  nanmean(Imelt_ym_p1,1)' nanmean(Smelt_ym_p1,1)'],'LineWidth',2.5); hold on; 
pl2 = plot(1:12, [-nanmean(SSN_ym_p1,1)' -nanmean(ET_ym_p1,1)'],'LineWidth',2.5);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 180]./256; 
plot(1:12, nanmean(nanmean(Pr_liq_ym_p1,1)' + nanmean(Imelt_ym_p1,1)'+ nanmean(Smelt_ym_p1,1)' ...
    -nanmean(SSN_ym_p1,1)' -nanmean(ET_ym_p1,1)',2),'--k','LineWidth',1.2)
xticks([1:1:12]); xlim([1 12]); ylim(ylim_runMo); text(8,270,[num2str(year(date_start_runoff_p1)) ' - ' num2str(year(date_end_runoff_p1))],'FontWeight','bold','FontSize',10)
set(gca,'XTickLabels',hydro_month_labels(1:1:12)); yline(0,'k','HandleVisibility','off','LineWidth',0.8)
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
lg2=legend(['Rain (' num2str(Rain_contrib_p1*100,2) '%)'],['Icemelt (' num2str(Icemelt_contrib_p1*100,2) '%)'],...
    ['Snowmelt (' num2str(Snowmelt_contrib_p1*100,2) '%)'],'Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg2,8,"points");  

nexttile(11,[1 2])
pl1 = plot(1:12, [nanmean(Pr_liq_ym_p2,1)'  nanmean(Imelt_ym_p2,1)' nanmean(Smelt_ym_p2,1)'],'LineWidth',2.5); hold on; 
pl2 = plot(1:12, [-nanmean(SSN_ym_p2,1)' -nanmean(ET_ym_p2,1)'],'LineWidth',2.5);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 180]./256; 
plot(1:12, nanmean(nanmean(Pr_liq_ym_p2,1)' + nanmean(Imelt_ym_p2,1)'+ nanmean(Smelt_ym_p2,1)' ...
    -nanmean(SSN_ym_p2,1)' -nanmean(ET_ym_p2,1)',2),'--k','LineWidth',1.2)
xticks([1:1:12]); xlim([1 12]); ylim([-100 300]); text(8,270,[num2str(year(date_start_runoff_p2)) ' - ' num2str(year(date_end_runoff_p2))],'FontWeight','bold','FontSize',10)
set(gca,'XTickLabels',hydro_month_labels(1:1:12)); yline(0,'k','HandleVisibility','off','LineWidth',0.8)
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
set(gca,'YAxisLocation','right')
lg2=legend(['Rain (' num2str(Rain_contrib_p2*100,2) '%)'],['Icemelt (' num2str(Icemelt_contrib_p2*100,2) '%)'],...
    ['Snowmelt (' num2str(Snowmelt_contrib_p2*100,2) '%)'],'Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg2,8,"points");  

exportgraphics(fi2,[dir_fig_runoff '\Elevation_fluxes_elev_2periods_lines_hypso.png'],'Resolution',300,'BackgroundColor','none')
% close(gcf)

% Save seasonal runoff componts

Runoff.Pr_liq_ym_el = Pr_liq_ym_el;
Runoff.Imelt_ym_el = Imelt_ym_el;
Runoff.Smelt_ym_el =Smelt_ym_el;
Runoff.SSN_ym_el =SSN_ym_el;
Runoff.ET_ym_el = ET_ym_el;
Runoff.ELs  = ELs';

Runoff.Pr_liq_ym = Pr_liq_ym;
Runoff.Imelt_ym = Imelt_ym;
Runoff.Smelt_ym =Smelt_ym;
Runoff.SSN_ym =SSN_ym;
Runoff.ET_ym = ET_ym;
Runoff.date_ym = date_ym;

save([dir_fig_runoff '\Runoff_components_' glacier '_monthly_' ...
    datestr(date_ym(1),'yyyy-mm') '_' datestr(date_ym(end),'yyyy-mm') '.mat'],'Runoff');

% Analyze changes in evapotranspiratio
figure
pl1 = plot(1:12, nanmean(ET_ym_p1,1)','LineWidth',2.5); grid on; hold on;
pl2 = plot(1:12,  nanmean(ET_ym_p2,1)','LineWidth',2.5); grid on; hold on;
xticks([1:1:12]); xlim([1 12]);
set(gca,'XTickLabels',hydro_month_labels(1:1:12)); 
ylabel('Evapotranspiration (mm)','FontSize',12); grid on;
legend([num2str(year(date_start_runoff_p1)) ' - ' num2str(year(date_end_runoff_p1))],[num2str(year(date_start_runoff_p2)) ' - ' num2str(year(date_end_runoff_p2))])

%% Seasonal air temperature trends

Hydro_year_label = strcat(string(num2str(Years_no_seas(1:end-1),'%02d ')), '/', ...
    string(num2str((1+Years_no_seas(1:end-1))-100*floor(Years_no_seas(1:end-1)/100))));

middle_year = year(date_start_runoff_p2);
ind_change_Ta = find(Years_no_seas==middle_year);

fi2 = figure('Renderer', 'painters', 'Position',[50 50 1.0474e+03 682.6667]);
tiledlayout(2,2,'TileSpacing','compact');
nexttile
plot(Years_no_seas(1:end-1),nanmean(Ta_yearly_ym(:,1:3),2),'-sqk','MarkerFaceColor','k','MarkerSize',5); grid on; hold on;
plot([Years_no_seas(1) Years_no_seas(ind_change_Ta)], [nanmean(Ta_yearly_ym(1:ind_change_Ta,1:3),'all') nanmean(Ta_yearly_ym(1:ind_change_Ta,1:3),'all')],'--k')
plot([Years_no_seas(ind_change_Ta) Years_no_seas(end)], [nanmean(Ta_yearly_ym(ind_change_Ta:end,1:3),'all') nanmean(Ta_yearly_ym(ind_change_Ta:end,1:3),'all')],'--k')
xline(Years_no_seas(ind_change_Ta),'r','LineWidth',0.6); %ylim([0 1800]);
xlim([Years_no_seas(1) Years_no_seas(end)]); ylabel('Air temperature [°C]','FontSize',12); xticks(Years_no_f2(1):2:Years_no_f2(end));
% ylim([-10.5 -5.5]); 
% text(2000, -6, 'Oct-Dec','FontWeight','bold','FontSize',12)
set(gca,'XTickLabels',Hydro_year_label(1:2:end))
title('Oct-Dec')

nexttile
plot(Years_no_seas(1:end-1),nanmean(Ta_yearly_ym(:,4:6),2),'-sqk','MarkerFaceColor','k','MarkerSize',5); grid on; hold on;
plot([Years_no_seas(1) Years_no_seas(ind_change_Ta)], [nanmean(Ta_yearly_ym(1:ind_change_Ta,4:6),'all') nanmean(Ta_yearly_ym(1:ind_change_Ta,4:6),'all')],'--k')
plot([Years_no_seas(ind_change_Ta) Years_no_seas(end)], [nanmean(Ta_yearly_ym(ind_change_Ta:end,4:6),'all') nanmean(Ta_yearly_ym(ind_change_Ta:end,4:6),'all')],'--k')
xline(Years_no_seas(ind_change_Ta),'r','LineWidth',0.6); %ylim([0 1800]);
xlim([Years_no_seas(1) Years_no_seas(end)]); ylabel('Air temperature [°C]','FontSize',12); xticks(Years_no_f2(1):2:Years_no_f2(end))
%ylim([-14.3 -9.3]); 
% text(2000, -9.8, 'Jan-Mar','FontWeight','bold','FontSize',12)
set(gca,'XTickLabels',Hydro_year_label(1:2:end))
title('Jan-Mar')

nexttile
plot(Years_no_seas(1:end-1),nanmean(Ta_yearly_ym(:,7:9),2),'-sqk','MarkerFaceColor','k','MarkerSize',5); grid on; hold on;
plot([Years_no_seas(1) Years_no_seas(ind_change_Ta)], [nanmean(Ta_yearly_ym(1:ind_change_Ta,7:9),'all') nanmean(Ta_yearly_ym(1:ind_change_Ta,7:9),'all')],'--k')
plot([Years_no_seas(ind_change_Ta) Years_no_seas(end)], [nanmean(Ta_yearly_ym(ind_change_Ta:end,7:9),'all') nanmean(Ta_yearly_ym(ind_change_Ta:end,7:9),'all')],'--k')
xline(Years_no_seas(ind_change_Ta),'r','LineWidth',0.6); %ylim([0 1800]);
xlim([Years_no_seas(1) Years_no_seas(end)]); ylabel('Air temperature [°C]','FontSize',12); xticks(Years_no_f2(1):2:Years_no_f2(end))
%ylim([-2.5 2.5]); 
% text(2000, 2, 'Apr-Jun','FontWeight','bold','FontSize',12)
title('Apr-Jun')
set(gca,'XTickLabels',Hydro_year_label(1:2:end))

nexttile
plot(Years_no_seas(1:end-1),nanmean(Ta_yearly_ym(:,10:12),2),'-sqk','MarkerFaceColor','k','MarkerSize',5); grid on; hold on;
plot([Years_no_seas(1) Years_no_seas(ind_change_Ta)], [nanmean(Ta_yearly_ym(1:ind_change_Ta,10:12),'all') nanmean(Ta_yearly_ym(1:ind_change_Ta,10:12),'all')],'--k')
plot([Years_no_seas(ind_change_Ta) Years_no_seas(end)], [nanmean(Ta_yearly_ym(ind_change_Ta:end,10:12),'all') nanmean(Ta_yearly_ym(ind_change_Ta:end,10:12),'all')],'--k')
xline(Years_no_seas(ind_change_Ta),'r','LineWidth',0.6); %ylim([0 1800]);
xlim([Years_no_seas(1) Years_no_seas(end)]); ylabel('Air temperature [°C]','FontSize',12); xticks(Years_no_f2(1):2:Years_no_f2(end))
%ylim([-2 2]); 
% text(2000, 7, 'Jul-Sep','FontWeight','bold','FontSize',12)
title('Jul-Sep')
set(gca,'XTickLabels',Hydro_year_label(1:2:end))
exportgraphics(fi2,[dir_fig '\Seasonal_catchment_temperature_timeseries.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Seasonal differences in snow fraction per elevation changes 


mid_year_id = find(Hydro_years==middle_year);

SR_el_ym = (Pr_ym_el - Pr_liq_ym_el)./Pr_ym_el;

fi3 =figure('Renderer', 'painters', 'Position',[237.6667 94.3333 599.3333 512.6667]);
tiledlayout(2,2,'TileSpacing','compact')
nexttile
p1 = plot(ELs, squeeze(nanmean(SR_el_ym(:,1:3,:),2)),'','HandleVisibility','off','LineWidth',0.4); grid on; hold on;
for ii = 1:size(squeeze(nanmean(SR_el_ym(:,1:3,:),2)),1)
    if Years_el(ii) < Years_el(mid_year_id)
        p1(ii,1).Color = [color_1970 0.25];
    else 
        p1(ii,1).Color = [color_2000 0.25];
    end
%     if ii == 1; p1(ii,1).HandleVisibility = "on"; end
end 
plot(ELs, nanmean(squeeze(nanmean(SR_el_ym(1:mid_year_id,1:3,:),2)),1),'LineWidth',1.3,'Color',color_1970,'HandleVisibility','off')
plot(ELs, nanmean(squeeze(nanmean(SR_el_ym(mid_year_id+1:end,1:3,:),2)),1),'LineWidth',1.3,'Color',color_2000,'HandleVisibility','off')
area(ELs, 3*(Hypso_el./nansum(Hypso_el)), 'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeColor','none')
area(ELs, 3*(Hypso_gla_el./nansum(Hypso_el)), 'FaceColor',[121, 208, 235]./255,'FaceAlpha',0.5,'EdgeColor','none')
ylim([0 1])
yticks([0:0.2:1]); 
view(90,-90); %set(gca,'XTickLabels',[]); 
lg1 = legend('Catchment hyspometry','Glacier hypsometry','Location','northwest','box','off');
lg1.FontSize = 8; lg1.NumColumns = 1; 
xlim([min(ELs) max(ELs)]); xlabel('Elevation [m a.s.l.]')
text(5400, 0.2, 'Oct-Dec','FontWeight','bold','FontSize',9)
lg1.Position = [0.19    0.75    0.2686    0.0634];

nexttile
p2 = plot(ELs, squeeze(nanmean(SR_el_ym(:,4:6,:),2)),'','HandleVisibility','off','LineWidth',0.4); grid on; hold on;
for ii = 1:size(squeeze(nanmean(SR_el_ym(:,4:6,:),2)),1)
    if Years_el(ii) < Years_el(mid_year_id)
        p2(ii,1).Color = [color_1970 0.25];
    else 
        p2(ii,1).Color = [color_2000 0.25];
    end
    if ii == 1; p2(ii,1).HandleVisibility = "on"; end
end 
plot(ELs, nanmean(squeeze(nanmean(SR_el_ym(1:mid_year_id,4:6,:),2)),1),'LineWidth',1.3,'Color',color_1970,'HandleVisibility','on')
plot(ELs, nanmean(squeeze(nanmean(SR_el_ym(mid_year_id+1:end,4:6,:),2)),1),'LineWidth',1.3,'Color',color_2000,'HandleVisibility','on')
area(ELs, 3*(Hypso_el./nansum(Hypso_el)), 'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off')
area(ELs, 3*(Hypso_gla_el./nansum(Hypso_el)), 'FaceColor',[121, 208, 235]./255,'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
ylim([0 1])
yticks([0:0.2:1]); 
view(90,-90); set(gca,'XTickLabels',[]); 
xlim([min(ELs) max(ELs)])
% lg2.Position = [0.14    0.8095    0.3200    0.0878];
text(5400, 0.2, 'Jan-Mar','FontWeight','bold','FontSize',9)

nexttile
p3 = plot(ELs, squeeze(nanmean(SR_el_ym(:,7:9,:),2)),'','HandleVisibility','off','LineWidth',0.4); grid on; hold on;
for ii = 1:size(squeeze(nanmean(SR_el_ym(:,7:9,:),2)),1)
    if Years_el(ii) < Years_el(mid_year_id)
        p3(ii,1).Color = [color_1970 0.25];
    else 
        p3(ii,1).Color = [color_2000 0.25];
    end
%     if ii == 1; p1(ii,1).HandleVisibility = "on"; end
end 
plot(ELs, nanmean(squeeze(nanmean(SR_el_ym(1:mid_year_id,7:9,:),2)),1),'LineWidth',1.3,'Color',color_1970,'HandleVisibility','off')
plot(ELs, nanmean(squeeze(nanmean(SR_el_ym(mid_year_id+1:end,7:9,:),2)),1),'LineWidth',1.3,'Color',color_2000,'HandleVisibility','off')
area(ELs, 3*(Hypso_el./nansum(Hypso_el)), 'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeColor','none')
area(ELs, 3*(Hypso_gla_el./nansum(Hypso_el)), 'FaceColor',[121, 208, 235]./255,'FaceAlpha',0.5,'EdgeColor','none')
ylabel('Snowfall fraction [-]','FontSize',11); ylim([0 1])
yticks(0:0.2:1); 
view(90,-90);  xlabel('Elevation [m a.s.l.]')
xlim([min(ELs) max(ELs)])
text(5400, 0.2, 'Apr-Jun','FontWeight','bold','FontSize',9)

nexttile
p4 = plot(ELs, squeeze(nanmean(SR_el_ym(:,10:12,:),2)),'','HandleVisibility','off','LineWidth',0.4); grid on; hold on;
for ii = 1:size(squeeze(nanmean(SR_el_ym(:,10:12,:),2)),1)
    if Years_el(ii) < Years_el(mid_year_id)
        p4(ii,1).Color = [color_1970 0.25];
    else 
        p4(ii,1).Color = [color_2000 0.25];
    end
     if ii == 1; p4(ii,1).HandleVisibility = "on"; end
end 
plot(ELs, nanmean(squeeze(nanmean(SR_el_ym(1:mid_year_id,10:12,:),2)),1),'LineWidth',1.3,'Color',color_1970,'HandleVisibility','on')
plot(ELs, nanmean(squeeze(nanmean(SR_el_ym(mid_year_id+1:end,10:12,:),2)),1),'LineWidth',1.3,'Color',color_2000,'HandleVisibility','on')
area(ELs, 3*(Hypso_el./nansum(Hypso_el)), 'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeColor','none','HandleVisibility','off')
area(ELs, 3*(Hypso_gla_el./nansum(Hypso_el)), 'FaceColor',[121, 208, 235]./255,'FaceAlpha',0.5,'EdgeColor','none','HandleVisibility','off')
ylabel('Snowfall fraction [-]','FontSize',11); ylim([0 1])
yticks(0:0.2:1); 
view(90,-90); set(gca,'XTickLabels',[]); 
xlim([min(ELs) max(ELs)])
text(5400, 0.2, 'Jul-Sep','FontWeight','bold','FontSize',9)
lg2 = legend('Hydrological years',[num2str(Hydro_years(1)) '-' num2str(Hydro_years(mid_year_id)) ' average'],...
    [num2str(Hydro_years(mid_year_id)) '-' num2str(Hydro_years(end)) ' average'],...
    'Catchment hyspometry','Glacier hypsometry','Location','southeast');
lg2.FontSize = 8.5; lg2.NumColumns = 1; lg2.Position; 
exportgraphics(fi3,[dir_fig '\Seasonal_catchment_snowfallFraction_altitudinal.png'],'Resolution',300,'BackgroundColor','none')

%% Compare runoff at the catchment-scale, with and without avalanches

if strcmp(glacier,'Kyzylsu')

sim_noav = 'ERA5Land_160824_1999_2023_PG00_Tmod0_2000mmSWEcap_noaval';
path_out_noav = ['C:\Users\jouberto\Desktop\T&C\Post-processing\Figures\Kyzylsu\' sim_noav];

Runoff_noav = load([path_out_noav '/Runoff/Runoff_Kyzylsu_monthly.mat'],'Runoff'); 
Runoff_noav = Runoff_noav.Runoff;

Pr_liq_y_el_noav = nanmean(squeeze(nansum(Runoff_noav.Pr_liq_ym_el,2)),1)'; 
Imelt_y_el_noav  = nanmean(squeeze(nansum(Runoff_noav.Imelt_ym_el,2)),1)'; 
Smelt_y_el_noav  = nanmean(squeeze(nansum(Runoff_noav.Smelt_ym_el,2)),1)';  
SSN_y_el_noav  = nanmean(squeeze(nansum(Runoff_noav.SSN_ym_el,2)),1)'; 
ET_y_el_noav  = nanmean(squeeze(nansum(Runoff_noav.ET_ym_el,2)),1)'; 

Pr_liq_y_el_av = nanmean(squeeze(nansum(Pr_liq_ym_el,2)),1)'; 
Imelt_y_el_av  = nanmean(squeeze(nansum(Imelt_ym_el,2)),1)'; 
Smelt_y_el_av  = nanmean(squeeze(nansum(Smelt_ym_el,2)),1)';  
SSN_y_el_av  = nanmean(squeeze(nansum(SSN_ym_el,2)),1)'; 
ET_y_el_av  = nanmean(squeeze(nansum(ET_ym_el,2)),1)'; 


fi2 = figure('Renderer', 'painters', 'Position',[-1.5577e+03 -225.6667 885.3667 784.6667]);
tiledlayout(5,4,'TileSpacing','tight')

nexttile(1,[3 1])
pl1 = plot(ELs, [Pr_liq_y_el_noav  Imelt_y_el_noav Smelt_y_el_noav].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
pl2 = plot(ELs, [-SSN_y_el_noav -ET_y_el_noav].*c_fact'.*10^-6,'LineWidth',2.5);
pl3 = plot(ELs, (Pr_liq_y_el_noav+Imelt_y_el_noav+Smelt_y_el_noav-SSN_y_el_noav-ET_y_el_noav).*c_fact'.*10^-6,'--k','LineWidth',1.1);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
yline(0,'Color',[0.2 0.2 0.2])
text(5500,2,{'Without','avalanches'},'FontWeight','bold','FontSize',10)
view(90,-90);grid on; ylim([-12 17])
ylabel('m^3 w.e. (*10^6)'); 
xlim([min(ELs) max(ELs)])
xlabel('Elevation [m a.s.l.]','FontSize',12); 
lg4 = legend('Rain ','Icemelt ','Snowmelt', 'Sublimation','Evapotranspiration','Net','location','northwest'); fontsize(lg4,9,"points");  
lg4.NumColumns = 6; lg4.Position = [0.095    0.936    0.6860    0.0246];

nexttile(2,[3 1])
pl1 = plot(ELs, [Pr_liq_y_el_av  Imelt_y_el_av Smelt_y_el_av].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
pl2 = plot(ELs, [-SSN_y_el_av -ET_y_el_av].*c_fact'.*10^-6,'LineWidth',2.5);
pl3 = plot(ELs, (Pr_liq_y_el_av+Imelt_y_el_av+Smelt_y_el_av-SSN_y_el_av-ET_y_el_av).*c_fact'.*10^-6,'--k','LineWidth',1.1);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
yline(0,'Color',[0.2 0.2 0.2])
text(5500,2,{'With', 'avalanches'},'FontWeight','bold','FontSize',10)
view(90,-90);grid on; ylim([-12 17]); xlim([min(ELs) max(ELs)])
%xlabel('Elevation [m a.s.l.]','FontSize',12); 
ylabel('m^3 w.e. (*10^6)'); 
set(gca,'XAxisLocation','top'); set(gca,'XTickLabel',[]);

nexttile(3,[3 1])
pl1 = plot(ELs, [Pr_liq_y_el_av-Pr_liq_y_el_noav  Imelt_y_el_av-Imelt_y_el_noav Smelt_y_el_av-Smelt_y_el_noav].*c_fact'.*10^-6,'LineWidth',2.5); hold on; 
pl2 = plot(ELs, [-SSN_y_el_av+SSN_y_el_noav -ET_y_el_av+ET_y_el_noav].*c_fact'.*10^-6,'LineWidth',2.5);
pl3 = plot(ELs, (Pr_liq_y_el_av+Imelt_y_el_av+Smelt_y_el_av-SSN_y_el_av-ET_y_el_av - (Pr_liq_y_el_noav+Imelt_y_el_noav+Smelt_y_el_noav-SSN_y_el_noav-ET_y_el_noav)).*c_fact'.*10^-6,'--k','LineWidth',1.1);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 256]./256; pl2(2).LineWidth = 2;
yline(0,'Color',[0.2 0.2 0.2])
text(5500,1,'Differences','FontWeight','bold','FontSize',10)
view(90,-90);grid on; ylim([-10 10]); xlim([min(ELs) max(ELs)])
ylabel('m^3 w.e. (*10^6)'); set(gca,'XTickLabel',[]);
set(gca,'XAxisLocation','top'); 

nexttile(4, [3 1])
b1 = area(ELs,lc_el,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.8,'LineStyle','none');
for i = 1:length(b1)
    set(b1(i), 'FaceColor', colors(i,:)) 
end
ylabel('Area [km^2]'); grid on;
xlim([min(ELs) max(ELs)]); view(-90,90); ylim([0 20]); set(gca,'XAxisLocation','top')
lg3 = legend(lc_labels,'FontSize',8.5,'Box','off','Location','NorthEast');
lg3.Position = [0.728    0.780    0.1502    0.1468];

xlabel('Elevation [m a.s.l.]','FontSize',12); 

nexttile(13,[2 4])
pl1 = plot(1:12, [nanmean(Pr_liq_ym,1)'  nanmean(Imelt_ym,1)' nanmean(Smelt_ym,1)'],'LineWidth',2); hold on; 
pl2 = plot(1:12, [-nanmean(SSN_ym,1)' -nanmean(ET_ym,1)'],'LineWidth',2);
pl1(2).Color = [0.65 0.65 0.65 1]; pl1(3).Color = [0.45 0.85 0.9 1]; pl1(1).Color = [0 0.6 1 1];
pl2(1).Color = [240 100 10 180]./256; pl2(2).Color =[0.4660 0.6740 0.1880];
plot(1:12, nanmean(nanmean(Pr_liq_ym,1)' + nanmean(Imelt_ym,1)'+ nanmean(Smelt_ym,1)' ...
    -nanmean(SSN_ym,1)' -nanmean(ET_ym,1)',2),'-k','LineWidth',1.2)
plot([0 1],[0 1],'--','Color',[0.8 0.8 0.8],'LineWidth',2)
pl1_noav = plot(1:12, [nanmean(Runoff_noav.Pr_liq_ym,1)'  nanmean(Runoff_noav.Imelt_ym,1)' nanmean(Runoff_noav.Smelt_ym,1)'],'LineWidth',2,'LineStyle','--'); 
pl2_noav = plot(1:12, [-nanmean(Runoff_noav.SSN_ym,1)' -nanmean(Runoff_noav.ET_ym,1)'],'LineWidth',2,'LineStyle','--');
pl1_noav(2).Color = [0.65 0.65 0.65 1]; pl1_noav(3).Color = [0.45 0.85 0.9 1]; pl1_noav(1).Color = [0 0.6 1 1];
pl2_noav(1).Color = [240 100 10 180]./256; pl2_noav(2).Color =[0.4660 0.6740 0.1880];
plot(1:12, nanmean(nanmean(Runoff_noav.Pr_liq_ym,1)' + nanmean(Runoff_noav.Imelt_ym,1)'+ nanmean(Runoff_noav.Smelt_ym,1)' ...
    -nanmean(Runoff_noav.SSN_ym,1)' -nanmean(Runoff_noav.ET_ym,1)',2),'--k','LineWidth',1.2)
xticks([1:1:12]); xlim([1 12]); ylim([-100 300]); 
set(gca,'XTickLabels',hydro_month_labels(1:1:12)); yline(0,'k','HandleVisibility','off','LineWidth',0.8)
ylabel('Runoff contribution (mm)','FontSize',12); grid on;
set(gca,'YAxisLocation','right')
lg2=legend('Rain','Icemelt','Snowmelt','Sublimation','Evapotranspiration','Net','Without avalanches','location','northwest'); fontsize(lg2,9,"points");  
exportgraphics(fi2,[dir_fig_runoff '\Elevation_fluxes_elev_2periods_lines_hypso_aval_comp.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

% Number for the paper

Total_runoff_av = sum(nanmean(nanmean(Pr_liq_ym,1)' + nanmean(Imelt_ym,1)'+ nanmean(Smelt_ym,1)' ...
    -nanmean(SSN_ym,1)' -nanmean(ET_ym,1)',2));

Total_runoff_noav = sum(nanmean(nanmean(Runoff_noav.Pr_liq_ym,1)' + nanmean(Runoff_noav.Imelt_ym,1)'+ nanmean(Runoff_noav.Smelt_ym,1)' ...
    -nanmean(Runoff_noav.SSN_ym,1)' -nanmean(Runoff_noav.ET_ym,1)',2));

Total_snowmelt_av = sum(nanmean(Smelt_ym,1));
Total_snowmelt_noav = sum(nanmean(Runoff_noav.Smelt_ym,1));

Total_icemelt_av = sum(nanmean(Imelt_ym,1));
Total_icemelt_noav = sum(nanmean(Runoff_noav.Imelt_ym,1));

100*(Total_runoff_noav-Total_runoff_av)./Total_runoff_noav
100*(Total_snowmelt_noav-Total_snowmelt_av)./Total_snowmelt_noav
100*(Total_icemelt_noav-Total_icemelt_av)./Total_icemelt_noav

end 
