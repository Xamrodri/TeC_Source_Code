
if ~exist('date_seas','var')
    date_seas = date_m;
end 

Year_no = Years_no;
date_modis = modis_fsc_40.dateTime;
Years_modis = unique(year(date_modis));

fsc_modis_tt = timetable(date_modis, fsc_modis_filt_40,'VariableNames',{'fsc'});
fsc_modis_mm = retime(fsc_modis_tt,"monthly",@nanmean);

% MODIS

[fsc_MODIS_ym, fsc_MODIS_y_std, se_MODIS_ym, se_MODIS_ym_std] = deal(NaN(length(Years_no),12));
[fsc_MODIS_y, fsc_MODIS_win, fsc_MODIS_sum, fsc_MODIS_spr, fsc_MODIS_aut, fsc_MODIS_y_std] = deal(NaN(length(Years_no),1));

for yy = 1:length(Years_no)
    for m = 1:12

    ind_ym = (year(date_modis) == Years_no(yy)) & (month(date_modis) == m);

    fsc_MODIS_ym(yy,m) = nanmean(fsc_modis_filt_40(ind_ym));
    fsc_MODIS_y_std(yy,m) = nanstd(fsc_modis_filt_40(ind_ym));

    se_MODIS_ym(yy,m) = nanmean(se_modis_filt_40(ind_ym));
    se_MODIS_ym_std(yy,m) = nanstd(se_modis_filt_40(ind_ym));    

    end
    fsc_MODIS_win(yy) = nanmean(fsc_modis_filt_40(year(date_modis) == Years_no(yy) & ismember(month(date_modis), [12 1 2])));
    fsc_MODIS_spr(yy) = nanmean(fsc_modis_filt_40(year(date_modis) == Years_no(yy) & ismember(month(date_modis), [3 4 5])));
    fsc_MODIS_aut(yy) = nanmean(fsc_modis_filt_40(year(date_modis) == Years_no(yy) & ismember(month(date_modis), [9 10 11])));
    fsc_MODIS_sum(yy) = nanmean(fsc_modis_filt_40(year(date_modis) == Years_no(yy) & ismember(month(date_modis), [6 7 8])));
    fsc_MODIS_y(yy) = nanmean(fsc_modis_filt_40(year(date_modis) == Years_no(yy)));
    fsc_MODIS_y_std(yy) = nanstd(fsc_modis_filt_40(year(date_modis) == Years_no(yy))); 
end 

% Landsat/Sentinel-2 

Year_no_LS = unique(year(LS_sno.dateTime));

[fsc_LS_ym, fsc_LS_y_std, se_LS_ym, se_LS_ym_std] = deal(NaN(length(Year_no_LS),12));
[fsc_LS_y, fsc_LS_y_std] = deal(NaN(length(Year_no_LS),1));

for yy = 1:length(Year_no_LS)
    for m = 1:12

    ind_ym = (year(LS_sno.dateTime) == Year_no_LS(yy)) & (month(LS_sno.dateTime) == m);

    fsc_LS_ym(yy,m) = nanmean(LS_sno.fsnow_hybrid(ind_ym));
    fsc_LS_y_std(yy,m) = nanstd(LS_sno.fsnow_hybrid(ind_ym));

    se_LS_ym(yy,m) = nanmean(LS_sno.se_hybrid(ind_ym));
    se_LS_ym_std(yy,m) = nanstd(LS_sno.se_hybrid(ind_ym));    

    end

    if nansum(isnan(fsc_LS_ym(yy,:))) < 4
       fsc_LS_y(yy) = nanmean(LS_sno.fsnow_hybrid(year(LS_sno.dateTime) == Year_no_LS(yy)));
       fsc_LS_y_std(yy) = nanstd(LS_sno.fsnow_hybrid(year(LS_sno.dateTime) == Year_no_LS(yy))); 
    end 
end 

% T&C

[fsc_tc_ym, fsc_tc_y_std, se_tc_ym, se_tc_ym_std] = deal(NaN(length(Years_no),12));
[fsc_tc_y, fsc_tc_win, fsc_tc_sum, fsc_tc_spr, fsc_tc_aut, fsc_tc_y_std] = deal(NaN(length(Years_no),1));

for yy = 1:length(Years_no)
    for m = 1:12

    ind_ym = (year(Date_d) == Years_no(yy)) & (month(Date_d) == m);

    fsc_tc_ym(yy,m) = nanmean(scas(ind_ym));
    fsc_tc_y_std(yy,m) = nanstd(scas(ind_ym));

    se_tc_ym(yy,m) = nanmean(se(ind_ym));
    se_tc_ym_std(yy,m) = nanstd(se(ind_ym));    

    end
   fsc_tc_spr(yy) = nanmean(scas(year(Date_d) == Years_no(yy) & ismember(month(Date_d), [3 4 5])) );
   fsc_tc_aut(yy) = nanmean(scas(year(Date_d) == Years_no(yy) & ismember(month(Date_d), 9:11)) );
   fsc_tc_win(yy) = nanmean(scas(year(Date_d) == Years_no(yy) & ismember(month(Date_d), [12 1 2])) );
   fsc_tc_sum(yy) = nanmean(scas(year(Date_d) == Years_no(yy) & ismember(month(Date_d), 6:8)) );
   fsc_tc_y(yy) = nanmean(scas(year(Date_d) == Years_no(yy)));
   fsc_LS_y_std(yy) = nanstd(scas(year(Date_d) == Years_no(yy))); 
end 

%% Compute mean snow depth per elevation band per hydrological years

period_start = find(month(date_seas) == 9 & day(date_seas) == 1,1);
period_end = find(month(date_seas) == 8,1,'last');

Years_no_seas = unique(year(date_seas(period_start:period_end)));


Hydro_year_label = strcat(string(num2str(Years_no_seas(1:end-1)-100*floor(Years_no_seas(1:end-1)/100),'%02d ')), '/', ...
    string(num2str(Years_no_seas(2:end)-100*floor(Years_no_seas(2:end)/100),'%02d ')));

Month_seas = month(date_seas(period_start:period_end));


dEL=100; % width of elevation bins
ELs = nanmin(DTM,[],'all'):dEL:nanmax(DTM,[],'all');

 [SWE_y_el] = deal(NaN(numel(Years_no_seas)-1, 12, numel(ELs)));

  [SWE_y_el_spr, SWE_y_el_aut, SWE_y_el_win, SWE_y_el_sum] = deal(NaN(numel(Years_no_seas)-1, numel(ELs)));

hydro_month = [9,10,11,12,1,2,3,4,5,6,7,8];
hydro_month_labels = {'Sep','Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug'};

for yy = 1:length(Years_no_seas)-1
    for mm = 1:12

      if hydro_month(mm) > 8
        ind_start = find(year(date_m) == Years_no_seas(yy) & month(date_m) == hydro_month(mm),1,'first');
        ind_end =   find(year(date_m) == Years_no_seas(yy) & month(date_m) == hydro_month(mm),1,'last');

        ind_start_d = find(year(Date_d) == Years_no_seas(yy) & month(Date_d) == hydro_month(mm),1,'first');
        ind_end_d =   find(year(Date_d) == Years_no_seas(yy) & month(Date_d) == hydro_month(mm),1,'last');

        date_ym(yy,mm) = datetime(Years_no_seas(yy), hydro_month(mm),15);

      else
        ind_start = find(year(date_m) == (Years_no_seas(yy)+1) & month(date_m) == hydro_month(mm),1,'first');
        ind_end =   find(year(date_m) == (Years_no_seas(yy)+1) & month(date_m) == hydro_month(mm),1,'last');

        ind_start_d = find(year(Date_d) == (Years_no_seas(yy)+1) & month(Date_d) == hydro_month(mm),1,'first');
        ind_end_d =   find(year(Date_d) == (Years_no_seas(yy)+1) & month(Date_d) == hydro_month(mm),1,'last');

        date_ym(yy,mm) = datetime(Years_no_seas(yy)+1, hydro_month(mm),15);
      end 

       SWE_y_map = nanmean(SWEm_map(:,:,ind_start:ind_end),3); SWE_y_map(MASK ~=1) = NaN;

        for iel = 1:numel(ELs)
              cur=(DTM<(ELs(iel)+dEL/2))&(DTM>=(ELs(iel)-dEL/2)); %current section of DEM
              SWE_y_el(yy,mm,iel) = nanmean(SWE_y_map(cur));
        end 

    end

    SWE_y_el_spr(yy,:) = nanmean(SWE_y_el(yy,[7 8 9],:),2);
    SWE_y_el_aut(yy,:) = nanmean(SWE_y_el(yy,[1 2 3],:),2);
    SWE_y_el_win(yy,:) = nanmean(SWE_y_el(yy,[4 5 6],:),2);
    SWE_y_el_sum(yy,:) = nanmean(SWE_y_el(yy,[10 11 12],:),2);
end

SWE_y_el_spr_mean = nanmean(SWE_y_el_spr,1);
SWE_y_el_aut_mean = nanmean(SWE_y_el_aut,1);
SWE_y_el_win_mean = nanmean(SWE_y_el_win,1);
SWE_y_el_sum_mean = nanmean(SWE_y_el_sum,1); 
SWE_y_el_seas_mean = [SWE_y_el_aut_mean' SWE_y_el_win_mean' SWE_y_el_spr_mean' SWE_y_el_sum_mean'];

SWE_y_el_spr_std = nanstd(SWE_y_el_spr,1);
SWE_y_el_aut_std = nanstd(SWE_y_el_aut,1);
SWE_y_el_win_std = nanstd(SWE_y_el_win,1);
SWE_y_el_sum_std = nanstd(SWE_y_el_sum,1);
SWE_y_el_seas_std = [SWE_y_el_aut_std' SWE_y_el_win_std' SWE_y_el_spr_std' SWE_y_el_sum_std'];

season_labels = {'Sep - Nov','Dec - Feb','Mar - May','Jun - Aug'};

fi2 = figure('Renderer', 'painters', 'Position',[223.6667 77 607.3333 516.6666]) ;
tiledlayout(2,2,"TileSpacing","tight")
for ii = 1:4
nexttile
errorbar(ELs./1000,SWE_y_el_seas_mean(:,ii),SWE_y_el_seas_std(:,ii),'r'); hold on; grid on;
plot(ELs./1000,SWE_y_el_seas_mean(:,ii),'k','LineWidth',1.2)
if ii == 3; legend('1 \sigma','Mean','Location','northwest'); end
view(90,-90); hold on; grid on; ylim([0 700])
xlabel('Elevation [km]'); ylabel('Mean SWE [mm w.e.]'); yticks(0:100:700)
text(4.7,450,season_labels{ii},'FontWeight','bold'); xlim([2 5])
end 
exportgraphics(fi2,[dir_fig '\Snow_cover\SWE_seasonal_elevation_profile.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Final figure

fi3 =figure('Renderer', 'painters', 'Position',  [209 146.3333 522.6667 408.6667]);
tiledlayout(2,2,'TileSpacing','tight')
nexttile
plot(Years_no, fsc_MODIS_aut,'-k'); grid on; hold on;
plot(Years_no,fsc_tc_aut,'--k');
set(gca,'XAxisLocation','top')
ylabel('Snow cover fraction [-]','FontSize',10)
lg1 = legend('MODIS','Model','location','south');
lg1.NumColumns = 2;
ylim([0 1]); xlim([2000 2023])
text(2005,0.8,'Sep - Nov','FontWeight','bold','FontSize',9)
nexttile
plot(Years_no, fsc_MODIS_win,'-k','HandleVisibility','off');  grid on; hold on;
plot(Years_no,fsc_tc_win,'--k','HandleVisibility','off'); 
plot(Years_no, fsc_MODIS_y,'-r','HandleVisibility','on');  grid on; hold on;
plot(Years_no,fsc_tc_y,'--r','HandleVisibility','off'); 
legend('Annual','Location','south')

set(gca,'XAxisLocation','top')
ylim([0 1]); xlim([2000 2023]); set(gca,'YTickLabels',[])
text(2005,0.8,'Dec - Feb','FontWeight','bold','FontSize',9)
nexttile
plot(Years_no, fsc_MODIS_spr,'-k');  grid on; hold on;
plot(Years_no,fsc_tc_spr,'--k'); 
ylim([0 1]); xlim([2000 2023])
ylabel('Snow cover fraction [-]','FontSize',10)
text(2005,0.6,'Mar - May','FontWeight','bold','FontSize',9)
nexttile
plot(Years_no, fsc_MODIS_sum,'-k'); grid on; hold on;
plot(Years_no,fsc_tc_sum,'--k'); set(gca,'YTickLabels',[])
ylim([0 1]); xlim([2000 2023])
text(2005,0.6,'Jun - Aug','FontWeight','bold','FontSize',9)
exportgraphics(fi3,[dir_fig '\Snow_cover\' glacier '_snow_cover_seasonal_trends_MODIS_comp.png'],'Resolution',300,'BackgroundColor','none')
