% Attributes the glacier mass changes between two periods to changes in
% precipitation, precipitation phase, snowmelt and icemelt amount
% This is done at the hydrological year level.

period_start = find(month(date_attr) == 10 & day(date_attr) == 1,1);
period_end = find(month(date_attr) == 9,1,'last');

period_middle = 2012;

seasons = [10 11 12 1; 2 3 4 5; 6 7 8 9]; 
seasons_HY = [1 2 3 4; 5 6 7 8; 9 10 11 12]; 
seasons_labels = {'Oct-Jan','Feb-May','Jun-Sep'};

Years_no_seas = unique(year(date_attr(period_start:period_end)));

Hydro_year_label = strcat(string(num2str(Years_no_seas(1:end-1)-100*floor(Years_no_seas(1:end-1)/100),'%02d ')), '/', ...
    string(num2str(Years_no_seas(2:end)-100*floor(Years_no_seas(2:end)/100),'%02d ')));

Month_seas = month(date_attr(period_start:period_end));

[Ta_yearly_ym, Pr_ym, Pr_sno_ym, SMB_ym, Imelt_ym, Smelt_ym, ESN_ym] = deal(NaN(numel(Years_no_seas)-1, 12));

hydro_month = [10,11,12,1,2,3,4,5,6,7,8,9];
hydro_month_labels = {'Oct','Nov','Dec','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep'};

% Compute monthly variable for all hydrological years

for yy = 1:length(Years_no_seas)-1
    for mm = 1:12

      if hydro_month(mm) > 8
        ind_start = find(year(date_m) == Years_no_seas(yy) & month(date_m) == hydro_month(mm),1,'first');
        ind_end =   find(year(date_m) == Years_no_seas(yy) & month(date_m) == hydro_month(mm),1,'last');
      else
        ind_start = find(year(date_m) == (Years_no_seas(yy)+1) & month(date_m) == hydro_month(mm),1,'first');
        ind_end =   find(year(date_m) == (Years_no_seas(yy)+1) & month(date_m) == hydro_month(mm),1,'last');
      end 

%     Ta_yearly_ym(yy,mm) = nanmean(TA_map(:,:,ind_start:ind_end),'all');

    Pr_y_map = nansum(PRECIP_map(:,:,ind_start:ind_end),3); Pr_y_map(GLA_ID ~= gla_id_smb) = NaN;
    Pr_sno_y_map = nansum(PSNOW_map(:,:,ind_start:ind_end),3); Pr_sno_y_map(GLA_ID ~= gla_id_smb) = NaN;

    Pr_ym(yy,mm) = nanmean(Pr_y_map,'all'); 
    Pr_sno_ym(yy,mm) = nanmean(Pr_sno_y_map,'all'); 

    SMB_mapi = nansum(PSNOW_map(:,:,ind_start:ind_end),3) + nansum(AVA_map(:,:,ind_start:ind_end),3).*1000 ...
        - nansum(SMSm_map(:,:,ind_start:ind_end),3) - nansum(ESN_map(:,:,ind_start:ind_end),3) ...
        - nansum(SMG_map(:,:,ind_start:ind_end),3) - nansum(EICE_map(:,:,ind_start:ind_end),3);
    SMB_mapi(GLA_ID ~= gla_id_smb) = NaN; % GMB at main glacier
    SMB_ym(yy,mm) = nanmean(SMB_mapi,'all');

    Imelt_y_map = nansum(SMG_map(:,:,ind_start:ind_end),3); Imelt_y_map(GLA_ID ~= gla_id_smb) = NaN;
    Imelt_ym(yy,mm) = nanmean(Imelt_y_map,'all');

    Smelt_y_map = nansum(SMSm_map(:,:,ind_start:ind_end),3); Smelt_y_map(GLA_ID ~= gla_id_smb) = NaN;
    Smelt_ym(yy,mm) = nanmean(Smelt_y_map,'all');

    ESN_y_map = nansum(ESN_map(:,:,ind_start:ind_end),3); ESN_y_map(GLA_ID ~= gla_id_smb) = NaN;
    ESN_ym(yy,mm) = nanmean(ESN_y_map,'all');

    SSN_y_map = nansum(SSN_map(:,:,ind_start:ind_end),3); SSN_y_map(GLA_ID ~= gla_id_smb) = NaN;
    SSN_ym(yy,mm) = nanmean(SSN_y_map,'all');   
    end
end 

% Compute aggregated variables per season for the two periods of interest

ind_p1 = find(Years_no_seas(1:end-1) < period_middle);
ind_p2 = find(Years_no_seas(1:end-1) >= period_middle);

for ss = 1:size(seasons_HY,1)
    
    psnow_p1_s(ss,:) = nanmean(nansum(Pr_sno_ym(ind_p1,seasons_HY(ss,:)),2));
    psnow_p2_s(ss,:) = nanmean(nansum(Pr_sno_ym(ind_p2,seasons_HY(ss,:)),2));

    ptot_p1_s(ss,:) = nanmean(nansum(Pr_ym(ind_p1,seasons_HY(ss,:)),2));
    ptot_p2_s(ss,:) = nanmean(nansum(Pr_ym(ind_p2,seasons_HY(ss,:)),2));   

    rsnow_p1(ss,:) = psnow_p1_s(ss,:)./ptot_p1_s(ss,:);
    rsnow_p2(ss,:) = psnow_p2_s(ss,:)./ptot_p2_s(ss,:);

    Imelt_p1(ss,:) = nanmean(nansum(Imelt_ym(ind_p1,seasons_HY(ss,:)),2));
    Imelt_p2(ss,:) = nanmean(nansum(Imelt_ym(ind_p2,seasons_HY(ss,:)),2));

    Smelt_p1(ss,:) = nanmean(nansum(Smelt_ym(ind_p1,seasons_HY(ss,:)),2));
    Smelt_p2(ss,:) = nanmean(nansum(Smelt_ym(ind_p2,seasons_HY(ss,:)),2));

    ESN_p1(ss,:) = nanmean(nansum(ESN_ym(ind_p1,seasons_HY(ss,:)),2));
    ESN_p2(ss,:) = nanmean(nansum(ESN_ym(ind_p2,seasons_HY(ss,:)),2));

end 

daccum_precip_s = (ptot_p2_s - ptot_p1_s).*rsnow_p1;
daccum_phase_s = ptot_p2_s.*(rsnow_p2- rsnow_p1);
diff_icemelt_s = Imelt_p2 - Imelt_p1;
diff_snowmelt_s = Smelt_p2 - Smelt_p1;
diff_esn = ESN_p2 - ESN_p1;

diff_all_s = [daccum_precip_s daccum_phase_s -diff_icemelt_s -diff_snowmelt_s -diff_esn];

pie_label = {{'Precipitation';'amount change'}, {'Precipitation';'phase change'}, 'Icemelt \uparrow','Snowmelt \uparrow', 'Sublimation \uparrow'};
pie_legend = {'Precip_{amount}', 'Precip_{phase}', 'Icemelt', 'Snowmelt','Sublimation'};

d_smb_s = nansum(diff_all_s,2);

%% Pie charts

fi6 = figure('Renderer', 'painters', 'Position',[166 120.3333 631 551.6667]);
tiledlayout(2,2,'TileSpacing','compact')    

for ss = 1:size(seasons_HY,1)
    nexttile
    x_neg = diff_all_s(ss,:); x_neg(x_neg>0)=0.001; 
    h1 = pie(abs(x_neg));
    h1(1,1).FaceColor = [0.8039    0.6314    1.0000];
    h1(1,3).FaceColor = [0 0.6 1];
    h1(1,5).FaceColor = [0.85 0.85 0.85];
    h1(1,7).FaceColor = [0 1 1];
    h1(1,9).FaceColor = [0.65 0.95 1]; 
    delete(findobj(h1,'Type','text'))%option 2
    title(seasons_labels{ss})
    text(-0.45,-1.2,[num2str(nansum(x_neg),3) ' mm w.e.'],'FontSize',9)
end
nexttile
x_neg = nansum(diff_all_s,1); x_neg(x_neg>0)=0.001; 
h1 = pie2(abs(x_neg),pie_label);
set(h1(2:2:end),'FontSize',9)
    h1(1,1).FaceColor = [0.8039    0.6314    1.0000];
    h1(1,3).FaceColor = [0 0.6 1];
    h1(1,5).FaceColor = [0.85 0.85 0.85];
    h1(1,7).FaceColor = [0 1 1];
    h1(1,9).FaceColor = [0.65 0.95 1]; 
    xlabel([num2str(nansum(x_neg),3) ' mm w.e.'])
    text(-0.45,-1.6,[num2str(nansum(x_neg),3) ' mm w.e.'],'FontSize',9)
exportgraphics(fi6,[dir_fig '\PieChart_GMBchange_attribution.png'],'Resolution',300,'BackgroundColor','none')



