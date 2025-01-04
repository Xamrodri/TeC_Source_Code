Hydro_years = year(date_m(find(month(date_m) == 10 & day(date_m) == 1,1))):1:year(date_m(find(month(date_m) == 9 & day(date_m) == 1,1,'last')));

[Ta_hy, SMB_hy, Pr_hy, Psno_hy, Psno_spring, Psno_winter, Psno_mon, Pliq_hy, ...
    Smg_hy, Sms_hy, Esn_hy, AVA_hy] = deal(NaN(1,numel(Hydro_years)-1));


if length(Hydro_years) > 1

for yy = 1:length(Hydro_years)-1

    ind_start = find(year(date_m) == Hydro_years(yy) & month(date_m) == 10 & day(date_m) == 1,1);
    ind_end =   find(year(date_m) == Hydro_years(yy+1) & month(date_m) == 9 & day(date_m) == 1,1);
    ind_period = ind_start:ind_end;

    ind_octfeb = month(date_m) > 9 | month(date_m) < 3;
    ind_marmay = month(date_m) > 2 & month(date_m) < 6;
    ind_junsep = month(date_m) > 5 & month(date_m) < 10;

    if ~exist('AVA_map','var'); AVA_map = PSNOW_map.*0; end % For multipoints analysis
    if ~exist('EICE_map','var'); EICE_map = PSNOW_map.*0; end % For multipoints analysis

%     SMB_map = ICE_map(:,:,ind_end) + SWEm_map(:,:,ind_end) - ICE_map(:,:,ind_start) - SWEm_map(:,:,ind_start);
    SMB_map = nansum(PSNOW_map(:,:,ind_start:ind_end),3) + nansum(AVA_map(:,:,ind_start:ind_end),3) ...
        - nansum(SMSm_map(:,:,ind_start:ind_end),3) - nansum(ESN_map(:,:,ind_start:ind_end),3) ...
        - nansum(SMG_map(:,:,ind_start:ind_end),3) - nansum(EICE_map(:,:,ind_start:ind_end),3);

    SMB_map(GLA_ID ~= gla_id_smb) = NaN;

    SMB_hy(yy) = nanmean(SMB_map,'all');

    Ta_hy_map = nanmean(TA_map(:,:,ind_start:ind_end),3);
    Ta_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Ta_hy(yy) = nanmean(Ta_hy_map,'all');

    Pr_hy_map = nansum(PRECIP_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps
    Pr_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Pr_hy(yy) = nanmean(Pr_hy_map,'all');

    Psno_hy_map = nansum(PSNOW_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps
    Psno_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Psno_hy(yy) = nanmean(Psno_hy_map,'all');

    Psno_spring_map = nansum(PSNOW_map(:,:,ind_period(ind_marmay(ind_period))),3); % Annual precipitation sum maps
    Psno_spring_map(GLA_ID ~= gla_id_smb) = NaN; 
    Psno_spring(yy) = nanmean(Psno_spring_map,'all');

    Psno_winter_map = nansum(PSNOW_map(:,:,ind_period(ind_octfeb(ind_period))),3); % Annual precipitation sum maps
    Psno_winter_map(GLA_ID ~= gla_id_smb) = NaN; 
    Psno_winter(yy) = nanmean(Psno_winter_map,'all');   

    Psno_mon_map = nansum(PSNOW_map(:,:,ind_period(ind_junsep(ind_period))),3); % Annual precipitation sum maps
    Psno_mon_map(GLA_ID ~= gla_id_smb) = NaN; 
    Psno_mon(yy) = nanmean(Psno_mon_map,'all');

    Pliq_hy_map = nansum(PRAIN_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps
    Pliq_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Pliq_hy(yy) = nanmean(Pliq_hy_map,'all');

    Smg_hy_map = nansum(SMG_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps
    Smg_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Smg_hy(yy) = nanmean(Smg_hy_map,'all');

    Sms_hy_map = nansum(SMSm_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps
    Sms_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Sms_hy(yy) = nanmean(Sms_hy_map,'all');

    Esn_hy_map = nansum(ESN_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps
    Esn_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Esn_hy(yy) = nanmean(Esn_hy_map,'all');    

    if exist('AVA_map','var')
    AVA_hy_map = nansum(AVA_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps
    AVA_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    AVA_hy(yy) = nanmean(AVA_hy_map,'all');  
    else
    AVA_hy(yy) = 0;
    end 

end 
end 

% if length(Hydro_years) <= 2

    if ~exist('AVA_map','var'); AVA_map = PSNOW_map.*0; end % For multipoints analysis
    if ~exist('EICE_map','var'); EICE_map = PSNOW_map.*0; end % For multipoints analysis

%     SMB_map = ICE_map(:,:,ind_end) + SWEm_map(:,:,ind_end) - ICE_map(:,:,ind_start) - SWEm_map(:,:,ind_start);
    SMB_map = nansum(PSNOW_map,3) + nansum(AVA_map,3) ...
        - nansum(SMSm_map,3) - nansum(ESN_map,3) ...
        - nansum(SMG_map,3) - nansum(EICE_map,3);

    SMB_map(GLA_ID ~= gla_id_smb) = NaN;

    SMB_hy = nanmean(SMB_map,'all');

    Ta_hy_map = nanmean(TA_map,3);
    Ta_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Ta_hy = nanmean(Ta_hy_map,'all');

    Pr_hy_map = nansum(PRECIP_map,3); % Annual precipitation sum maps
    Pr_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Pr_hy = nanmean(Pr_hy_map,'all');

    Psno_hy_map = nansum(PSNOW_map,3); % Annual precipitation sum maps
    Psno_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Psno_hy = nanmean(Psno_hy_map,'all');

    Pliq_hy_map = nansum(PRAIN_map,3); % Annual precipitation sum maps
    Pliq_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Pliq_hy = nanmean(Pliq_hy_map,'all');

    Smg_hy_map = nansum(SMG_map,3); % Annual precipitation sum maps
    Smg_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Smg_hy = nanmean(Smg_hy_map,'all');

    Sms_hy_map = nansum(SMSm_map,3); % Annual precipitation sum maps
    Sms_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Sms_hy = nanmean(Sms_hy_map,'all');

    Esn_hy_map = nansum(ESN_map,3); % Annual precipitation sum maps
    Esn_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    Esn_hy = nanmean(Esn_hy_map,'all');    

    if exist('AVA_map','var')
    AVA_hy_map = nansum(AVA_map,3); % Annual precipitation sum maps
    AVA_hy_map(GLA_ID ~= gla_id_smb) = NaN; 
    AVA_hy = nanmean(AVA_hy_map,'all');  
    else
    AVA_hy = 0;
    end 


% Find rows and columns containing only zeros
SMB_flipped = flipud(SMB_map);
zero_rows = all(isnan(SMB_flipped), 2);
zero_cols = all(isnan(SMB_flipped), 1);

% Find the minimum and maximum indices of rows and columns containing ones
min_row = find(~zero_rows, 1, 'first');
max_row = find(~zero_rows, 1, 'last');
min_col = find(~zero_cols, 1, 'first');
max_col = find(~zero_cols, 1, 'last');

buffer = 6;
min_row = max(1, min_row - buffer);
max_row = min(size(SMB_flipped, 1), max_row + buffer);
min_col = max(1, min_col - buffer);
max_col = min(size(SMB_flipped, 2), max_col + buffer);

% Crop the array based on minimum and maximum indices
cropped_SMB = SMB_flipped(min_row:max_row, min_col:max_col);

cropped_x = demLons(min_row:max_row, min_col:max_col);
cropped_y = demLats(min_row:max_row, min_col:max_col);

%% For the whole period
fi3 =figure('Renderer', 'painters', 'Position', [209 251.6667 906.6667 416]);
tiledlayout(1,2,'TileSpacing','compact')
nexttile
imagesc(cropped_x(1,:), cropped_y(:,1),cropped_SMB.*0.001./length(Hydro_years),'AlphaData',~isnan(cropped_SMB)); hold on;
cb = colorbar;
ylabel(cb,'Glacier surface mass balance [m w.e./yr]','FontSize',11)
set(gca,'YDir','normal'); 
%set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
title([datestr(date_m(end),'dd-mmm-yyyy') ' - ' datestr(date_m(1),'dd-mmm-yyyy')]); 
clim([-5 5])
xlabel({['Glacier mean SMB: ' num2str(round(SMB_hy./(1000.*length(Hydro_years)),2)) ' m w.e./yr'],['RGI ID: ' num2str(gla_id_smb)]},'FontSize',10)
colormap(flipud(redblue))
set(gca,'Color',[0.8 0.8 0.8])
nexttile
bar(1:6, [Psno_hy -Smg_hy -Sms_hy -Esn_hy AVA_hy SMB_hy]./(1000.*length(Hydro_years))); grid on;
ylabel('Mass balance components [m w.e./yr]','FontSize',11)
set(gca,'YAxisLocation','right')
set(gca,'XTickLabels', {'Snowfall','Icemelt','Snowmelt','Evap/Subli','Avalanches','net SMB'})
exportgraphics(fi3,[dir_fig '\GMB\' glacier '_SMB_component_summary_' num2str(gla_id_smb) '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)