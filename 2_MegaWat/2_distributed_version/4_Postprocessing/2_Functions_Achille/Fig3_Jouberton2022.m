%% Figure 3 of Jouberton et al. 2022

% needs:
% - Monthly or daily maps of: air temperature, sold and total
% precipitation, avalanches, snowmelt, icemelt, evaporation (Ice and snow)
% - Corresponding time vector (date_m)

if ~exist('gla_fig3_id','var')
    gla_fig3_id = gla_id;
end 

if ~exist('gla_fig3_label','var')
    gla_fig3_label = glacier;
end 

if ~exist('date_f3','var')
    date_f3 = date_m;
end 

period_start = min(find(month(date_f3) == 10 & day(date_f3) == 1,1), find(month(date_f3) == 11 & day(date_f3) == 1,1));
period_end = find(month(date_f3) == 9,1,'last');

Years_f3 = unique(year(date_f3(period_start:period_end)));
Years_f3_all = Years_f3;
Years_f3 = Years_f3(1:end-1);

[PrecS_spring_y, PrecS_monsoon_y, PrecS_octfeb_y, PrecS_y, Prec_spring_y,...
  Prec_monsoon_y, Prec_octfeb_y, Prec_y, Smg_spring_y, Smg_monsoon_y, Smg_octfeb_y,...
  Smg_y, Sms_spring_y, Sms_monsoon_y, Sms_octfeb_y, Sms_y,...
  ESN_spring_y, ESN_monsoon_y, ESN_octfeb_y, ESN_y, Net_y, AVA_y, ...
  PrecS_octdec_y, PrecS_janmar_y, PrecS_aprjun_y, PrecS_julsep_y, ...
  Prec_octdec_y, Prec_janmar_y, Prec_aprjun_y, Prec_julsep_y] = deal(NaN(1,numel(Years_f3)));

month_m = date_m-calmonths(1);

for yy = 1:length(Years_f3)

    autumn_idx = year(month_m) == Years_f3(yy) & ismember(month(month_m), [10 11]);    
    spring_idx = year(month_m) == Years_f3(yy)+1 & ismember(month(month_m), [3 4 5]);   
    monsoon_idx = year(month_m) == Years_f3(yy)+1 & ismember(month(month_m), 6:9); 
%   summer_idx = year(month_m) == Years_f3(yy)+1 & ismember(month_m, [6 7 8]); 
    winter_idx = (year(month_m) == Years_f3(yy) & ismember(month(month_m), [12])) | ...
        (year(month_m) == Years_f3(yy)+1 & ismember(month(month_m), [1 2]));        

    Oct_Dec_idx = year(month_m) == Years_f3(yy) & ismember(month(month_m), [10 11 12]);    
    Jan_Mar_idx = year(month_m) == Years_f3(yy)+1 & ismember(month(month_m), [1 2 3]);
    Apr_Jun_idx = year(month_m) == Years_f3(yy)+1 & ismember(month(month_m), [4 5 6]);
    Jul_Sep_idx = year(month_m) == Years_f3(yy)+1 & ismember(month(month_m), [7 8 9]);
    
    % Snowfalls
    Psno_spring = nansum(PSNOW_map(:,:,spring_idx),3); Psno_spring(GLA_ID ~= gla_fig3_id) = NaN;
    PrecS_spring_y(yy) = nanmean(Psno_spring,'all');

    Psno_monsoon = nansum(PSNOW_map(:,:,monsoon_idx),3); Psno_monsoon(GLA_ID ~= gla_fig3_id) = NaN;
    PrecS_monsoon_y(yy) = nanmean(Psno_monsoon,'all');

    PrecS_octfeb = nansum(PSNOW_map(:,:,autumn_idx),3) +  nansum(PSNOW_map(:,:,winter_idx),3); PrecS_octfeb(GLA_ID ~= gla_fig3_id) = NaN;
    PrecS_octfeb_y(yy) = nanmean(PrecS_octfeb,'all');

    PrecS_octdec = nansum(PSNOW_map(:,:,Oct_Dec_idx),3); PrecS_octdec(GLA_ID ~= gla_fig3_id) = NaN;
    PrecS_octdec_y(yy) = nanmean(PrecS_octdec,'all');

    PrecS_janmar = nansum(PSNOW_map(:,:,Jan_Mar_idx),3); PrecS_janmar(GLA_ID ~= gla_fig3_id) = NaN;
    PrecS_janmar_y(yy) = nanmean(PrecS_janmar,'all');

    PrecS_aprjun = nansum(PSNOW_map(:,:,Apr_Jun_idx),3); PrecS_aprjun(GLA_ID ~= gla_fig3_id) = NaN;
    PrecS_aprjun_y(yy) = nanmean(PrecS_aprjun,'all');

    PrecS_julsep = nansum(PSNOW_map(:,:,Jul_Sep_idx),3); PrecS_julsep(GLA_ID ~= gla_fig3_id) = NaN;
    PrecS_julsep_y(yy) = nanmean(PrecS_julsep,'all');    

    PrecS_y(yy) = PrecS_octdec_y(yy) + PrecS_janmar_y(yy) + PrecS_aprjun_y(yy) + PrecS_julsep_y(yy);

    %Precipitations 

    Prec_spring = nansum(PRECIP_map(:,:,spring_idx),3); Prec_spring(GLA_ID ~= gla_fig3_id) = NaN;
    Prec_spring_y(yy) = nanmean(Prec_spring,'all');

    Prec_monsoon = nansum(PRECIP_map(:,:,monsoon_idx),3); Prec_monsoon(GLA_ID ~= gla_fig3_id) = NaN;
    Prec_monsoon_y(yy) = nanmean(Prec_monsoon,'all');

    Prec_octfeb = nansum(PRECIP_map(:,:,autumn_idx),3) + nansum(PRECIP_map(:,:,winter_idx),3); Prec_octfeb(GLA_ID ~= gla_fig3_id) = NaN;
    Prec_octfeb_y(yy) = nanmean(Prec_octfeb,'all');

    Prec_octdec = nansum(PRECIP_map(:,:,Oct_Dec_idx),3); Prec_octdec(GLA_ID ~= gla_fig3_id) = NaN;
    Prec_octdec_y(yy) = nanmean(Prec_octdec,'all');

    Prec_janmar = nansum(PRECIP_map(:,:,Jan_Mar_idx),3); Prec_janmar(GLA_ID ~= gla_fig3_id) = NaN;
    Prec_janmar_y(yy) = nanmean(Prec_janmar,'all');

    Prec_aprjun = nansum(PRECIP_map(:,:,Apr_Jun_idx),3); Prec_aprjun(GLA_ID ~= gla_fig3_id) = NaN;
    Prec_aprjun_y(yy) = nanmean(Prec_aprjun,'all');

    Prec_julsep = nansum(PRECIP_map(:,:,Jul_Sep_idx),3); Prec_julsep(GLA_ID ~= gla_fig3_id) = NaN;
    Prec_julsep_y(yy) = nanmean(Prec_julsep,'all');    
    
    Prec_y(yy) = Prec_spring_y(yy) + Prec_monsoon_y(yy) + Prec_octfeb_y(yy);

    %Glacier melt

    Smg_spring = nansum(SMG_map(:,:,spring_idx),3); Smg_spring(GLA_ID ~= gla_fig3_id) = NaN;
    Smg_spring_y(yy) = nanmean(Smg_spring,'all');

    Smg_monsoon = nansum(SMG_map(:,:,monsoon_idx),3); Smg_monsoon(GLA_ID ~= gla_fig3_id) = NaN;
    Smg_monsoon_y(yy) = nanmean(Smg_monsoon,'all');

    Smg_octfeb = nansum(SMG_map(:,:,autumn_idx),3)+ nansum(SMG_map(:,:,winter_idx),3); Smg_octfeb(GLA_ID ~= gla_fig3_id) = NaN;
    Smg_octfeb_y(yy) = nanmean(Smg_octfeb,'all');

    Smg_y_map =  nansum(SMG_map(:,:,(spring_idx+autumn_idx+monsoon_idx+winter_idx)==1),3);
    Smg_y_map(GLA_ID ~= gla_fig3_id) = NaN;
    Smg_y(yy) = nanmean(Smg_y_map,'all');

    % Snow melt
    
    Sms_spring = nansum(SMSm_map(:,:,spring_idx),3); Sms_spring(GLA_ID ~= gla_fig3_id) = NaN;
    Sms_spring_y(yy) = nanmean(Sms_spring,'all');

    Sms_monsoon = nansum(SMSm_map(:,:,monsoon_idx),3); Sms_monsoon(GLA_ID ~= gla_fig3_id) = NaN;
    Sms_monsoon_y(yy) = nanmean(Sms_monsoon,'all');

    Sms_octfeb = nansum(SMSm_map(:,:,autumn_idx),3) + nansum(SMSm_map(:,:,winter_idx),3); Sms_octfeb(GLA_ID ~= gla_fig3_id) = NaN;
    Sms_octfeb_y(yy) = nanmean(Sms_octfeb,'all');

    Sms_y(yy) = Sms_spring_y(yy) + Sms_monsoon_y(yy) + Sms_octfeb_y(yy);
   
    % ESN
    ESN_spring = nansum(EICE_map(:,:,spring_idx)+ESN_map(:,:,spring_idx),3); ESN_spring(GLA_ID ~= gla_fig3_id) = NaN;
    ESN_spring_y(yy) = nanmean(ESN_spring,'all');

    ESN_monsoon = nansum(EICE_map(:,:,monsoon_idx)+ESN_map(:,:,monsoon_idx),3); ESN_monsoon(GLA_ID ~= gla_fig3_id) = NaN;
    ESN_monsoon_y(yy) = nanmean(ESN_monsoon,'all');

    ESN_octfeb = nansum(ESN_map(:,:,autumn_idx),3) +nansum(EICE_map(:,:,autumn_idx),3) + ...
        nansum(ESN_map(:,:,winter_idx),3)+nansum(EICE_map(:,:,winter_idx),3); 
    ESN_octfeb(GLA_ID ~= gla_fig3_id) = NaN;
    ESN_octfeb_y(yy) = nanmean(ESN_octfeb,'all');

    %Avalanches

    AVA_y_map =  nansum(AVA_map(:,:,(spring_idx+autumn_idx+monsoon_idx+winter_idx)==1),3);
    AVA_y_map(GLA_ID ~= gla_fig3_id) = NaN;
    AVA_y(yy) = nanmean(AVA_y_map,'all');

    ESN_y(yy) = ESN_spring_y(yy) + ESN_monsoon_y(yy) + ESN_octfeb_y(yy);
   
    Net_y(yy) = PrecS_y(yy) - Smg_y(yy) - Sms_y(yy) - ESN_y(yy) + AVA_y(yy);
end 

SMB_annual = table(Years_f3, Net_y','VariableNames',{'Year','SMB'});
save([dir_fig '\GMB\SMB_' glacier '_annual.mat'],'SMB_annual');

%% Figure 2 and 3 combined with monsoon months

Xnet_y = [PrecS_monsoon_y; PrecS_spring_y; PrecS_octfeb_y; -Smg_y; -Sms_y; -ESN_y];

Hydro_year_label = strcat(string(num2str(Years_f3_all(1:end-1),'%02d ')), '/', ...
    string(num2str(Years_f3_all(2:end)-100*floor(Years_f3_all(2:end)/100),'%02d ')));

Hydro_year_label_nan = Hydro_year_label;
Hydro_year_label_nan(2:2:end) = ""; 

Hug_indices = find(Years_f3==2000):find(Years_f3==2020);
KH9_SRTM_indices = find(Years_f3==1974):find(Years_f3==2000);

Hug_period_i = '2000-01-01_2020-01-01'; % 2000-2020
GMB_HUG_MAIN_mean = Hug_gmb.(['RGI_' num2str(gla_fig3_id)]).dmdtda(categorical(cellstr(Hug_period_i)) == Hug_gmb.(['RGI_' num2str(gla_fig3_id)]).period);
GMB_HUG_MAIN_err = Hug_gmb.(['RGI_' num2str(gla_fig3_id)]).err_dmdtda(categorical(cellstr(Hug_period_i)) == Hug_gmb.(['RGI_' num2str(gla_fig3_id)]).period);
GMB_HUG_MAIN_tc = nanmean(Net_y(Years_f3>1998 & Years_f3 < 2019));
GMB_KH9_SRTM_tc = nanmean(Net_y(Years_f3>1973 & Years_f3 < 2000));

if sum(ismember([2000,2018], Years_f3))==2
disp(['GMB over the period 1999-2018: ' num2str(0.001*nanmean(Net_y(Years_f3>1998 & Years_f3 < 2018)),2) ...
    ' ' char(177) ' '  num2str(0.001*nanstd(Net_y(Years_f3>1998 & Years_f3 < 2018)),2)  ' m w.e./yr'])

disp(['GMB over the period 2018-2023: ' num2str(0.001*nanmean(Net_y(Years_f3>2017)),2) ...
   ' ' char(177) ' '   num2str(0.001*nanstd(Net_y(Years_f3>2017)),2) ' m w.e./yr'])
end 

net_lab_period = {'Icemelt','Snowmelt','Solid precipitation','Net mass balance', ' ','Total precipitation'};
fontsz = 11; p_year = Years_f3';

fi5 = figure('Renderer', 'painters', 'Position',[150.3333 99.6667 686.6667 508.6666]);
% stack positive parts
pos_x = Xnet_y'; 
pos_x(pos_x<0) = nan; 
a1 = area(pos_x.*0.001,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.5,'LineStyle','none');
colors = [colorbrewer.qual.Dark2{1, 8}(2,:)./255; colorbrewer.qual.Dark2{1, 8}(1,:)./255; colorbrewer.qual.Dark2{1, 8}(3,:)./255;...
            0.75 0.75 0.75; 0.42 0.95 1; 0.3 0.3 0.3];
for i = 1:length(a1)
    set(a1(i), 'FaceColor', colors(i,:)) 
end
hold on; 
% stack negative parts
neg_x = Xnet_y'; 
neg_x(neg_x>0) = nan; 
a2 = area( neg_x.*0.001,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.7,'HandleVisibility','off','LineStyle','none');
colors = [colorbrewer.qual.Dark2{1, 8}(2,:)./255; colorbrewer.qual.Dark2{1, 8}(1,:)./255; colorbrewer.qual.Dark2{1, 8}(3,:)./255;...
            0.75 0.75 0.75; 0.42 0.95 1; 0.3 0.3 0.3];
for i = 1:length(a1)
    set(a2(i), 'FaceColor', colors(i,:)) 
end
movmean_tot = movmean(Net_y,10);

% Fig 2. part
b1 = bar(1:length(Years_f3), Net_y.*0.001,'FaceColor',[0.4 0.4 0.4],'DisplayName','Individual years',...
    'BarWidth',0.3); hold on; grid on;
b1.FaceColor = [0.4 0.4 0.4]; b1.LineWidth = 1.3;

plot(1:5,movmean_tot(1:5).*0.001,'--r','LineWidth',1.5,'HandleVisibility','off')
p1 = plot(5:length(Years_f3)-5,movmean_tot(5:end-5).*0.001,'r','LineWidth',1.5);
yline(0,'k','LineWidth',1.3,'HandleVisibility','off')
plot((length(Years_f3)-6):(length(Years_f3)-1),movmean_tot(end-5:end).*0.001,'--r','LineWidth',1.5,'HandleVisibility','off')
ylabel('Glacier mass balance (m w.e. a^{-1})','FontSize',fontsz);

for ii = 1:(length(Years_f3)-1)
    patch([ii-0.5 ii+0.5 ii+0.5 ii-0.5], [1.93 1.93 2 2],Ta_devia_col(ii,:),'HandleVisibility','off');
end 

if isempty(Hug_indices); Hug_indices = [0 0]; end
    s1 = shadedErrorBar(Hug_indices, Hug_indices.*0+GMB_HUG_MAIN_mean, Hug_indices.*0+GMB_HUG_MAIN_err,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(3,:)./255 0.65],'LineWidth',1.4,'DisplayName','Hugonnet 2021','LineStyle','-'}); hold on
    p2 = plot(Hug_indices, Hug_indices.*0+GMB_HUG_MAIN_tc.*0.001,'Color',colorbrewer.qual.Dark2{1, 8}(3,:)./255,'LineWidth',2,'DisplayName','Hugonnet 2021 (modelled)','LineStyle','--','HandleVisibility','off'); hold on
 
if strcmp(glacier,'Parlung4') && ~isempty(KH9_SRTM_indices)
    shadedErrorBar(KH9_SRTM_indices, KH9_SRTM_indices.*0 - 0.18, KH9_SRTM_indices.*0+0.14,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(8,:)./255 0.80],'LineWidth',1.3,'DisplayName','KH-9 - SRTM','LineStyle','-'}); hold on
    plot(KH9_SRTM_indices, KH9_SRTM_indices.*0 + GMB_KH9_SRTM_tc.*0.001,'Color',[colorbrewer.qual.Dark2{1, 8}(8,:)./255],'LineWidth',2,'DisplayName','KH-9 - SRTM (modelled)','LineStyle','--','HandleVisibility','off'); hold on
elseif strcmp(glacier,'Kyzylsu') && ~isempty(KH9_SRTM_indices)
    shadedErrorBar(KH9_SRTM_indices, KH9_SRTM_indices.*0 +0.05, KH9_SRTM_indices.*0+0.2,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(8,:)./255 0.80],'LineWidth',1.3,'DisplayName','KH-9 - SRTM','LineStyle','-'}); hold on
    plot(KH9_SRTM_indices, KH9_SRTM_indices.*0 + GMB_KH9_SRTM_tc.*0.001,'Color',[colorbrewer.qual.Dark2{1, 8}(8,:)./255],'LineWidth',2,'DisplayName','KH-9 - SRTM (modelled)','LineStyle','--','HandleVisibility','off'); hold on
end

lg2 = legend(a1,'Jun-Sep snowfall','Mar-May snowfall','Oct-Feb snowfall','Box','off');
lg2.NumColumns = 3; lg2.Location = 'northwest'; fontsize(lg2,11,"points")

set(gca,'XTick',p_year(1:1:end)'-p_year(1)+1) ;set(gca,'XTickLabel',Hydro_year_label_nan); 
xtickangle(45); ylim([-2.5 2]); a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',fontsz);
xlim([0.9 length(Years_f3)+0.1])
title({[gla_fig3_label ' Glacier'], ' '})
ah1=axes('position',get(gca,'position'),'visible','off');
lg3 = legend(ah1,[a2(1,4:end) b1 p1 s1.mainLine p2],'Icemelt','Snowmelt','Evap/Sublim','Annual GMB','GMB (10yr mov-avg)',...
    '2000-2019 geodetic MB','2000-2019 modelled GMB','Box','off');%,'Total precip');
lg3.NumColumns = 2; lg3.Location = 'south'; fontsize(lg3,11,"points")
exportgraphics(fi5,[dir_fig '\T&C_Figure3_Fig2_combined_' gla_fig3_label '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Figure 2 and 3 combined with 4 seasons for snowfall

Xnet_y = [PrecS_octdec_y; PrecS_janmar_y; PrecS_aprjun_y; PrecS_julsep_y; AVA_y; -Smg_y; -Sms_y; -ESN_y];

Hydro_year_label = strcat(string(num2str(Years_f3_all(1:end-1),'%02d ')), '/', ...
    string(num2str(Years_f3_all(2:end)-100*floor(Years_f3_all(2:end)/100),'%02d ')));

Hug_indices = find(Years_f3==2000):find(Years_f3==2020);
KH9_SRTM_indices = find(Years_f3==1974):find(Years_f3==2000);

Hug_period_i = '2000-01-01_2020-01-01'; % 2000-2020
GMB_HUG_MAIN_mean = Hug_gmb.(['RGI_' num2str(gla_fig3_id)]).dmdtda(categorical(cellstr(Hug_period_i)) == Hug_gmb.(['RGI_' num2str(gla_fig3_id)]).period);
GMB_HUG_MAIN_err = Hug_gmb.(['RGI_' num2str(gla_fig3_id)]).err_dmdtda(categorical(cellstr(Hug_period_i)) == Hug_gmb.(['RGI_' num2str(gla_fig3_id)]).period);
GMB_HUG_MAIN_tc = nanmean(Net_y(Years_f3>1998 & Years_f3 < 2019));
GMB_KH9_SRTM_tc = nanmean(Net_y(Years_f3>1973 & Years_f3 < 2000));

if sum(ismember([2000,2018], Years_f3))==2
disp(['GMB over the period 1999-2017: ' num2str(0.001*nanmean(Net_y(Years_f3>1998 & Years_f3 < 2018)),2) ...
    ' +/- ' num2str(0.001*nanstd(Net_y(Years_f3>1998 & Years_f3 < 2018)),2) ' m w.e./yr'])

disp(['GMB over the period 2018-2023: ' num2str(0.001*nanmean(Net_y(Years_f3>2017)),2) ...
    '+/-' num2str(0.001*nanstd(Net_y(Years_f3>2017)),2) ' m w.e./yr'])
end

net_lab_period = {'Icemelt','Snowmelt','Solid precipitation','Net mass balance', ' ','Total precipitation'};
p_year = Years_f3';
fontsz = 12; 

fi5 = figure('Renderer', 'painters', 'Position',[50 50 779.3333 632]);
% stack positive parts
pos_x = Xnet_y'; 
pos_x(pos_x<0) = nan; 
a1 = area(pos_x.*0.001,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.7,'LineStyle','none');
colors_pos = [colorbrewer.qual.Dark2{1, 8}(3,:)./255;...
           0.4 0.2 0.8; colorbrewer.qual.Dark2{1, 8}(1,:)./255; colorbrewer.qual.Dark2{1, 8}(2,:)./255; 0.8 0.8 0.8; 0.75 0.75 0.75; 0.42 0.95 1; 0.3 0.3 0.3];
for i = 1:length(a1); set(a1(i), 'FaceColor', colors_pos(i,:)); end
hold on; 
% stack negative parts
neg_x = Xnet_y'; 
neg_x(neg_x>0) = nan; 
a2 = area( neg_x.*0.001,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.7,'HandleVisibility','off','LineStyle','none');
colors_neg = [colorbrewer.qual.Dark2{1, 8}(8,:)./255; colorbrewer.qual.Dark2{1, 8}(2,:)./255; colorbrewer.qual.Dark2{1, 8}(1,:)./255; colorbrewer.qual.Dark2{1, 8}(3,:)./255;...
            0.4 0.2 0.8; 0.75 0.75 0.75; 0.42 0.95 1; 0.3 0.3 0.3];
for i = 1:length(a1); set(a2(i), 'FaceColor', colors_neg(i,:)); end
movmean_tot = movmean(Net_y,5);

% Fill the area gaps on the left and right margins
a3 = area([0 1], [neg_x(1,:); neg_x(1,:)].*0.001,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.7,'HandleVisibility','off','LineStyle','none');
for i = 1:length(a3); set(a3(i), 'FaceColor', colors_neg(i,:),'FaceAlpha',0.5); end

a4 = area([size(neg_x,1) size(neg_x,1)+1], [neg_x(end,:); neg_x(end,:)].*0.001,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.7,'HandleVisibility','off','LineStyle','none');
for i = 1:length(a4); set(a4(i), 'FaceColor', colors_neg(i,:),'FaceAlpha',0.5); end

a5 = area([0 1], [pos_x(1,:); pos_x(1,:)].*0.001,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.7,'HandleVisibility','off','LineStyle','none');
for i = 1:length(a5); set(a5(i), 'FaceColor', colors_pos(i,:),'FaceAlpha',0.5); end

a6 = area([size(pos_x,1) size(pos_x,1)+1], [pos_x(end,:); pos_x(end,:)].*0.001,'BaseValue', 0,'FaceColor','flat','FaceAlpha',0.7,'HandleVisibility','off','LineStyle','none');
for i = 1:length(a6); set(a6(i), 'FaceColor', colors_pos(i,:),'FaceAlpha',0.5); end

% Fig 2. part
b1 = bar(1:length(Years_f3), Net_y.*0.001,'FaceColor',[0.4 0.4 0.4],'DisplayName','Individual years',...
    'BarWidth',0.3); hold on; grid on;
b1.FaceColor = [0.4 0.4 0.4]; b1.LineWidth = 1.3;

plot(1:2,movmean_tot(1:2).*0.001,'--r','LineWidth',1.5,'HandleVisibility','off')
p1 = plot(2:length(Years_f3)-1,movmean_tot(2:end-1).*0.001,'r','LineWidth',1.5);
yline(0,'k','LineWidth',1.3,'HandleVisibility','off')
plot((length(Years_f3)-1):(length(Years_f3)),movmean_tot(end-1:end).*0.001,'--r','LineWidth',1.5,'HandleVisibility','off')
ylabel('Glacier mass balance (m w.e. a^{-1})','FontSize',fontsz);

for ii = 1:length(Years_f3)
    patch([ii-0.5 ii+0.5 ii+0.5 ii-0.5], [2.43 2.43 2.5 2.5],Ta_devia_col(ii,:),'HandleVisibility','off');
end 

if isempty(Hug_indices); Hug_indices = [0 0]; end
    s1 = shadedErrorBar(Hug_indices, Hug_indices.*0+GMB_HUG_MAIN_mean, Hug_indices.*0+GMB_HUG_MAIN_err,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(3,:)./255 0.65],'LineWidth',1.4,'DisplayName','Hugonnet 2021','LineStyle','-'}); hold on
    p2 = plot(Hug_indices, Hug_indices.*0+GMB_HUG_MAIN_tc.*0.001,'Color',colorbrewer.qual.Dark2{1, 8}(3,:)./255,'LineWidth',2,'DisplayName','Hugonnet 2021 (modelled)','LineStyle','--','HandleVisibility','off'); hold on
 
if strcmp(glacier,'Parlung4') && ~isempty(KH9_SRTM_indices)
    shadedErrorBar(KH9_SRTM_indices, KH9_SRTM_indices.*0 - 0.18, KH9_SRTM_indices.*0+0.14,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(8,:)./255 0.80],'LineWidth',1.3,'DisplayName','KH-9 - SRTM','LineStyle','-'}); hold on
    plot(KH9_SRTM_indices, KH9_SRTM_indices.*0 + GMB_KH9_SRTM_tc.*0.001,'Color',[colorbrewer.qual.Dark2{1, 8}(8,:)./255],'LineWidth',2,'DisplayName','KH-9 - SRTM (modelled)','LineStyle','--','HandleVisibility','off'); hold on
elseif strcmp(glacier,'Kyzylsu') && ~isempty(KH9_SRTM_indices)
    shadedErrorBar(KH9_SRTM_indices, KH9_SRTM_indices.*0 +0.05, KH9_SRTM_indices.*0+0.2,'lineProps',{'Color',[colorbrewer.qual.Dark2{1, 8}(8,:)./255 0.80],'LineWidth',1.3,'DisplayName','KH-9 - SRTM','LineStyle','-'}); hold on
    plot(KH9_SRTM_indices, KH9_SRTM_indices.*0 + GMB_KH9_SRTM_tc.*0.001,'Color',[colorbrewer.qual.Dark2{1, 8}(8,:)./255],'LineWidth',2,'DisplayName','KH-9 - SRTM (modelled)','LineStyle','--','HandleVisibility','off'); hold on
end

lg2 = legend(a1(1:4),'Oct-Dec','Jan-Mar','Apr-Jun','Jul-Sep','Avalanches','Box','off');
lg2.NumColumns = 2; lg2.Location = 'northeast'; fontsize(lg2,11,"points")
title(lg2,'Snowfall')

Hydro_year_label_nan = Hydro_year_label;
Hydro_year_label_nan(2:2:end) = ""; 

set(gca,'XTick',p_year(1:1:end)'-p_year(1)+1) ;set(gca,'XTickLabel',Hydro_year_label_nan); 
xtickangle(45); ylim([-2.5 2.5]); a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',fontsz);
xlim([0.7 length(Years_f3)+0.3])
title({[gla_fig3_label ' Glacier'], ' '})
Ax = gca;
Ax.Box = 'off';
xline(max(xlim),'HandleVisibility','off')
yline(max(ylim),'HandleVisibility','off')

ah1=axes('position',get(gca,'position'),'visible','off');
lg3 = legend(ah1,[b1 p1 s1.mainLine p2 a2(1,6:end)],'Annual GMB','GMB (5yr mov-avg)',...
    '2000-2019 geodetic MB','2000-2019 modelled GMB','Icemelt','Snowmelt','Evap/Sublim','Box','off');%,'Total precip');
lg3.NumColumns = 2; lg3.Position = [0.15 0.12 0.46 0.125]; fontsize(lg3,11,"points")

ah2=axes('position',get(gca,'position'),'visible','off');
lg4 = legend(ah2,a1(5),'Avalanches','Box','off');
lg4.NumColumns = 1; lg4.Position = [0.09 0.83 0.30 0.05]; fontsize(lg4,11,"points")

exportgraphics(fi5,[dir_fig '\T&C_Figure3_Fig2_combined_4seasons_' gla_fig3_label '.png'],'Resolution',300,'BackgroundColor','none')

%% Number for the papers

if strcmp(glacier,'Kyzylsu')

disp(['Snowfall contribution to GMB in Oct-Dec: ' num2str(100*nanmean(PrecS_octdec_y./PrecS_y),2) '% (' num2str(nanmean(PrecS_octdec_y),3) ' mm)'] )
disp(['Snowfall contribution to GMB in Jan-Mar: ' num2str(100*nanmean(PrecS_janmar_y./PrecS_y),2) '% (' num2str(nanmean(PrecS_janmar_y),3) ' mm)'] )
disp(['Snowfall contribution to GMB in Apr-Jun: ' num2str(100*nanmean(PrecS_aprjun_y./PrecS_y),2) '% (' num2str(nanmean(PrecS_aprjun_y),3) ' mm)'] )
disp(['Snowfall contribution to GMB in Jul-Sep: ' num2str(100*nanmean(PrecS_julsep_y./PrecS_y),2) '% (' num2str(nanmean(PrecS_julsep_y),3) ' mm)'] )
disp(['Avalanches contribution to GMB: ' num2str(nanmean(AVA_y),3) ' +/- '  num2str(nanstd(AVA_y),2) ' mm w.e.'] )

Icemelt_gla_p1 = nanmean(Smg_y(Years_f3>1998 & Years_f3 < 2018));
Icemelt_gla_p2 = nanmean(Smg_y(Years_f3 > 2017));
disp(['Increase in icemelt contribution to GMB : ' num2str(Icemelt_gla_p2-Icemelt_gla_p1,3) ' mm w.e.'])

Snowmelt_gla_p1 = nanmean(Sms_y(Years_f3>1998 & Years_f3 < 2018));
Snowmelt_gla_p2 = nanmean(Sms_y(Years_f3 > 2017));

AVAL_p1 = nanmean(AVA_y(Years_f3>1998 & Years_f3 < 2018))*0.001; AVAL_std_p1 = nanstd(AVA_y(Years_f3>1998 & Years_f3 < 2018))*0.001;
AVAL_p2 = nanmean(AVA_y(Years_f3 > 2017))*0.001; AVAL_std_p2 = nanstd(AVA_y(Years_f3 > 2017))*0.001;

ET_gla_p1 = nanmean(ESN_y(Years_f3>1998 & Years_f3 < 2018));
ET_gla_p2 = nanmean(ESN_y(Years_f3 > 2017));
ET_gla_all = nanmean(ESN_y);

disp(['Avalanches contribution to GMB in 2000-2018: ' num2str(AVAL_p1,2) ' +/- '  num2str(AVAL_std_p1,2) ' mm w.e.'] )
disp(['Avalanches contribution to GMB in 2018-2023: ' num2str(AVAL_p2,2) ' +/- '  num2str(AVAL_std_p2,2) ' mm w.e.'] )

disp(['ET glacier mass losses relative to annual snowfall: ' num2str(100*nanmean(ESN_y)./nanmean(PrecS_y),2) ' %'])
disp(['Mass removed from glacier budget by evapotranspiration (including sublimation): ' num2str(nanmean(ESN_y),3) ' mm w.e./yr'])

Snow_OctDec_p1 = nanmean(PrecS_octdec_y(Years_f3>1998 & Years_f3 < 2018));
Snow_OctDec_p2 = nanmean(PrecS_octdec_y(Years_f3 > 2017));
Prec_OctDec_p1 = nanmean(Prec_octdec_y(Years_f3>1998 & Years_f3 < 2018));
Prec_OctDec_p2 = nanmean(Prec_octdec_y(Years_f3 > 2017));

Snow_JanMar_p1 = nanmean(PrecS_janmar_y(Years_f3>1998 & Years_f3 < 2018));
Snow_JanMar_p2 = nanmean(PrecS_janmar_y(Years_f3 > 2017));
Prec_JanMar_p1 = nanmean(Prec_janmar_y(Years_f3>1998 & Years_f3 < 2018));
Prec_JanMar_p2 = nanmean(Prec_janmar_y(Years_f3 > 2017));

Snow_AprJun_p1 = nanmean(PrecS_aprjun_y(Years_f3>1998 & Years_f3 < 2018));
Snow_AprJun_p2 = nanmean(PrecS_aprjun_y(Years_f3 > 2017));
Prec_AprJun_p1 = nanmean(Prec_aprjun_y(Years_f3>1998 & Years_f3 < 2018));
Prec_AprJun_p2 = nanmean(Prec_aprjun_y(Years_f3 > 2017));

Snow_JulSep_p1 = nanmean(PrecS_julsep_y(Years_f3>1998 & Years_f3 < 2018));
Snow_JulSep_p2 = nanmean(PrecS_julsep_y(Years_f3 > 2017));
Prec_JulSep_p1 = nanmean(Prec_julsep_y(Years_f3>1998 & Years_f3 < 2018));
Prec_JulSep_p2 = nanmean(Prec_julsep_y(Years_f3 > 2017));

Snow_gla_p1 = Snow_OctDec_p1 + Snow_JanMar_p1 + Snow_AprJun_p1 + Snow_JulSep_p1;
Snow_gla_p2 = Snow_OctDec_p2 + Snow_JanMar_p2 + Snow_AprJun_p2 + Snow_JulSep_p2;

Prec_gla_p1 = Prec_OctDec_p1 + Prec_JanMar_p1 + Prec_AprJun_p1 + Prec_JulSep_p1;
Prec_gla_p2 = Prec_OctDec_p2 + Prec_JanMar_p2 + Prec_AprJun_p2 + Prec_JulSep_p2;

disp(['Decrease in snow contribution to GMB : ' num2str(Snow_gla_p2-Snow_gla_p1,3) 'mm (' num2str(100*(Snow_gla_p2-Snow_gla_p1)/Snow_gla_p1,2) ' %)'])
disp(['Decrease in snow contribution to GMB in Oct-Dec : ' num2str(Snow_OctDec_p2-Snow_OctDec_p1,3) 'mm (' num2str(100*(Snow_OctDec_p2-Snow_OctDec_p1)/Snow_OctDec_p1,2) ' %)'])
disp(['Decrease in snow contribution to GMB in Jan-Mar : ' num2str(Snow_JanMar_p2-Snow_JanMar_p1,3) 'mm (' num2str(100*(Snow_JanMar_p2-Snow_JanMar_p1)/Snow_JanMar_p1,2) ' %)'])
disp(['Decrease in snow contribution to GMB in Apr-Jun : ' num2str(Snow_AprJun_p2-Snow_AprJun_p1,3) 'mm (' num2str(100*(Snow_AprJun_p2-Snow_AprJun_p1)/Snow_AprJun_p1,2) ' %)'])
disp(['Decrease in snow contribution to GMB in Jul-Sep : ' num2str(Snow_JulSep_p2-Snow_JulSep_p1,3) 'mm (' num2str(100*(Snow_JulSep_p2-Snow_JulSep_p1)/Snow_JulSep_p1,2) ' %)'])

disp(['Snowfall fraction in 2000-2017 and 2018-2023: ' num2str(Snow_gla_p1./Prec_gla_p1,2) ' and ' num2str(Snow_gla_p2./Prec_gla_p2,2)])
disp(['Snowfall fraction in 2000-2017 and 2018-2023 in Oct-Dec: ' num2str(Snow_OctDec_p1./Prec_OctDec_p1,2) ' and ' num2str(Snow_OctDec_p2./Prec_OctDec_p2,2)])
disp(['Snowfall fraction in 2000-2017 and 2018-2023 in Jan-Mar: ' num2str(Snow_JanMar_p1./Prec_JanMar_p1,2) ' and ' num2str(Snow_JanMar_p2./Prec_JanMar_p2,2)])
disp(['Snowfall fraction in 2000-2017 and 2018-2023 Apr-Jun: ' num2str(Snow_AprJun_p1./Prec_AprJun_p1,2) ' and ' num2str(Snow_AprJun_p2./Prec_AprJun_p2,2)])
disp(['Snowfall fraction in 2000-2017 and 2018-2023 Jul-Sep: ' num2str(Snow_JulSep_p1./Prec_JulSep_p1,2) ' and ' num2str(Snow_JulSep_p2./Prec_JulSep_p2,2)])

end 

%% Comparison of run with and without avalanches

if strcmp(glacier,'Kyzylsu') && ~exist('distributed','var')

sim_noav = 'ERA5Land_160824_1999_2023_PG00_Tmod0_2000mmSWEcap_noaval';
path_out_noav = ['C:\Users\jouberto\Desktop\T&C\Post-processing\Figures\Kyzylsu\' sim_noav];

SMB_noav = load([path_out_noav '/GMB/SMB_Kyzylsu_annual.mat'],'SMB_Kyz_annual'); SMB_noav = SMB_noav.SMB_Kyz_annual;

fi5 = figure('Renderer', 'painters', 'Position',[289 428.3333 654 230.6667]);
b1 = bar(1:length(Years_f3), [Net_y'.*0.001 SMB_noav.SMB.*0.001],'BarWidth',1); hold on; grid on;
set(gca,'XTick',p_year(1:1:end)'-p_year(1)+1) ;set(gca,'XTickLabel',Hydro_year_label_nan); 
xtickangle(45); ylim([-1.5 1]); a = get(gca,'XTickLabel'); set(gca,'XTickLabel',a,'fontsize',9);
xlim([0.7 length(Years_f3)+0.3])
legend('With avalanches','Without avalanches');
title('Kyzylsu Glacier'); ylabel('Glacier-wide SMB [m w.e./yr]','FontSize',10)
exportgraphics(fi5,[dir_fig '\T&C_SMB_timseries_comp_aval.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
% Numbers for the paper
nanmean(Net_y'.*0.001) 
nanstd(Net_y'.*0.001) 

disp(['Kyzylsu 2000-2023 SMB with avalanches: ' num2str(nanmean(Net_y'.*0.001),2) ' ' char(177) ' '  num2str(nanstd(Net_y'.*0.001),2) ' m w.e./yr' ])
disp(['Kyzylsu 2000-2023 SMB without avalanches: ' num2str(nanmean(SMB_noav.SMB.*0.001),2) ' ' char(177) ' '  num2str(nanstd(SMB_noav.SMB.*0.001),2) ' m w.e./yr' ])

disp(['Kyzylsu 2000-2018 vs 2018-2023 difference with avalanches: ' num2str(nanmean(Net_y(1:19)'.*0.001) - nanmean(Net_y(20:24)'.*0.001),2) ' m w.e./yr' ])
disp(['Kyzylsu 2000-2018 vs 2018-2023 difference without avalanches: ' num2str(nanmean(SMB_noav.SMB(1:19)'.*0.001) - nanmean(SMB_noav.SMB(20:24)'.*0.001),2) ' m w.e./yr' ])

end