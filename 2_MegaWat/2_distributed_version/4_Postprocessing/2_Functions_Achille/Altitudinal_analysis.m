% Compute mean major fluxes per elevation band (snowmelt, fractional snow
% cover, precipitation ratio, submlimation

if ~exist('date_el','var')
    date_el = date_m; % 'catchment','glaciers';
end 

if ~exist('volumes','var')
    volumes = 1; % 'catchment','glaciers';
end 

% Where to store altitudinal figures

dir_fig_alti = [dir_fig '\Altitudinal_figures'];
if ~exist(dir_fig_alti, 'dir'); mkdir(dir_fig_alti); end 

dEL=100; % width of elevation bins
ELs = nanmin(DTM,[],'all'):dEL:nanmax(DTM,[],'all');

Hydro_years = year(date_el(find(month(date_el) == 9 & day(date_el) == 1,1))):1:year(date_el(find(month(date_el) == 9 & day(date_el) == 1,1,'last')));
Years_el = Hydro_years;
[fsc_el_y, SR_el_y, Pr_liq_el_y, Pr_sno_el_y, Smelt_el_y, Imelt_el_y ,ESN_el_y] = deal(NaN(numel(ELs),numel(Years_el)));

if ~exist('middle_year','var')
    middle_year = Hydro_years(floor(numel(Hydro_years)./2)); % For period comparison analysis
end 

month_m = date_m-calmonths(1);

for yy = 1:length(Hydro_years)-1

    ind_start = find(year(month_m) == Hydro_years(yy) & month(month_m) == 9 & day(month_m) == 1,1);
    ind_end =   find(year(month_m) == Hydro_years(yy+1) & month(month_m) == 9 & day(month_m) == 1,1);

    Pr_y_map = nansum(PRECIP_map(:,:,ind_start:ind_end),3); Pr_y_map(MASK ~=1) = NaN;
    Pr_sno_y_map = nansum(PSNOW_map(:,:,ind_start:ind_end),3); Pr_sno_y_map(MASK ~=1) = NaN;
    Pr_liq_y_map = nansum(PRAIN_map(:,:,ind_start:ind_end),3); Pr_liq_y_map(MASK ~=1) = NaN;
    
    Smelt_y_map = nansum(SMSm_map(:,:,ind_start:ind_end),3); Smelt_y_map(MASK ~=1) = NaN;
    Imelt_y_map = nansum(SMG_map(:,:,ind_start:ind_end),3); Imelt_y_map(MASK ~=1) = NaN;
    Esn_y_map = nansum(ESN_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps
    SWE_y_map = nanmean(SWEm_map(:,:,ind_start:ind_end),3); % Annual precipitation sum maps

    for iel = 1:numel(ELs)
        cur=(DTM<(ELs(iel)+dEL/2))&(DTM>=(ELs(iel)-dEL/2)); %current section of DEM
        Hypso_el(iel) = nansum(cur==1,'all').*sim_res.*sim_res;
        Hypso_gla_el(iel) = nansum(cur==1 & (GLA_ID > 0),'all').*sim_res.*sim_res;
%         Pr_y(iel,yy) = nanmean(Pr_y_map(cur)); 
        Pr_sno_el_y(iel,yy) = nanmean(Pr_sno_y_map(cur)); 
        Pr_liq_el_y(iel,yy) = nanmean(Pr_liq_y_map(cur)); 
        SR_el_y(iel,yy) = nanmean(Pr_sno_y_map(cur)./Pr_y_map(cur)); 
        Smelt_el_y(iel,yy) = nanmean(Smelt_y_map(cur));
        Imelt_el_y(iel,yy) = nanmean(Imelt_y_map(cur));
        ESN_el_y(iel,yy) = nanmean(Esn_y_map(cur));
    end 
end 

Pr_el_y = Pr_liq_el_y + Pr_sno_el_y;
%%

if Hydro_years(1) < 2000
    mid_year_id = find(Hydro_years==2000);
else
    mid_year_id = find(Hydro_years==2018); floor(size(Pr_sno_el_y,2)./2);
end 
cmap = [redblue(numel(Years_el)) ones(numel(Years_el),1).*0.3];

if volumes == 1
    c_fact = Hypso_el.*0.001; volume_lab = 'volumetric'; unit_lab = 'm^3';
else 
    c_fact = 1; volume_lab = 'flux'; unit_lab = 'mm w.e.';
end 

fi3 =figure('Renderer', 'painters', 'Position', [209 251.6667 906.6667 416]);
tiledlayout(1,3,"TileSpacing","compact")
nexttile
p1 = plot(ELs, Pr_sno_el_y.*c_fact','LineWidth',0.7); grid on; hold on;
for ii = 1:numel(Years_el)
    p1(ii,1).Color = cmap(ii,:);
end 
plot(ELs, nanmean(Pr_sno_el_y(:,1:mid_year_id).*c_fact',2),'LineWidth',1.2,'Color',[0 0 1])
plot(ELs, nanmean(Pr_sno_el_y(:,mid_year_id+1:end).*c_fact',2),'LineWidth',1.2,'Color',[1 0 0])
plot(ELs, ESN_el_y.*c_fact','HandleVisibility','off','LineWidth',0.8,'Color',[0.8 0.8 0.8]); grid on; hold on;
plot(ELs, nanmean(ESN_el_y.*c_fact',2), 'LineWidth',1,'Color',[0.3 0.3 0.3],'HandleVisibility','off')
view(90,-90); ylabel(['Snowfall [' unit_lab ']'],'FontSize',11); xlabel('Elevation [m a.s.l.]','FontSize',11);
nexttile
p1 = plot(ELs, SR_el_y,'','HandleVisibility','off','LineWidth',0.8); grid on; hold on;
for ii = 1:numel(Years_el)
    p1(ii,1).Color = cmap(ii,:);
    if ii == 1; p1(ii,1).HandleVisibility = "on"; end
end 
plot(ELs, nanmean(SR_el_y(:,1:mid_year_id),2),'LineWidth',1.2,'Color',[0 0 1])
plot(ELs, nanmean(SR_el_y(:,mid_year_id+1:end),2),'LineWidth',1.2,'Color',[1 0 0])
ylabel('Snowfall fraction [-]','FontSize',11);
legend('Individual hydrological years',[num2str(Hydro_years(1)) '-' num2str(Hydro_years(mid_year_id)) ' average'],...
    [num2str(Hydro_years(mid_year_id+1)) '-' num2str(Hydro_years(end)) ' average'])
view(90,-90); set(gca,'XTickLabels',[]); 
colormap(cmap(:,1:3))
cb = colorbar('location','northoutside'); cb.Ticks = [0 1]; cb.TickLabels = {num2str(Years_el(1)), num2str(Years_el(end))};
pos = cb.Position;
pos(4) = 0.2*pos(4); pos(2) = pos(2) +0.07;
cb.Position = pos; cb.Box = 0;

nexttile
p1 = plot(ELs, Smelt_el_y.*c_fact','LineWidth',0.8); grid on; hold on;
for ii = 1:numel(Years_el)
    p1(ii,1).Color = cmap(ii,:);
end 
plot(ELs, nanmean(Smelt_el_y(:,1:mid_year_id).*c_fact',2),'LineWidth',1.2,'Color',[0 0 1])
plot(ELs, nanmean(Smelt_el_y(:,mid_year_id+1:end).*c_fact',2),'LineWidth',1.2,'Color',[1 0 0])
view(90,-90); set(gca,'XAxisLocation','top'); 
ylabel(['Snowmelt [' unit_lab ']'],'FontSize',11); xlabel('Elevation [m a.s.l.]','FontSize',11)
exportgraphics(fi3,[dir_fig_alti '\' glacier '_altitudinal_snowfallmelt_' volume_lab '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Same but with sublimation instead of snowfall ratio

cmap = [redblue(numel(Years_el)) ones(numel(Years_el),1).*0.3];
% volumes = 1;

% if volumes == 1
     c_fact = Hypso_el.*0.001; volume_lab = 'volumetric'; unit_lab = 'm^3';
% else 
 %   c_fact = 1; volume_lab = 'flux'; unit_lab = 'mm w.e.';
% end 

fi3 =figure('Renderer', 'painters', 'Position', [209 251.6667 906.6667 416]);
tiledlayout(1,3,"TileSpacing","compact")
nexttile
p1 = plot(ELs, Pr_liq_el_y.*c_fact','LineWidth',0.8); grid on; hold on;
for ii = 1:numel(Years_el)
    p1(ii,1).Color = cmap(ii,:);
end 
plot(ELs, nanmean(Pr_liq_el_y(:,1:mid_year_id).*c_fact',2),'LineWidth',1,'Color',[0 0 1])
plot(ELs, nanmean(Pr_liq_el_y(:,mid_year_id+1:end).*c_fact',2),'LineWidth',1,'Color',[1 0 0])
view(90,-90); ylabel(['Rainfall [' unit_lab ']'],'FontSize',11); xlabel('Elevation [m a.s.l.]','FontSize',11);
ylim([0 16*10^6])
nexttile
p2 = plot(ELs, Imelt_el_y.*c_fact','HandleVisibility','off','LineWidth',0.8); grid on; hold on;
for ii = 1:numel(Years_el)
    p2(ii,1).Color = cmap(ii,:);
    if ii == 1; p2(ii,1).HandleVisibility = "on"; end
end 
plot(ELs, nanmean(Imelt_el_y(:,1:mid_year_id).*c_fact',2),'LineWidth',1,'Color',[0 0 1])
plot(ELs, nanmean(Imelt_el_y(:,mid_year_id+1:end).*c_fact',2),'LineWidth',1,'Color',[1 0 0])

ylabel(['Icemelt [' unit_lab ']'],'FontSize',11);
legend('Individual hydrological years',[num2str(Hydro_years(1)) '-' num2str(Hydro_years(mid_year_id)) ' average'],...
    [num2str(Hydro_years(mid_year_id+1)) '-' num2str(Hydro_years(end)) ' average'])
view(90,-90); set(gca,'XTickLabels',[]); 
colormap(cmap(:,1:3))
cb = colorbar('location','northoutside'); cb.Ticks = [0 1]; cb.TickLabels = {num2str(Years_el(1)), num2str(Years_el(end))};
pos = cb.Position;
pos(4) = 0.2*pos(4); pos(2) = pos(2) +0.07;
cb.Position = pos; cb.Box = 0;
ylim([0 16*10^6])
nexttile
p1 = plot(ELs, Smelt_el_y.*c_fact','LineWidth',0.8); grid on; hold on;
for ii = 1:numel(Years_el)
    p1(ii,1).Color = cmap(ii,:);
end 
plot(ELs, nanmean(Smelt_el_y(:,1:mid_year_id).*c_fact',2),'LineWidth',1,'Color',[0 0 1])
plot(ELs, nanmean(Smelt_el_y(:,mid_year_id+1:end).*c_fact',2),'LineWidth',1,'Color',[1 0 0])
plot(ELs, ESN_el_y.*c_fact','HandleVisibility','off','LineWidth',0.8,'Color',[0.8 0.8 0.8]); grid on; hold on;
plot(ELs, nanmean(ESN_el_y.*c_fact',2), 'LineWidth',1,'Color',[0.3 0.3 0.3],'HandleVisibility','off')
view(90,-90); set(gca,'XAxisLocation','top');
ylim([0 16*10^6])
ylabel(['Snowmelt [' unit_lab ']'],'FontSize',11); xlabel('Elevation [m a.s.l.]','FontSize',11)
exportgraphics(fi3,[dir_fig_alti '\' glacier '_altitudinal_snowfallsublmimation_' volume_lab '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Figure xx Kyzylsu paper : Snowfall and snowfall fraction + hypsometry

mid_year_id = find(Hydro_years==middle_year);

cmap = [redblue(numel(Years_el)) ones(numel(Years_el),1).*0.3];
color_2000 = (colorbrewer.qual.Set1{1, 9}(1,:))./255;
color_1970 = (colorbrewer.qual.Set1{1, 9}(2,:))./255;

fi3 =figure('Renderer', 'painters', 'Position',[237.6667 221.6667 524.0000 385.3333]);
tiledlayout(1,2,"TileSpacing","tight")
nexttile
p1 = plot(ELs, Pr_sno_el_y,'LineWidth',0.4,'LineStyle','-','HandleVisibility','off'); grid on; hold on;
for ii = 1:numel(Years_el)
    if Years_el(ii) < Years_el(mid_year_id)
        p1(ii,1).Color = [color_1970 0.25];
    else 
        p1(ii,1).Color = [color_2000 0.25];
    end
    if ii == 1; p1(ii,1).HandleVisibility = "on"; end
end 
ylim([0 2000]); yticks([0:500:2000])
plot(ELs, nanmean(Pr_sno_el_y(:,1:mid_year_id),2),'LineWidth',1.4,'Color',color_1970,'LineStyle','-','HandleVisibility','on')
plot(ELs, nanmean(Pr_sno_el_y(:,mid_year_id+1:end),2),'LineWidth',1.4,'Color',color_2000,'LineStyle','-','HandleVisibility','on')
view(90,-90); ylabel(' Snowfall [mm]','FontSize',11); xlabel('Elevation [m a.s.l.]','FontSize',11);
xlim([min(ELs) max(ELs)])
lg1 = legend('Hydrological years',[num2str(Hydro_years(1)) '-' num2str(Hydro_years(mid_year_id)) ' average'],...
    [num2str(Hydro_years(mid_year_id)) '-' num2str(Hydro_years(end)) ' average'],...
    'Catchment hyspometry','Glacier hypsometry','Location','southeast');
lg1.FontSize = 8.5; lg1.NumColumns = 3; 

nexttile
p1 = plot(ELs, SR_el_y,'','HandleVisibility','off','LineWidth',0.4); grid on; hold on;
for ii = 1:numel(Years_el)
    if Years_el(ii) < Years_el(mid_year_id)
        p1(ii,1).Color = [color_1970 0.25];
    else 
        p1(ii,1).Color = [color_2000 0.25];
    end
%     if ii == 1; p1(ii,1).HandleVisibility = "on"; end
end 
plot(ELs, nanmean(SR_el_y(:,1:mid_year_id),2),'LineWidth',1.3,'Color',color_1970,'HandleVisibility','off')
plot(ELs, nanmean(SR_el_y(:,mid_year_id+1:end),2),'LineWidth',1.3,'Color',color_2000,'HandleVisibility','off')

area(ELs, 3*(Hypso_el./nansum(Hypso_el)), 'FaceColor',[0.6 0.6 0.6],'FaceAlpha',0.3,'EdgeColor','none')
area(ELs, 3*(Hypso_gla_el./nansum(Hypso_el)), 'FaceColor',[121, 208, 235]./255,'FaceAlpha',0.5,'EdgeColor','none')

ylabel('Snowfall fraction [-]','FontSize',11); ylim([0 1])
lg2 = legend('Catchment hyspometry','Glacier hypsometry','Location','northwest');
lg2.FontSize = 8.5; lg2.NumColumns = 1; lg2.Position
yticks([0:0.2:1]); 
view(90,-90); set(gca,'XTickLabels',[]); 
xlim([min(ELs) max(ELs)])
lg1.Position = [0.12   0.94    0.8467    0.0484];
lg2.Position = [0.56    0.8095    0.3200    0.0878];
exportgraphics(fi3,[dir_fig_alti '\' glacier '_altitudinal_snowfraction_hypsometry.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
