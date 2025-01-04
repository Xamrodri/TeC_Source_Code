pSMB_period = [];

if ~exist('gla_id_smbs','var')
    gla_id_smbs = gla_id; % In case no glaciers IDs other than the main glaciers were provided
end 

if ~exist('gla_names_smb','var')
    gla_names_smb{1} = glacier;
end 

if exist('spatialize','var') && spatialize == 1 && ~exist('AVA_map','var') % cluster simulation 
    AVA_map = SSN_map.*0;
end 

% Spatial limits for SMB maps
if ~exist('y_lim_smb','var') || ~exist('x_lim_smb','var')
    y_lim_smb = [ymin_map ymax_map]; x_lim_smb =[xmin_map xmax_map];
end 

% Temporal limits for a subperiod of glacier SMB
if ~exist('date_start_smb','var') || ~exist('date_end_smb','var')
    date_start_smb = date_m(1); 
    date_end_smb = date_m(end); 
end 

% Temporal limits for two subperiod of glacier SMB profile comparison

if ~exist('date_start_smb1','var') || ~exist('date_end_smb1','var') || ~exist('date_start_smb2','var') || ~exist('date_end_smb2','var')
    date_start_smb1 = date_m(1);
    date_end_smb1 = date_m(floor(size(date_m,1)./2));

    date_start_smb2 = date_m(ceil(size(date_m,1)./2));
    date_end_smb2 = date_m(end);
end 

% GMB of the whole period for all glaciers in the catchment
GMB_period = nansum(PSNOW_map,3) + nansum(AVA_map,3) - nansum(SMSm_map,3) - nansum(ESN_map,3) ...
    - nansum(SMG_map,3) - nansum(EICE_map,3);
GMB_period(GLA_ID == 0) = NaN;
GMB_period_mean = nanmean(GMB_period.*0.001./(size(SMSm_map,3)/12),'all'); % -0.56m w.e./yr for 2015-2020 ERA5-Land

%GMB for the glacier of interest

GMB_period_maingla_mean = NaN(numel(gla_id_smbs),1);

for gg = 1:numel(gla_id_smbs)
    GMB_period_maingla = GMB_period;
    GMB_period_maingla(GLA_ID ~= gla_id_smbs(gg)) = NaN; 
    GMB_period_maingla_mean(gg) = nanmean(GMB_period_maingla.*0.001./(size(PSNOW_map,3)/12),'all'); % -0.56m w.e./yr for 2015-2020 ERA5-Land
end 

% Compute SMB elevation profile for the main glacier of interest

DEM_MAIN = flipud(DEM.Z);
DEM_MAIN(GLA_ID ~= gla_id) = NaN; % Necessary to compute main-glacier SMB elevation profile

dEL=50; % width of elevation bins
ELs_tc = nanmin(DEM_MAIN,[],'all'):dEL:nanmax(DEM_MAIN,[],'all');

for iel = 1:numel(ELs_tc)
    cur=(DEM_MAIN<(ELs_tc(iel)+dEL/2))&(DEM_MAIN>=(ELs_tc(iel)-dEL/2)); %current section of DEM
    pSMB_period(iel) = nanmean(GMB_period(cur & GLA_ID == gla_id)).*0.001; 
end 

if exist('SMB_full_period_mp','var') && year(Date(end)) > 2015
    for iel = 1:numel(ELs_tc)
        cur=(DEM_MAIN<(ELs_tc(iel)+dEL/2))&(DEM_MAIN>=(ELs_tc(iel)-dEL/2)); %current section of DEM
        pSMB_period_mp(iel) = nanmean(SMB_full_period_mp(cur & GLA_ID == gla_id)); 
        Hypso_gla_el(iel) = nansum(cur==1 & (GLA_ID == gla_id),'all').*sim_res.*sim_res;
    end 
end 

%%%% GMB of the Hugonnet period %%%%

Hugonnet_years = [2000 2005 2010 2015 2020];
Hugonnet_periods = {'2000-2005','2005-2010','2010-2015','2015-2020','2000-2020'};
Hugonnet_periods_xlx = {'2000-01-01_2005-01-01','2005-01-01_2010-01-01',...
    '2010-01-01_2015-01-01','2015-01-01_2020-01-01','2000-01-01_2020-01-01'};

% For several Hugonnet periods if they exist
Is_Hug_2000_2005 = nansum(ismember([Hugonnet_years(1) Hugonnet_years(2)],Years_no)) == 2;
Is_Hug_2005_2010 = nansum(ismember([Hugonnet_years(2) Hugonnet_years(3)],Years_no)) == 2;
Is_Hug_2010_2015 = nansum(ismember([Hugonnet_years(3) Hugonnet_years(4)],Years_no)) == 2;
Is_Hug_2015_2020 = nansum(ismember([Hugonnet_years(4) Hugonnet_years(5)],Years_no)) == 2;
Is_Hug_2000_2020 = nansum(ismember([Hugonnet_years(1) Hugonnet_years(5)],Years_no)) == 2;

Is_Hug = [Is_Hug_2000_2005 Is_Hug_2005_2010 Is_Hug_2010_2015 Is_Hug_2015_2020 Is_Hug_2000_2020];
Ind_Is_Hug = find(Is_Hug == 1);

GMB_Hug_Main_mean =  deal(NaN(numel(gla_id_smbs), length(Ind_Is_Hug)));

if ~isempty(Ind_Is_Hug)
for hh = 1:numel(Ind_Is_Hug)

    period_i = Hugonnet_periods{Ind_Is_Hug(hh)};
    year_hug_start = str2num(extractBefore(period_i,'-'));
    year_hug_end = str2num(extractAfter(period_i,'-'));
    
    m_start_hug = find(datetime(year_hug_start,2,1) == datetime(year(date_m), month(date_m), day(date_m)),1); % beginning of Hugonnet period
    m_end_hug = find(datetime(year_hug_end,1,1) == datetime(year(date_m), month(date_m), day(date_m)),1); % ending of Hugonnet period

    GMB_Hug = nansum(PSNOW_map(:,:,m_start_hug:m_end_hug),3) + nansum(AVA_map(:,:,m_start_hug:m_end_hug),3) ...
        - nansum(SMSm_map(:,:,m_start_hug:m_end_hug),3) - nansum(ESN_map(:,:,m_start_hug:m_end_hug),3) ...
        - nansum(SMG_map(:,:,m_start_hug:m_end_hug),3) - nansum(EICE_map(:,:,m_start_hug:m_end_hug),3);
    
    GMB_Hug_tc(:,:,hh) = GMB_Hug./(year_hug_end-year_hug_start); 
    GMB_Hug_Main = GMB_Hug_tc(:,:,hh);
    GMB_Hug_Main(GLA_ID == 0) = NaN;
    GMB_Hug_mean(hh) = nanmean(GMB_Hug_Main,'all').*0.001; 

    for gg = 1:numel(gla_id_smbs) 
        GMB_Hug_map = GMB_Hug_Main;
        GMB_Hug_map(GLA_ID ~= gla_id_smbs(gg)) = NaN;
        GMB_Hug_Main_mean(gg,hh) = nanmean(GMB_Hug_map,'all').*0.001;
    end 

% Compute simulated SMB profile for the Hugonnet period  
    for iel = 1:numel(ELs_tc)
        cur=(DEM_MAIN<(ELs_tc(iel)+dEL/2))&(DEM_MAIN>=(ELs_tc(iel)-dEL/2)); %current section of DEM
        pSMB_tc(iel)=nanmean(GMB_Hug_Main(cur)); % mean elevation change within that elevation band
        pSMB_std_tc(iel)=nanstd(GMB_Hug_Main(cur)); % mean elevation change within that elevation band 
%         pSMB_period(iel) = nanmean(GMB_period_maingla(cur)); 
    end 
end
end

%% Compute the Hugonnet 2021 mean GMB and uncertainty for the glacier ID provided

[GMB_HUG_MAIN_mean, GMB_HUG_MAIN_err] = deal(NaN(numel(gla_id_smbs), length(Ind_Is_Hug)));

for ii = 1:length(Ind_Is_Hug)
for gg = 1:numel(gla_id_smbs) 
    Hug_period_i = Hugonnet_periods_xlx{Ind_Is_Hug(ii)};
    GMB_HUG_MAIN_mean(gg,ii) = Hug_gmb.(['RGI_' num2str(gla_id_smbs(gg))]).dmdtda(categorical(cellstr(Hug_period_i)) == Hug_gmb.(['RGI_' num2str(gla_id_smbs(gg))]).period);
    GMB_HUG_MAIN_err(gg,ii) = Hug_gmb.(['RGI_' num2str(gla_id_smbs(gg))]).err_dmdtda(categorical(cellstr(Hug_period_i)) == Hug_gmb.(['RGI_' num2str(gla_id_smbs(gg))]).period);
end
end 

%% Glacier mass balance figure for the whole period

GMB_period(GLA_ID == 0) = 0;
cm_redblue = flipud(redblue);
cm_rdylblue = cbrewer('div', 'RdYlBu', 50); cm_rdylblue(cm_rdylblue<0)=0; cm_rdylblue(cm_rdylblue>1)=1;

SMB_show = flipud(GMB_period).*0.001./nYears;
SMB_show(~flipud(MASK) | flipud(GLA_ID)==0) = NaN;

tot_height_lat = y_lim_smb(2) - y_lim_smb(1);
x_start = y_lim_smb(1) + 0.1*(x_lim_smb(2) - x_lim_smb(1));
y_start = y_lim_smb(1) + 1*(y_lim_smb(2) - y_lim_smb(1));

% For the whole period
fi3 =figure('Renderer', 'painters', 'Position', [360 278 468.3333 420]);
% imagesc(demLons(1,:), demLats(:,1),flipud(GMB_period).*0.001./nYears,'AlphaData',flipud(MASK) & flipud(GLA_ID)>0); hold on;
imageschs(DEM,SMB_show,'colormap',cm_rdylblue,'caxis',[-5 5],'exaggerate',0.4); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.1 0.1 0.1],'LineWidth',0.8); end 
plot(catchShp_utm.X,catchShp_utm.Y,'Color',[0 0 0 0.4])
cb = colorbar;
ylabel(cb,'Glacier surface mass balance [m w.e./yr]','FontSize',12)
set(gca,'YDir','normal'); set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
title([datestr(date_m(end),'dd-mmm-yyyy') ' - ' datestr(date_m(1),'dd-mmm-yyyy')]); clim([-5 5])
xlabel(['All glacier mean SMB: ' num2str(round(GMB_period_mean,2)) ' m w.e./yr'],'FontSize',10)
for gg = 1:numel(gla_id_smbs) 
    text(x_start, y_start-gg*(tot_height_lat/20), [gla_names_smb{gg} ' SMB: ' num2str(round(GMB_period_maingla_mean(gg),2)) ...
    ' m w.e./yr'],'FontSize',7)
end
set(gca,'Color',[0.8 0.8 0.8])
ylim(y_lim_smb); xlim(x_lim_smb);
exportgraphics(fi3,[dir_fig '\GMB\' glacier '_period_GMB_' datestr(date_m(1),'yyyy-mm-dd') '-' datestr(date_m(end),'yyyy-mm-dd')...
    '.png'],'Resolution',300,'BackgroundColor','none')
SMB_full_period = flipud(GMB_period).*0.001./nYears; % m w.e./yr;
save([dir_fig '\GMB\SMB_full_period_' datestr(date_m(1),'yyyy-mm-dd') '-' datestr(date_m(end),'yyyy-mm-dd')],'SMB_full_period');
% close(gcf)

%% Glacier mass balance figure for a sub-period of interest

date_smb = date_m;% - calmonths(1);

m_start = find(date_start_smb == datetime(year(date_smb), month(date_smb), day(date_smb)),1); % beginning of GMB period
m_end = find(date_end_smb == datetime(year(date_smb), month(date_smb), day(date_smb)),1); % ending of GMB period

GMB_subperiod = nansum(PSNOW_map(:,:,m_start:m_end),3) + nansum(AVA_map(:,:,m_start:m_end),3) ...
     - nansum(SMSm_map(:,:,m_start:m_end),3) - nansum(ESN_map(:,:,m_start:m_end),3) ...
     - nansum(SMG_map(:,:,m_start:m_end),3) - nansum(EICE_map(:,:,m_start:m_end),3);

GMB_subperiod(GLA_ID == 0) = NaN;
GMB_subperiod_mean = nanmean(GMB_subperiod.*0.001./(numel(m_start:m_end)/12),'all'); % -0.56m w.e./yr for 2015-2020 ERA5-Land

SMB_show = flipud(GMB_subperiod).*0.001./(numel(m_start:m_end)/12);
SMB_show(~flipud(MASK) | flipud(GLA_ID)==0) = NaN;

% SMB of main glaciers
for gg = 1:numel(gla_id_smbs)
    GMB_subperiod_maingla = GMB_subperiod;
    GMB_subperiod_maingla(GLA_ID ~= gla_id_smbs(gg)) = NaN; 
    GMB_subperiod_maingla_mean(gg) = nanmean(GMB_subperiod_maingla.*0.001./(numel(m_start:m_end)/12),'all'); % -0.56m w.e./yr for 2015-2020 ERA5-Land
end 

tot_height_lat = y_lim_smb(2) - y_lim_smb(1);
x_start = x_lim_smb(1) + 0.06*(x_lim_smb(2) - x_lim_smb(1));
y_start = y_lim_smb(1) + 1*(y_lim_smb(2) - y_lim_smb(1));

% For the subperiod period

fi3 =figure('Renderer', 'painters', 'Position', [360 278 468.3333 420]);
imageschs(DEM,SMB_show,'colormap',cm_rdylblue,'caxis',[-5 5],'exaggerate',0.4,'colorbar',1); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.1 0.1 0.1],'LineWidth',0.8); end 
plot(catchShp_utm.X,catchShp_utm.Y,'Color',[0 0 0 0.4])
cb = colorbar; ylabel(cb,'Glacier surface mass balance [m w.e./yr]','FontSize',12)
cb.Location = 'westoutside'; %cb.Position = [0.26 0.073 0.50 0.026];
set(gca,'YDir','normal'); set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
title([datestr(date_start_smb,'dd-mmm-yyyy') ' - ' datestr(date_end_smb,'dd-mmm-yyyy') ', ' PP_dataset]); clim([-7 7])
% xlabel(['All glacier mean SMB: ' num2str(round(GMB_subperiod_mean,2)) ' m w.e./yr'],'FontSize',10)
% for gg = 1:1%numel(gla_id_smbs) 
%     text(x_start, y_start-3*gg*(tot_height_lat/20), [gla_names_smb{gg} ' SMB: ' num2str(round(GMB_subperiod_maingla_mean(gg),2)) ...
%     ' m w.e./yr'],'FontSize',8)
% end
% plot([7.085*10^5 7.105*10^5],[4.3213*10^6 4.3213*10^6],'k','LineWidth',1.3);
% text(7.088*10^5, 4.3218*10^6, '2 km', 'FontSize',10,'FontWeight','bold')
xlabel(['All (' gla_names_smb{1} ') glacier mean SMB: ' num2str(round(GMB_subperiod_mean,2)) ' (' num2str(round(GMB_subperiod_maingla_mean(1),2)) ') m w.e./yr'],'FontSize',10)
set(gca,'Color',[0.8 0.8 0.8]); axis image
ylim(y_lim_smb); xlim(x_lim_smb);
set(gca,'YAxisLocation','right')
exportgraphics(fi3,[dir_fig '\GMB\' glacier '_period_GMB_' datestr(date_start_smb,'yyyy-mm-dd') '-' datestr(date_end_smb,'yyyy-mm-dd')...
    '.png'],'Resolution',300,'BackgroundColor','none')

SMB_subperiod = SMB_show;

save([dir_fig '\GMB\SMB_distributed_' datestr(date_start_smb,'yyyy-mm-dd') '-' datestr(date_end_smb,'yyyy-mm-dd')],'SMB_subperiod');

%% For several Hugonnet periods if they exist

%plot key glaciers and find their shpfile

if ~isempty(Ind_Is_Hug)

for ii = 1:length(Ind_Is_Hug)

fi2 = figure('Renderer', 'painters', 'Position', [199.6667 211 898.0000 440.6666]);
tiledlayout(1,2,'TileSpacing','compact')
nexttile
imageschs(DEM,SMB_show,'colormap',cm_rdylblue,'caxis',[-5 5],'exaggerate',0.4,'colorbar',1); hold on;
cb = colorbar;
ylabel(cb,'Glacier surface mass balance [m w.e./yr]','FontSize',12)
set(gca,'YDir','normal'); set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
for gg = 1:numel(gla_id_smbs) 
    text(x_start, y_start-gg*(tot_height_lat/20), [gla_names_smb{gg} ' SMB: ' num2str(round(GMB_Hug_Main_mean(gg,ii),2)) ...
    ' m w.e./yr'],'FontSize',7)
end
plot(catchShp_utm.X,catchShp_utm.Y,'Color',[0 0 0 0.4])
for ss = 1:size(gla_shp,2); if ismember(gla_shp(ss).gridval, gla_id_smbs); plot(gla_shp(ss).X,gla_shp(ss).Y,'Color',[0.2 0.2 0.2],'LineWidth',0.6); end; end
set(gca,'Color',[0.8 0.8 0.8]); 
ylim(y_lim_smb); xlim(x_lim_smb);
title([glacier ' - ' Hugonnet_periods{Ind_Is_Hug(ii)} ' - T&C'])

nexttile
Hug_period_i = Hugonnet_periods_xlx{Ind_Is_Hug(ii)};
imageschs(GMB_dem_all,GMB_hug_all(:,:,Ind_Is_Hug(ii)),'colormap',cm_rdylblue,'caxis',[-5 5],'exaggerate',0.4,'colorbar',1); hold on;
set(gca,'Color',[0.6 0.6 0.6]); 
cb = colorbar;
ylabel(cb,'Surface elevation change [m/yr]','FontSize',12); hold on;
plot(catchShp_utm.X,catchShp_utm.Y,'Color',[0 0 0 0.4])
for ss = 1:size(gla_shp,2); if ismember(gla_shp(ss).gridval, gla_id_smbs); plot(gla_shp(ss).X,gla_shp(ss).Y,'Color',[0.2 0.2 0.2],'LineWidth',0.6); end; end
set(gca,'YDir','normal'); set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
for gg = 1:numel(gla_id_smbs) 
    text(x_start, y_start-gg*(tot_height_lat/20), [gla_names_smb{gg} ' SMB: ' num2str(round(GMB_HUG_MAIN_mean(gg,ii),2)) ...
    ' ' char(177) ' ' num2str(round(GMB_HUG_MAIN_err(gg,ii),2)) ' m w.e./yr'],'FontSize',7)
end
ylim(y_lim_smb); xlim(x_lim_smb);
title([glacier ' - ' Hugonnet_periods{Ind_Is_Hug(ii)} ' - Hugonnet 2021'])
exportgraphics(fi2,[dir_fig '\GMB\' glacier '_HugonnetGMB_comp' Hugonnet_periods{Ind_Is_Hug(ii)}...
    '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
end 
end

%% Compare SMB per elevation bands between two periods

date_smb = date_m;% - calmonths(1);

m_start_p1 = find(date_start_smb1 == datetime(year(date_smb), month(date_smb), day(date_smb)),1); % beginning of GMB period
m_end_p1 = find(date_end_smb1 == datetime(year(date_smb), month(date_smb), day(date_smb)),1); % ending of GMB period

m_start_p2 = find(date_start_smb2 == datetime(year(date_smb), month(date_smb), day(date_smb)),1); % beginning of GMB period
m_end_p2 = find(date_end_smb2 == datetime(year(date_smb), month(date_smb), day(date_smb)),1); % ending of GMB period

GMB_p1 = nansum(PSNOW_map(:,:,m_start_p1:m_end_p1),3) + nansum(AVA_map(:,:,m_start_p1:m_end_p1),3) ...
     - nansum(SMSm_map(:,:,m_start_p1:m_end_p1),3) - nansum(ESN_map(:,:,m_start_p1:m_end_p1),3) ...
     - nansum(SMG_map(:,:,m_start_p1:m_end_p1),3) - nansum(EICE_map(:,:,m_start_p1:m_end_p1),3);

GMB_p2 = nansum(PSNOW_map(:,:,m_start_p2:m_end_p2),3) + nansum(AVA_map(:,:,m_start_p2:m_end_p2),3) ...
     - nansum(SMSm_map(:,:,m_start_p2:m_end_p2),3) - nansum(ESN_map(:,:,m_start_p2:m_end_p2),3) ...
     - nansum(SMG_map(:,:,m_start_p2:m_end_p2),3) - nansum(EICE_map(:,:,m_start_p2:m_end_p2),3);

GMB_p1(GLA_ID ~= gla_id_smbs(1)) = NaN; GMB_p2(GLA_ID ~= gla_id_smbs(1)) = NaN;

nYears_p1 = yearfrac(date_start_smb1,date_end_smb1,0);
nYears_p2 = yearfrac(date_start_smb2,date_end_smb2,0);

for iel = 1:numel(ELs_tc)
    cur=(DEM_MAIN<(ELs_tc(iel)+dEL/2))&(DEM_MAIN>=(ELs_tc(iel)-dEL/2)); %current section of DEM
    pSMB_p1(iel) = nanmean(GMB_p1(cur & GLA_ID == gla_id)).*0.001./nYears_p1; 
    pSMB_p2(iel) = nanmean(GMB_p2(cur & GLA_ID == gla_id)).*0.001./nYears_p2; 
end 

fi2 =figure('Renderer', 'painters', 'Position', [360 278 468.3333 420]);
plot(ELs_tc, pSMB_p1,'LineWidth',1.2,'Color',(colorbrewer.qual.Set1{1, 9}(2,:))./255); grid on; hold on;
plot(ELs_tc, pSMB_p2,'LineWidth',1.2,'Color',(colorbrewer.qual.Set1{1, 9}(1,:))./255); 
view(90,-90); xlabel('Elevation [m a.s.l.]'); ylabel('SMB [m w.e./yr]')
legend([num2str(year(date_start_smb1)) '-' num2str(year(date_end_smb1))],[num2str(year(date_start_smb2)) '-' num2str(year(date_end_smb2))],'Location','NorthWest')
ylim([-6 6])
yline(0,'--k','HandleVisibility','off')
title([gla_names_smb{1} ' Glacier'])
exportgraphics(fi2,[dir_fig '\GMB\' glacier '_comb_elevSMB_' num2str(year(date_start_smb1)) '-' num2str(year(date_end_smb1)) '-' num2str(year(date_end_smb2)) '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)

%% Compare fully distributed and multipoints outputs

if exist('SMB_full_period_mp','var') && year(Date(end)) > 2015

c_fact = Hypso_gla_el.*0.001;
GMB_period_maingla_mean_mp = nanmean(SMB_full_period_mp(GLA_ID == gla_id),'all'); % -0.56m w.e./yr for 2015-2020 ERA5-Land
GMB_period_mean_mp = nanmean(SMB_full_period_mp(GLA_ID > 0),'all'); % -0.56m w.e./yr for 2015-2020 ERA5-Land


fi2 =figure('Renderer', 'painters', 'Position', [360 278 468.3333 420]);
tlc = tiledlayout(1,2,'TileSpacing','compact');
nexttile
plot(ELs_tc, pSMB_period_mp,'LineWidth',1.2); grid on; hold on;
plot(ELs_tc, pSMB_period./nYears,'LineWidth',1.2); 
view(90,-90); xlabel('Elevation [m a.s.l.]'); ylabel('SMB [m w.e./yr]')
nexttile
plot(ELs_tc, pSMB_period_mp.*c_fact,'LineWidth',1.2); grid on; hold on;
plot(ELs_tc, pSMB_period./nYears.*c_fact,'LineWidth',1.2); 
legend('Multipoints','Fully distributed','Location','NorthWest')
view(90,-90); xlabel('Elevation [m a.s.l.]'); ylabel('SMB [m^3 w.e.]')
title(tlc,[glacier ' SMB, ' datestr(date_m(1),'yyyy-mm-dd') ' to ' datestr(date_m(end),'yyyy-mm-dd')],'Fontsize',10)
exportgraphics(fi2,[dir_fig '\GMB\' glacier '_Multi_vs_fully_elevSMB_.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
%% Spatial differences in SMB

fi3 =figure('Renderer', 'painters', 'Position', [360 278 468.3333 420]);
imagesc(demLons(1,:), demLats(:,1),flipud(GMB_period).*0.001./nYears - flipud(SMB_full_period_mp),'AlphaData',flipud(MASK) & flipud(GLA_ID)>0); hold on;
for ii = 1:size(gla_shp,2), plot(gla_shp(ii).X,gla_shp(ii).Y,'Color',[0.7 0.7 0.7]); end 
plot(catchShp_utm.X,catchShp_utm.Y,'k')
cb = colorbar;
ylabel(cb,'\DeltaSMB  [m w.e./yr]','FontSize',12)
set(gca,'YDir','normal'); set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[]);
title(['Distributed - Multipoints, ' datestr(date_m(1),'dd-mmm-yyyy') ' - ' datestr(date_m(end),'dd-mmm-yyyy')]); clim([-5 5])
xlabel([glacier ' (all) glacier mean \DeltaSMB: ' num2str(round(GMB_period_maingla_mean(1) - GMB_period_maingla_mean_mp,2)) ...
    ' (' num2str(round(GMB_period_mean - GMB_period_mean_mp,2)) ') m w.e./yr'],'FontSize',10)
colormap(flipud(redblue))
set(gca,'Color',[0.8 0.8 0.8])
exportgraphics(fi3,[dir_fig '\GMB\' glacier '_DeltaGMB_FullyvsMulti_' datestr(date_m(end),'yyyy-mm-dd') '-' datestr(date_m(1),'yyyy-mm-dd')...
    '.png'],'Resolution',300,'BackgroundColor','none')
close(gcf)
end 

clear y_lim_smb x_lim_smb date_start_smb date_end_smb date_start_smb1 date_end_smb1 date_start_smb2 date_end_smb2

