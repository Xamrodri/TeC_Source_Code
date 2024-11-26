%% Derive the dh/dt parametrization from Hugonnet et al. 2021
% Achille Jouberton, 14.02.2021
% TC pre-processing step

clc; clear; close all

%% Paths
catchment_name = 'Kyzylsu';  % Langtang; Rolwaling; Kyzylsu
path_tc = 'C:\Users\jouberto\Desktop\PhD\T&C\';
path_out = [path_tc '\QGIS\TC_' catchment_name '\GMB\']; %where to save the Huss relationship;

%% Load DEMS, glacier mask and thickness from the TC input file.

%define spatial resolution
spatial_res = 100;

%define site name
SIMnm =  [ catchment_name '_' num2str(spatial_res) 'm'];

% Load T&C spatial pre-processing file
load(['C:\Users\jouberto\Desktop\PhD\T&C\OUTPUTS\dtm_' SIMnm '.mat'])

DEM_ini = flipud(DTM);
icethx = flipud(GLH);
DH_orig = flipud(DH);  % Filled Hugonnet dh/dt
GLA_ID = flipud(GLA_ID);
GLA_ID(isnan(DEM_ini)) = NaN;

% Get the glacier IDs
glacier_id_all = unique(GLA_ID); glacier_id_all(glacier_id_all == 0 | isnan(glacier_id_all)) = [];
Gla_nEl_nDH = table2struct(table(glacier_id_all,'VariableNames',{'Glacier_id'}));

% Derive the relationship between normalized elevation and normalized
% elevation change for each individual glacier

for i = 1:length(glacier_id_all)     % Langtang glacier: 1504121, i = 31

DEM = DEM_ini;
mask = GLA_ID == glacier_id_all(i) & ~isnan(DEM);

rho_g = 0.85; %glacier-wide specific gravity of ice at the glacier surface (see Huss et al 2013);
DH = DH_orig; DH(~mask) = NaN; % note that this is a rate!! per year; but works in any case 
mMB = nanmean(DH(mask)).*rho_g; %glacier-wide mass balance, in m w.e. per year
DEM(~mask) = NaN;

%%%%%%%%%%%% determine hypsometry %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gArea=sum(mask(:)); %just npixels, to normalize area!

dEL=10; %10m interval
ELs = min(DEM(mask)):dEL:max(DEM(mask));
pArea = NaN.*ELs;
pDH=pArea;
pTHX = pArea;

for iel=1:numel(ELs)
    cur=(DEM<ELs(iel)+dEL/2)&(DEM>=ELs(iel)-dEL/2); %current section of DEM
    pArea(iel)=sum(cur(:))./gArea;
    pTHX(iel)=nanmean(icethx(cur));
    pDH(iel)=nanmean(DH(cur));
end

%%%%%%%%%%% normalize and filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nELs=1-(ELs-min(DEM(mask)))/(max(DEM(mask))-min(DEM(mask)));

pDH2=movmedian(pDH,30,'omitnan');
nDH=pDH2./nanmin(pDH2);
nDH(nDH<0)=0;

if gArea*100*100*10^-6 > 3 % display curve for glacier > 3km2
figure
subplot(2,1,1)
plot(ELs,pDH);hold on
plot(ELs,pDH2);hold on
ax1=gca;
ylabel('Elevation change (m a~{-1})')
xlabel('Elevation (m)')
set(ax1,'Xdir','reverse')
grid on
legend('observed','smoothed')

ax2=subplot(2,1,2);
plot(ax2,nELs,nDH);
ylabel('Normalized elevation change (m a~{-1})')
xlabel('Normalized elevation')
set(ax2,'Ydir','reverse')
ylim([0,1])
grid on
end 

Gla_nEl_nDH(i).nEls = nELs; Gla_nEl_nDH(i).nDH = nDH;

if sum(isnan(nDH)) == length(nDH)
  Gla_nEl_nDH(i).nEls = NaN; Gla_nEl_nDH(i).nDH = NaN;  
end 

end 

% Save relationships for each glacier in a matlab structure
save([path_out 'Thinning_elev_relation'],'Gla_nEl_nDH');

