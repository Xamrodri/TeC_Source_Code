%% Apply the dh/dt - elevation relationship to redistribute glacier mass changes
%according to thinning profiles observed in Hugonnet et al. 2021
% Achille Jouberton, 14.02.2021

clc; clear; close all

%% Path 

catchment_name = 'Kyzylsu';  % Langtang; Rolwaling; Kyzylsu
spatial_res = 100;
path_tc = 'C:\Users\jouberto\Desktop\PhD\T&C\';
addpath([path_tc 'Huss_flowparam_TC']);
path_func = [path_tc 'Functions'];
path_relev = [path_tc '\QGIS\TC_' catchment_name '\GMB\']; 
addpath(path_func)

%% Load Langtang T&C GMB outputs  (will have to be modified by direct T&C GMB output)

addpath(genpath('C:\Users\jouberto\Desktop\PhD\T&C\Langtang\Pascal_spatial_outputs'));
Dir_tc = dir('C:\Users\jouberto\Desktop\PhD\T&C\Langtang\Pascal_spatial_outputs\*.mat');

for i = 1:size(Dir_tc,1)

    fn  = Dir_tc(i).name;
    fn_date = erase(fn,["OUTPUT_v20211221_SPATIAL_",".mat"]);
    hour_no(i) = str2num(fn_date);

load(fn,'ICE_D_spatial','ICE_spatial','SWE_spatial','IP_wc_spatial','SP_wc_spatial')

ICE_D = reshape(ICE_D_spatial, 291, 311);   %mm 
ICE = reshape(ICE_spatial, 291, 311);       %mm 
SWE = reshape(SWE_spatial,291,311);         %mm 
IP_wc = reshape(IP_wc_spatial,291,311);     %mm 
SP_wc = reshape(SP_wc_spatial,291,311);     %mm 

SMB_int = ICE + SWE + IP_wc + SP_wc;
SMB_int(291,:) = [];
SMB_int(:,306:311) = [];

ICE(291,:) = [];
ICE(:,306:311) = [];
SMB(:,:,i) = SMB_int;

end

GMB_1819 = flipud(SMB(:,:,26) - SMB(:,:,12));
gMB_all = GMB_1819*0.001;

%% Load inputs

SIMnm =  [ catchment_name '_' num2str(spatial_res) 'm'];
load([path_tc 'OUTPUTS\dtm_' SIMnm '.mat'],'GLH','DTM_orig','DTM','GLA_ID','DH') 
load([path_relev 'Thinning_elev_relation'],'Gla_nEl_nDH'); % Load all glaciers'curve of normalized elevation change per normalized elevation

icethx_ini = flipud(GLH);  % initial ice thickness
mask_gla_ini = icethx_ini > 0;  % initial glacier mask
DEM_ini = flipud(DTM);
GLA_ID = flipud(GLA_ID);   %  glacier IDs based on RGI 6.0
GLA_ID(isnan(DEM_ini)) = NaN; % remove glacier ID outside of catchment's mask

glacier_id_all = unique(GLA_ID); glacier_id_all(glacier_id_all == 0 | isnan(glacier_id_all)) = [];

rho_g = 0.85; %glacier-wide specific gravity of ice at the glacier surface (see Huss et al 2013); 

%% Double loop, the large one is for years (will be removed in T&C I suppose)
% The inner loop is to go through each glacier of the RGI within the domain

DEM = DEM_ini; 
THX = flipud(GLH);

% Loop for timestep (annually)
for h = 1:85  % The test was for 85 years
    % Start the loop for individual glaciers
for i = 1:length(glacier_id_all)    % Langtang glacier: 1504121, i = 31

  if isnan(Gla_nEl_nDH(i).nEls)
     continue
  else

mask = GLA_ID == glacier_id_all(i) & THX > 0; % mask of each individual glacier
mMB = -0.4 + 0.2-rand(1);%,nanmean(gMB_all(mask));% + 0.2-rand(1); % glacier-wide MB (in m w.e. per year)

THX_i = THX; THX_i(~mask) = NaN;  % ice thickness for the glacier i
relDH = Gla_nEl_nDH(i).nDH;

[newDHg,newDEM] = apply_Huss_Hug(relDH,mMB,DEM,THX_i);

% The following lines are needed to deal with ice disappearance

    THX_i(THX_i <= 0)=0;     
    cTHXt=THX_i+newDHg;  %
    t1=cTHXt<0;  % Identify areas where SMB is more negative that ice thickness left
    newDHg(t1)=-THX_i(t1); %Set the maximum SMB to the ice thickness left

%     cSMB = gMB_all; cSMB(~mask) = NaN;
%     t2=cSMB>0; 
%     newDHg(t2&(cTHXt<=0))=0; % Orginally made to ensure that positive SMB areas remain glacier, with thickness of xx m..
%     %but here set to 0 (the glacier can retreat even in zones of positive SMB).
%     newDHg=newDHg.*mask;

    %rescale cDH to conserve mass (loss) after above corrections
    % This is not strictly mass conversative, but it also makes sense since
    % the disappareance of ice reduces the full mass loss theorically caused by SMB
    % (sometimes the ice would disappear before the end of the timestep)
    
   % newDHg(~(t1|t2))=newDHg(~(t1|t2)).*mMB./nanmean(newDHg(mask&~(t1|t2)));
    newDHg(~(t1))=newDHg(~(t1)).*mMB./nanmean(newDHg(mask&~(t1)));
    
    newDHg(isnan(newDHg)) = 0;  
    newDHg(isinf(newDHg)) = 0; 
    THX(mask)=THX_i(mask)+ newDHg(mask);  % Updated ice thickness
    DEM(mask)=DEM(mask)+newDHg(mask);   % Updated surface DEM
    
   end
end
THX_y(:,:,h) = THX; % Ice thickness of each individual years.
h
end

%% Check everything worked more or less fine

figure
subplot(1,2,1)
imagesc(gMB_all,'AlphaData',~isnan(DEM_ini))
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])
caxis([-3,3]); set(gca,'Color',[0.7 0.7 0.7])
title('2018-2019 T&C SMB'); colorbar
subplot(1,2,2)
imagesc(THX_y(:,:,1) - icethx_ini,'AlphaData',~isnan(DEM_ini))
 colormap(flipud(redblue)); %ylabel(c,'[m w.e]','FontSize',11)
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])
caxis([-3,3]); set(gca,'Color',[0.7 0.7 0.7])%%
title('GMB after dh/dt correction'); c= colorbar;
ylabel(c,'[m w.e.]','Fontsize',11)

%% Ice thickness evolution difference

fi3 = figure('Position',[86.3333 310.3333 890 308.0000]);
a = subplot(1,3,1);
imagesc(icethx_ini,'AlphaData',icethx_ini > 0 & ~isnan(DEM_ini))
set(gca,'Color',[0.7 0.7 0.7]); caxis([0 300])
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])
set(a,'Position',[0.1300 0.1100 0.2134 0.8150])
title('Initial')

b = subplot(1,3,2);
imagesc(THX_y(:,:,35),'AlphaData',THX_y(:,:,35) > 0 & ~isnan(DEM_ini))
set(gca,'Color',[0.7 0.7 0.7]); caxis([0 300])
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])
set(b,'Position',[0.3808 0.1100 0.2134 0.8150])
title('2050')

c = subplot(1,3,3);
imagesc(THX,'AlphaData', THX > 0 & ~isnan(DEM_ini))
set(gca,'Color',[0.7 0.7 0.7]); caxis([0 300])
set(gca,'XTickLabels',[]); set(gca,'YTickLabels',[])
set(c,'Position', [0.6316 0.1100 0.2134 0.8150])
title('2100'); c = colorbar;
set(c,'Position',[0.8590 0.1095 0.0117 0.8159]); 
ylabel(c,'Ice thickness [m w.e.]','Fontsize',11)

print(fi3,[path_tc 'Figures\' catchment_name '_test_dhparam_GMB_100years'],'-dpng','-r300')
