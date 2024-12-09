%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%              POINT SCALE MODEL                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
% Created on Nov 25, 2024
% Author: MAXIMILIANO RODRIGUEZ
% Code originally from: ACHILLE JOUBERTON
% Area of Study: Pyrenees
% Region: Cinca Subcatchment
% Code explanation: This code launches the Point Scale version of TC model.
% Output variables are managed by the file OUTPUT_MANAGER_DIST_LABEL.
% Some code and names still refers to previous versions of the code.
%==========================================================================

%% Clear all
clc; clear all;

%% Directories
folder_path = 'M:/19_ISTA/1_TC/3_Model_Source/2_MegaWat/'; % Put here the path of where you downloaded the repository

%% Modelling period
x1s =  "01-Nov-2022 00:00:00"; % Starting point of the simulation
x2s =  "01-Jun-2023 23:00:00"; % Last timestep of the simulation

date_start = datetime(x1s);
date_end = datetime(x2s);

%% MODEL PARAMETERS
study_name = 'Pyrenees_pointscale';
% Time step for the model
dt=3600; % [s]
dth=1; % [h]

% Integration interval for Solar variables
% Hours or fraction before and after.
t_bef = 0.5; t_aft = 0.5;
%% Pixel selection
%==========================================================================
% Depends on the point to be modelled
% Select for which pixel to run the point-scale version of T&C
% 1: AWS_OnGlacier
% 2: Pluvio
%==========================================================================
point_id = 1;     
LOC = point_id;

%% Data storing  
%==========================================================================
% How to store outputs
% Recommended for long-term simulations (> 10-20 years, otherwise data volume is too important)
%==========================================================================
output_daily = 1; 

%% Study site details
SITE ='Cinca_Mid'; 
s = 2; % ID for catchment selection
FORCING = "ERA5Land";
UTM_zone = 31; % for Spain
DeltaGMT= 1; % for Spain
outlet_names = "Cinca_Mid_OUT";
outlet_name = char(outlet_names);

%% Load DEM  
dtm_file = char("dtm_Cinca_Mid_250m.mat"); 
res = str2num(extractBetween(string(dtm_file),[SITE '_'],'m.mat')); % simulation resolution [m]
disp(strcat('Model resolution: ',num2str(res)))

Tmod = 0; % temperature modification above clean ice [°C];
Pmod = 0; % factor Pmod preciptation modification, e.g. 0.3 means 30% more precipitation at highest elevation
Z_min = 3370; % lowest elevation for linear precipitation modification (min factor -> 0)
Z_max = 5000; % highest elevation for linear precipitation modification (max factor -> Pmod)

%% Precipitation phase parametrization
% 1 = 2-threshold, 
% 2 = Ding 2017, 
% 3 = single-threshold, 
% 4 = Pomeroy 2013,
% 5 = Wang 2019, 
% 6 = Jennings 2018

parameterize_phase.OPT_Pr_Part = 2; % Choice of the precipitation phase scheme
parameterize_phase.Tmax = 2; % Upper air temperature for dual temperature threshold
parameterize_phase.Tmin = 0; % Lower air temperature for dual temperature threshold
parameterize_phase.Tconst = 2; % Air temperature for constant thresholds

% Skin layer thickness:
hSTL = 0.003; %m

% Albedo scheme choice
Albsno_method = 5; % 3 doesn't work, 4 is Brock 2000, 5 is Ding 2017

%% Create the directory where model outputs will be stored
outlocation = [folder_path,'3_Pyrenees_PointScale/4_Outputs/'];
if ~exist(outlocation, 'dir'); mkdir(outlocation); addpath(genpath(outlocation)); end

% Saving initial conditions of the model
out = strcat(outlocation,'INIT_COND_', SITE ,'_MultiPoint.mat'); % file path initial conditions

%% Dependencies
addpath(genpath([folder_path,'1_Functions'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path,'5_Common_inputs'])); % Where are distributed model set-up files (needed ? yes to load dtm)
addpath(genpath([folder_path,'3_Pyrenees_PointScale/2_Forcing'])); % Where is located the meteorological forcing and Shading matrix 
addpath(genpath([folder_path, '3_Pyrenees_PointScale/3_Inputs'])); % Add path to Ca_Data


load(dtm_file); % Distributed maps pre-processing. Useful here to get the DTM and initial snow depth
DTM = DTM_orig; % Use the full DEM in case running POI outside of mask

% Precipitation vertical gradient
Pmod_S = MASK;
rate = Pmod/(Z_max-Z_min);
Pmod_S(DTM>Z_min) = 1+rate.*(DTM(DTM>Z_min)-Z_min);

% Load CO2 data
load('Ca_Data.mat');
Ca_all = Ca;
topo = 1;


%% Impose measured albedo on glacier areas
fn_alb_elev = [SITE '_Albedo_vs_elev.mat'];

if exist(fn_alb_elev,'file')>0
disp('Using measured glacier albedo')
load([SITE '_Albedo_vs_elev'])

Afirn = DTM.*0 + 0.28;
dem_inc = DTM <= Alb_el_tt{1,2};
Afirn(dem_inc) = Alb_el_tt{1,1};

for ii = 1:size(Alb_el_tt,1)
   if ii == size(Alb_el_tt,1)
       dem_inc = DTM >= Alb_el_tt{ii,2};
       Afirn(dem_inc) = Alb_el_tt{ii,1};
   else
       dem_inc = DTM > Alb_el_tt{ii,2} & DTM < Alb_el_tt{ii+1,2};
       Afirn(dem_inc) = 0.5*(Alb_el_tt{ii,1} + Alb_el_tt{ii+1,1});
   end 
end 
Afirn(Afirn > 0.4) = 0.4; %Limit bare ice albedo to 0.4, as above it's firn.
else 
    Afirn = DTM.*0 + 0.28;   
end

%% POIs
POI = readtable([SITE '_MultiPoints.txt']); %import table with points info
[POI.LAT, POI.LON] = utm2ll(POI.LON_UTM, POI.LAT_UTM, UTM_zone);

for loc = LOC  
 %% Get location for POI
id_location = char(string(POI.Name(loc))); %id
Lat = POI.LAT(loc);
Lon = POI.LON(loc);
Zbas = DTM(POI.ROW(loc),POI.COL(loc)); % Altitude
%ij = POI.idx(loc);
ij = sub2ind(size(DTM),POI.COL(loc),POI.ROW(loc)); % Location
[i, j] = ind2sub(size(DTM), ij); % Location

%% FORCING
%==========================================================================
%
%==========================================================================
load(strcat(id_location,'.mat')); % Load forcing table for the current POI

Date_all=forcing_all.Time; 

%define period and time zone info
x1=find(date_start == Date_all,1);
x2=find(date_end == Date_all,1);


%% Displaying modelling parameters
disp(['Site selected: ' SITE])
disp(['Forcing selected: ' char(FORCING)])
disp(['Running T&C for pixel: ' id_location])
disp(['Simulation period: ' datestr(date_start) ' to ' datestr(date_end)])

%% Fetch time and do date handling
Date = Date_all(x1:x2);
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

%% Carbon data
%==========================================================================
% Load carbon data and narrow down to period
%==========================================================================

d1 = find(abs(Date_CO2-datenum(Date(1)))<1/36);
d2 = find(abs(Date_CO2-datenum(Date(end)))<1/36);
Ca=Ca_all(d1:d2); 
clear d1 d2

Oa= 210000; % Intercellular Partial Pressure Oxygen [umolO2/mol]


%% Meteorological input

%narrow down period of forcing data
forcing = forcing_all(x1:x2,:);
NN= height(forcing);%%% time Step

%height of virtual station
zatm_hourly = repmat(2.00,height(forcing),1);
zatm_surface = [18 18 2 2 18 18];
zatm_hourly_on=0;

%load all forcing data
Ameas = zeros(NN,1);
N=forcing.LWIN; Latm=forcing.LWIN;

% Precipitation
Pr=forcing.PP;
Pr(isnan(Pr))=0;
Pr(Pr<0.001)=0;

if Pmod >0
  Pr = Pr.*Pmod_S(ij);
end

Pre=forcing.PRESS;    

% 2m air temperature
Ta=forcing.TA; 

% Apply Tmod 
if (GLH(ij)>0) && (DEB_MAP(ij) < 10) % glacier, but without debris
    Ta(Ta >0) = Ta(Ta >0) + Tmod;
end 

% Wind Speed
Ws=forcing.FF; Ws(Ws < 0.01) = 0.01;

% Relative humidity
if max(forcing.RH)<= 1
    U=forcing.RH;
else
    U=forcing.RH./100;
end

% Radiation
SAD1=forcing.SAD1;SAD2=forcing.SAD2; SAB1=forcing.SAB1;SAB2=forcing.SAB2;
PARB=forcing.PARB; PARD=forcing.PARD;
alpha=0; %switch for albedo
Ameas_t=0; %albedo
Aice_meas_on_hourly=zeros(height(forcing),1); %albedo
Asno_meas_on_hourly=zeros(height(forcing),1); %albedo

% esat/ea/Ds/Tdew
a=17.27; b=237.3;
esat=611*exp(a*Ta./(b+Ta)); % Vapour pressure at saturation (Pa)
ea=U.*esat;                 % Vapour pressure (Pa)
Ds= esat - ea;              % Vapor Pressure Deficit (Pa)
Ds(Ds<0)=0; 
xr=a*Ta./(b+Ta)+log10(U);
Tdew=b*xr./(a-xr);          % Presumed dewpoint temperature (°C)
clear a b xr;

%% DING PARAMETERIZATION
% Initial daily mean values for Ding parametrization
Ta_Ding_d = nanmean(Ta(1:24));
Pre_Ding_d = nanmean(Pre(1:24));
ea_Ding_d = nanmean(ea(1:24));

%% TOPOGRAPHY
num_cell=numel(DTM);
[m_cell,n_cell]=size(DTM);
MASK = MASK.*0+1;
MASKn=reshape(MASK,num_cell,1);

if topo == 1
    %load topography data and narrow down to period
    m = matfile([SITE '_ShF_' char(FORCING) '.mat']); % ShF matrix created during pre-processing step

    x1_top=find(date_start==m.Top_Date,1);
    x2_top=find(date_end==m.Top_Date,1); 

    ShF = double(squeeze(m.ShF_S(ij,:)));
    ShF = ShF(x1_top:x2_top)';
    
    rho_g = 0.35; %%% Spatial Albedo
    zeta_S = m.zeta_Sts(x1_top:x2_top,1);
    h_S = m.h_Sts(x1_top:x2_top,1);
    clear e1 e2 Top_Date zeta_Sts H_Sts ShF_S
    SvF = m.SvF_S(ij,1); % Sky view factor at pixel ij
    Ct = m.Ct_S(i,j);
    Slo_top = m.Slo_top_S(ij,1); % Slope at pixel ij
    Aspect = m.Aspect_S(ij,1); % Aspect at pixel ij
    clear SvF_S Ct_S Aspect_S

    cos_fst = cos(atan(Slo_top)).*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;

    clear zeta_S 
  
    SAB1(sin(h_S) <= 0.10)  =  0;
    SAB2(sin(h_S) <= 0.10)  =  0;
    SAD1(sin(h_S) <= 0.10)  =  0;
    SAD2(sin(h_S) <= 0.10)  =  0;
    PARB(sin(h_S) <= 0.10)  =  0;
    PARD(sin(h_S) <= 0.10)  =  0;
    SAD1 = SAD1.*SvF + Ct.*rho_g.*((SAB1./sin(h_S)).*cos_fst + (1-SvF).*SAD1);
    SAD2 = SAD2.*SvF + Ct.*rho_g.*((SAB2./sin(h_S)).*cos_fst + (1-SvF).*SAD2);
    PARD = PARD.*SvF + Ct.*rho_g.*((PARB./sin(h_S)).*cos_fst + (1-SvF).*PARD);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SAB1 =(SAB1./sin(h_S)).*cos_fst.*ShF;
    SAB2 =(SAB2./sin(h_S)).*cos_fst.*ShF;
    PARB = (PARB./sin(h_S)).*cos_fst.*ShF;

    %correSITEions, temporary
    SAB1(SAB1<0)=0;SAB2(SAB2<0)=0;PARB(PARB<0)=0;PARD(PARD<0)=0;
    SAD1(SAD1<0)=0;SAD2(SAD2<0)=0;
    SAB1(isnan(SAB1)) = 0;SAB2(isnan(SAB2)) = 0;SAD1(isnan(SAD1)) = 0;SAD2(isnan(SAD2)) = 0;
    PARB(isnan(PARB)) = 0;PARD(isnan(PARD)) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   SvF=1; %% Sky View factor = 1 if topography is not considered
   Slo_top_S = MASKn*0;
   Slo_top=Slo_top_S(ij);
end

%% LAND COVER
ksv=reshape(VEG_CODE,num_cell,1);
%%% 1 Fir (evergr.)
%%% 2 Larch (decid.)
%%% 3 Grass C3
%%% 4 Shrub (decid.)
%%% 5 Broadleaf evergreen
%%% 6 Broadleaf deciduous
%%% 7 Rock
ksv(ksv==8)=7; %%% 8 Ice = 7 Rock

%land cover partition
cc_max = 1; %% one vegetation types
switch ksv(ij)
    case 1
        %%%% LAND COVER PARTITION  Fir - evergreen
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [1 0 0 0 0 0]>0;  
    case 2
        %%%% LAND COVER PARTITION  Larch - deciduous
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];  
        cc=length(Ccrown);%% Crown area
        II = [0 1 0 0 0 0]>0;  
    case 3
        %%%% LAND COVER PARTITION  Grass C3
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 1 0 0 0]>0;  
    case 4
        %%%% LAND COVER PARTITION  Shrub dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.1; Ccrown = [0.9];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 1 0 0]>0;  
    case 5
        %%%% LAND COVER PARTITION  broadleaf evergreen vegetation dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 0 1 0]>0;  
    case 6
        %%%% LAND COVER PARTITION  broadleaf deciduous vegetation dec.
        Cwat = 0; Curb = 0.0 ; Crock = 0.0;
        Cbare = 0.0; Ccrown = [1.0];
        cc=length(Ccrown);%% Crown area
        II = [0 0 0 0 0 1]>0;  
    case 7
        %%%% LAND COVER PARTITION  Rock
        Cwat = 0; Curb = 0.0 ; Crock = 1.0;
        Cbare = 0.0; Ccrown = [0.0];
        cc=length(Ccrown);%% Crown area
        II = [ 0 0 0 1 0 0]>0;  


    otherwise
        disp('INDEX FOR SOIL VEGETATION PARAMETER INCONSISTENT')
        return
end

zatm = max(zatm_surface(II)); %choose correct atmospheric reference height

%%% SOIL
PSAN=reshape(PSAN/100,num_cell,1); Psan = PSAN(ij); % Soil sand content at pixel ij
PCLA=reshape(PCLA/100,num_cell,1); Pcla = PCLA(ij); % Soil clay content at pixel ij
PORG=reshape(PORG/100,num_cell,1); Porg= PORG(ij); % Soil organic content at pixel ij

ms=10 ; %% 11 ; 
SOIL_TH=reshape(SOIL_TH,num_cell,1);
ms_max = 10; %% Number of soil layers

%%% DEBRIS
DEB_MAP=reshape(DEB_MAP,num_cell,1);
md_max = 10; %% % Number of debris layers

%%% Initial snow depth
SNOWD=reshape(SNOWD,num_cell,1);

% Initial snow albedo
if ~exist('SNOWALB','var')
    SNOWALB = SNOWD;
    SNOWALB(SNOWD>0) = 0.6;
else 
    SNOWALB=reshape(SNOWALB,num_cell,1);
end 

%% SAVING INITIAL CONDITIONS AND PARAMETERS 
% (run this only once in MultiPoint model!)
if exist(out, 'file') == 2
load(out);
else
INIT_COND_v2(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,m.Slo_top_S,ksv,Ca,SNOWD,SNOWALB,out);
load(out);
end

%% RUN MODEL
%==========================================================================
% PARAM_IC: Define parameter file
% MAIN_FRAME: Contains the model
%==========================================================================
PARAM_IC = strcat(folder_path,'3_Pyrenees_PointScale/3_Inputs/MOD_PARAM_Multipoint.m');
MAIN_FRAME; % Launch the main frame of T&C. Most of the things happen in this line of code

%% Output manager
%==========================================================================
%
%==========================================================================
Date_R = char(Date); % make date usable for R

Param_t = table(Lat,Lon,Zbas,dbThick,'VariableNames',{'Lat','Lon','Zbas','dbThick'});
Param_t = [Param_t, struct2table(SnowIce_Param), struct2table(Deb_Par)];
Param_t = rows2vars(Param_t);
Param_t = renamevars(Param_t,{'OriginalVariableNames','Var1'},{'Parameter','Value'});

%%post-compute sublimation from ESN
SSN = ESN.*(Ts<0);

% Here I manually choose the T&C outputs I want to save at each point. 

Outputs_t = table(Date,EICE,ESN,SND,SWE,...
Ta,Ws,U,N,SAD1+SAD2+SAB1+SAB2,Pre,Pr,Pr_sno,ALB,Smelt,Imelt,SSN,ICE,ET,ros,'VariableNames',{ ...
'Date','EICE','ESN','SND','SWE',...
'Ta','Ws','U','N','Rsw',...
'Pre','Pr','Pr_sno','Albedo','Smelt','Imelt','SSN','ICE','ET','ros'});

% For daily outputs
% If I want to only save daily aggregated output
if output_daily == 1

Outputs_tt = table2timetable(Outputs_t);

Outputs_ds = retime(Outputs_tt,'daily',@nansum);
Outputs_dm = retime(Outputs_tt,'daily',@nanmean);

Outputs_d = Outputs_dm; 
Outputs_d.Pr = Outputs_ds.Pr;
Outputs_d.Pr_sno = Outputs_ds.Pr_sno;

Outputs_t = timetable2table(Outputs_d);

end 

writetable(Outputs_t, strcat(outlocation,id_location,'_results.txt'))
writetable(Param_t, strcat(outlocation,id_location,'_param.txt') )
end 

%% Run evaluation
%==========================================================================
%
%==========================================================================

path_evaluation = [folder_path '3_Pyrenees_PointScale/1_Evaluation_data/'];

%{
if point_id == 1

% On-glacier AWS data
AWS = load([path_evaluation 'AWS_Kyzylsu_15min_2021-06-28_2023-09-12.mat']);
AWS = AWS.AWS_tt_15;

SMB_aws = Outputs_tt.ICE.*0.001./0.91 + Outputs_tt.SND ...
    - Outputs_tt.ICE(find(Outputs_tt.Date>datetime(2021,6,26),1)).*0.001./0.91 ...
    - Outputs_tt.SND(find(Outputs_tt.Date>datetime(2021,6,26),1));

fi2 = figure('Renderer', 'painters', 'Position',[141.6667 244.3333 460.0000 370.6667]);
plot(AWS.Date(hour(AWS.Date) == 11), AWS.HS(hour(AWS.Date) == 11).*0.01 - AWS.HS(find(AWS.Date>datetime(2021,6,26),1))*0.01,'k','LineWidth',1.1); hold on; grid on;
plot(Outputs_tt.Date, SMB_aws,'Color',[1 0 0 0.7],'LineWidth',1.1)
% plot(SMB_aws_noDB.Time, SMB_aws_noDB.DH_AWS,'Color',[0 0 1 0.7],'LineWidth',1.1)
xlim([datetime(2021,7,1) datetime(2023,8,31)])
% legend('Measurement','Model','Model (without debris)','location','SouthWest')
legend('Measurement','Model','location','SouthWest')
title('On-glacier AWS'); ylabel('Surface elevation change [m]')
% exportgraphics(fi2,[outlocation '\T&C_AWS_validation_new.png'],'Resolution',300,'BackgroundColor','none')

end 

if point_id == 2

% Load pluviometer data 
Pluvio = load([path_evaluation '\Pluviostation_Kyzylsu_2021-07-04_2023-09-13.mat']);
Pluvio = Pluvio.Pluvio_all;
SND_meas = 3.90 - Pluvio.DIST_corr(hour(Pluvio.Date) == 12);
SND_meas(SND_meas>1.5) = NaN;


fi2 = figure('Renderer', 'painters', 'Position',[245.6667 411.6667 595.3333 286.3333]);
plot(Pluvio.Date(hour(Pluvio.Date) == 12), SND_meas,'k','LineWidth',1.1); hold on; grid on;
plot(Outputs_tt.Date, Outputs_tt.SND,'Color',[1 0 0 0.7],'LineWidth',1.1)
% plot(SMAP_Pluvio.Date, SMAP_Pluvio.SND_SMAP,'Color',[0 0 1 0.7],'LineWidth',0.8);
xlim([datetime(2021,7,1) datetime(2023,7,31)])
legend('Obs','T&C','location','NorthWest')
ylim([0 2]); ylabel('Snow depth [m]');
title('Pluviometer')
% exportgraphics(fi2,[outlocation '\T&C_SnowDepth_validation_AlbBrock.png'],'Resolution',300,'BackgroundColor','none')

end

%Max plot
%{
fi3 = figure('Renderer', 'painters', 'Position',[141.6667 244.3333 460.0000 370.6667]);
plot(out_my.Date, out_my.WaterTable); hold on; grid on;
ylim([0 2000])
%}

%}