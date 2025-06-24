%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C) - ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%
%%%%%%%%%%%%%%              POINT SCALE MODEL                   %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% AUTHOR INFO AND STUDY SITE
%==========================================================================
%{
Created on Nov 25, 2024
Author: MAXIMILIANO RODRIGUEZ
Code originally from: ACHILLE JOUBERTON
Area of Study: Apennines - Tiber River

Code explanation: 
     This code launches the Point Scale version of TC model.
     This function runs only one point. There is no for loop.
     
     

%}
%==========================================================================

%% Debugger
%Point = "VelinoCluster95";


%% As a Function 
function result=run_Point_Pro(root, outlocation, run_folder, Point, POI, ksv, date_start, date_end)


%% Directories
folder_path = [root '1_TC/3_Model_Source/2_MegaWat/']; % Put here the path of where you downloaded the repository
forc_path = [root '7_Forcing/2_Extractions_Point/3_Radiation_Partition/']; % Put here the path of where you downloaded the repository


%% Study site details
%==========================================================================
%{
Sites:
   Cinca_Mido
   Tiber
   Velino
Outlet_names:
    "Cinca_Mid_OUT"
    "Apennine_out",
    "Monte_Terminillo"
    "Velino_OUT"
%}
%==========================================================================

SITE = ["Velino"]; 
FORCING = "ERA5Land";
UTM_zone = 33; % for Italy
DeltaGMT= 1; % for Italy

%% MODEL PARAMETERS
%==========================================================================
%
%==========================================================================

% Time step for the model
dt=3600; % [s]
dth=1; % [h]

% Integration interval for Solar variables
% Hours or fraction before and after. Values obtained from the
% Automatic_Radiation_Partition_I
t_bef = 1.5; t_aft = -0.5;


%% Pixel selection
%==========================================================================
% Depends on the point to be modelled
% Select for which pixel to run the point-scale version of T&C
% 1: AWS_OnGlacier
% 2: Pluvio
%==========================================================================
point_id = 3;     
LOC = point_id;

%% Data storing  
%==========================================================================
% How to store outputs
% Recommended for long-term simulations (> 10-20 years, otherwise data volume is too important)
%==========================================================================
output_daily = 1; 


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

%% Initial condition
% Saving initial conditions of the model
out = strcat(outlocation, '4_INIT_COND/INIT_COND_', SITE ,'_MultiPoint_',Point,'.mat'); % file path initial conditions

%% Load DEM and geographical information
%==========================================================================
% Different dtm
% "dtm_Tiber_250m.mat"
% "dtm_Cinca_Mid_250m.mat"
%==========================================================================
dtm_file = ["dtm_Velino_250m.mat"]; 
res = 250; % simulation resolution [m]
disp(strcat('Model resolution: ',num2str(res)));
dtm_file_op = strcat(folder_path,'5_Common_inputs/',SITE,'/',dtm_file);

load(dtm_file_op); % Distributed maps pre-processing. Useful here to get the DTM and initial snow depth
DTM = DTM_orig; % Use the full DEM in case running POI outside of mask
DTM(isnan(DTM)) == 0; %%??

[m_cell,n_cell]=size(DTM);

% Precipitation vertical gradient
% Pmod_S = MASK;
% rate = Pmod/(Z_max-Z_min);
% Pmod_S(DTM>Z_min) = 1+rate.*(DTM(DTM>Z_min)-Z_min);

%% Load CO2 data
load('Ca_Data.mat');
Ca_all = Ca;
topo = 1;

%% FAKE SNOW (TRIAL)
%load(strcat(folder_path,'4_Preparation_files/Apennine_preparation_Achille_Method_Distributed/','7_SnowCover/4_Snow_Depth/Apennine_Init_Snow_Depth_virtual.mat'));
%load(strcat(folder_path,'4_Preparation_files/Apennine_preparation_Achille_Method_Distributed/','7_SnowCover/5_Snow_Albedo/Apennine_Init_Snow_Albedo_virtual.mat'));

%% Impose measured albedo on glacier areas
fn_alb_elev = strcat(SITE, '_Albedo_vs_elev.mat');

if exist(fn_alb_elev,'file')>0
disp('Using measured glacier albedo');
load([SITE '_Albedo_vs_elev']);

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



%% FOR LOOP for locations
%Loc=1
Names = string(POI.Name);
loc = find(Point == Names);  

%% Definitions of site
%% Crowns
cc = POI.NCrown(loc); 
II = cell2mat(POI.II(loc));
Cwat = POI.Cwat(loc); 
Curb = POI.Curb(loc); 
Crock = POI.Crock(loc);
Cbare = POI.Cbare(loc);
Ccrown = cell2mat(POI.Ccrown(loc));
zatm = POI.zatm(loc);
id_location = string(POI.Name(loc));
cc_max = POI.cc_max(loc);

%% Locations
Zbas = POI.Zbas(loc);
Lat = POI.LAT(loc);
Lon = POI.LON(loc);
ij = POI.ij(loc);


%% FORCING
%==========================================================================
% ERA5
%==========================================================================

% Years needed
years_forc = year(date_start):year(date_end);
forc_opened = struct(); 

%% Extraction of forcing data
for yy = years_forc

forc_file= [forc_path 'Forcing_ERA5_Land_' char(SITE) '_' char(string(yy)) '_corr_all.mat']; % Put here the path of where you downloaded the repository
load(forc_file); % Load forcing table for the current POI

fieldName = ['year_', char(string(yy))];
forc_opened.(fieldName) = forc_f.(Point);

end

%% Combination of data into a single table
% First, extract all the tables from the struct into a cell array
forc_Array = struct2cell(forc_opened);

% Now, use vertcat to combine all the tables in the cell array
forc = vertcat(forc_Array{:});
%forc = forc_f.(Point);

Date_all=forc.Date; 

%define period and time zone info
x1=find(date_start == Date_all,1);
x2=find(date_end == Date_all,1);


%% Displaying modelling parameters
disp(strcat("Site selected: ", SITE));
disp(['Forcing selected: ' char(FORCING)]);
disp(['Running T&C for pixel: ' char(id_location)]);
disp(['Simulation period: ' datestr(date_start) ' to ' datestr(date_end)]);

%% Fetch time and do date handling
Date = Date_all(x1:x2);
[YE,MO,DA,HO,MI,SE] = datevec(Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;
clear YE MO DA HO MI SE

%% Carbon data
%==========================================================================
% Load carbon data and narrow down to period
%==========================================================================

formattedDate_CO2 = datetime(Date_CO2, 'ConvertFrom', 'datenum');

d1 = find(formattedDate_CO2 == Date(1));
d2 = find(formattedDate_CO2 == Date(end));

%d1 = find(abs(Date_CO2-datenum(Date(1)))<1/36);
%d2 = find(abs(Date_CO2-datenum(Date(end)))<1/36);
Ca=Ca_all(d1:d2); 
clear d1 d2

Oa= 210000; % Intercellular Partial Pressure Oxygen [umolO2/mol]

% Debugging
% formattedDate = datetime(Date_CO2, 'ConvertFrom', 'datenum');
% disp(strcat('Launcher',num2str(size(Ca,2))))


%% Meteorological input
%==========================================================================
%{
FORCINGS: 
    Precipitation [mm]
    Air Pressure [Pa]
    Temperature [C]
    Wind Speed [m s-1]    
    Radiation [W m-2]    
    Dew Point temperature [C]

CALCULATIONS:
    Relative Humidity [-]
%}
%==========================================================================

% Period of forcing data
forcing = forc(x1:x2,:);
NN= height(forcing);%%% time Step

% Height of virtual station
zatm_hourly = repmat(2.00,height(forcing),1);

%categories    [fir_high    larch     grass_A   grass_B    shrub    BLever_high   BLdec_high NoVeg]  
zatm_surface = [18          18        2         2          2        18            18         NaN  ];

zatm_hourly_on = 0;

%% Precipitation

% Precipitation from ERA5Land
Pr=double(forcing.Total_Precipitation_HH);
Pr(isnan(Pr))=0;
Pr(Pr<0.01)=0;

%Precipitation from stations


%% Air pressure
% Air Pressure in mbar based on the PDF for variables and parameters of TC
% Pressure comes in Pa from ERA5
Pre=double(forcing.Pressure/100);    

%% Temperature
% 2m air temperature
Ta=double(forcing.Temperature);

%% Wind Speed
Ws=double(forcing.Wind_Speed);
Ws(Ws < 0.01) = 0.01;

%% Relative humidity
% Divided by 100 to set the number in the range 0-1
U=double(forcing.RH/100);

%% Longwave radiation
% N can be cloud cover [-] or longwave incoming radiation [W m-2]
% Latm=forcing.LW_rad_downward_HH; % Latm:Incoming long wave radiation [W m-2]
% N=ones(NN,1); % cloud cover [-]
% N = forcing.LW_rad_downward_HH;
N = double(forcing.N);

%% Radiation partition
SAD1=double(forcing.SAD1); SAD2=double(forcing.SAD2); 
SAB1=double(forcing.SAB1); SAB2=double(forcing.SAB2);
PARB=double(forcing.PARB); PARD=double(forcing.PARD);

%% Albedo parameters
%Ameas = ones(NN,1);
alpha = 0; % switch for albedo
%Ameas_t=0; % albedo
%Aice_meas_on_hourly = ones(height(forcing),1)/2; % albedo
%Asno_meas_on_hourly = ones(height(forcing),1)/2; % albedo

%% Vapor pressure - Dew Point temperature
%esat/ea/Ds/Tdew
esat = double(forc.es);   % Vapour pressure at saturation (Pa)
ea = double(forc.ea);     % Vapour pressure (Pa)
Ds = esat - ea;  % Vapor Pressure Deficit (Pa)
Ds(Ds<0) = 0; 
Tdew = double(forc.Dew_Point_Temp);

%a=17.27; b=237.3;
%clear a b xr;
%xr=a*Ta./(b+Ta)+log10(U);
%Tdew=b*xr./(a-xr);          % Presumed dewpoint temperature (°C)

%% DING PARAMETERIZATION
% Initial daily mean values for Ding parametrization
Ta_Ding_d = nanmean(Ta(1:24));
Pre_Ding_d = nanmean(Pre(1:24));
ea_Ding_d = nanmean(ea(1:24));

%% RADIATION AND TOPOGRAPHY
num_cell=numel(DTM);

%MASK = MASK.*0+1;
MASKn=reshape(MASK,num_cell,1);

if topo == 1
    %load topography data and narrow down to period 
    Topo_data = load([folder_path,'4_Preparation_files/4_GeoTerrain_MultiPoint/4_Results/1_Velino/',char(SITE),'_ShF.mat']); % ShF matrix created during pre-processing step

    Par_time = Topo_data.Par_time;
    Par_points = Topo_data.Par_points;

    %% Topography
    %string(POI.Name(loc))
    %% Debugging
    Par_time_table=Par_time.(Point); %Table to use
    Par_time_period = Par_time_table(Par_time_table.Time>=date_start & Par_time_table.Time<=date_end, :); %variables per period

    %ShF
    ShF = Par_time_period.ShF_S;
    
    rho_g = 0.35; %%% Spatial Albedo
    zeta_S = Par_time_period.zeta_Sts;
    h_S = Par_time_period.h_Sts;
    %clear Par_time

    Par_points_data = Par_points(Par_points.Name == string(POI.Name(loc)),:);
    SvF = Par_points_data.SvF_S; % Sky view factor at pixel ij
    Ct = Par_points_data.Ct_S;
    Slo_top = Par_points_data.Slo_top_S; % Slope at pixel ij
    Aspect = Par_points_data.Aspect; % Aspect at pixel ij
    %clear Points

    cos_fst = cos(atan(Slo_top)).*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
    cos_fst(cos_fst<0)=0;

    %clear zeta_S 
    
    %% Radiative parameters
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
    SAB1(SAB1<0)=0; SAB2(SAB2<0)=0;
    PARB(PARB<0)=0; PARD(PARD<0)=0;
    SAD1(SAD1<0)=0; SAD2(SAD2<0)=0;

    SAB1(isnan(SAB1))= 0; SAB2(isnan(SAB2)) = 0;
    SAD1(isnan(SAD1))= 0; SAD2(isnan(SAD2)) = 0;
    PARB(isnan(PARB))= 0; PARD(isnan(PARD)) = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
   SvF=1; %% Sky View factor = 1 if topography is not considered
   Slo_top_S = MASKn*0;
   Slo_top=Slo_top_S(ij);
end

%% SOIL 
%==========================================================================
% Vector is divided by 100 to put the numbers within 0-1
% Original PSAN, PCLA and PORG come with values between 0-100, g/100g
%==========================================================================
%imagesc(PSAN)
%imagesc(PCLA)

PSAN=reshape(PSAN,num_cell,1)/100; Psan = PSAN(ij); % Soil sand content at pixel ij
PCLA=reshape(PCLA,num_cell,1)/100; Pcla = PCLA(ij); % Soil clay content at pixel ij
PORG=reshape(PORG,num_cell,1)/100; Porg= PORG(ij); % Soil organic content at pixel ij

%Psan = 0.244
%Pcla = 0.334
%Porg = 0.0608

ms=10 ; %% 11 ; 
SOIL_TH=reshape(SOIL_TH,num_cell,1);
ms_max = 10; %% Number of soil layers

%% GLACIERS
%%% DEBRIS
% INIT_COND_v2 has md_max parameter
% Do not delete this code
DEB_MAP=reshape(DEB_MAP,num_cell,1);
md_max = 10; %% % Number of debris layers

%% SNOW
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
%==========================================================================
% In MultiPoint analysis INIT_COND_Tiber_MultiPoint.mat is created.
% This can cause problems with the initial conditions. Specially with Ca.
%==========================================================================
% (run this only once in MultiPoint model!)
%if exist(out, 'file') == 2
%load(out);
%else
if topo == 1
INIT_COND_v3(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Ca,SNOWD,SNOWALB,out, II, ij);
load(out);
else 
INIT_COND_v2(num_cell,m_cell,n_cell,...
   cc_max,ms_max,md_max,...
   MASKn,GLH,Slo_top_S,ksv,Ca,SNOWD,SNOWALB,out);
load(out);
end
%end

%% Debugger
%disp(Nreserve)

%% RUN MODEL
%==========================================================================
% PARAM_IC: Define parameter file
% MAIN_FRAME: Contains the model
%==========================================================================
PARAM_IC = strcat(folder_path,'3_PointScale_version/3_Inputs/MOD_PARAM_Multipoint_Pro.m');
MAIN_FRAME_Pro; % Launch the main frame of T&C. Most of the things happen in this line of code

%% Post-compute calculations
%==========================================================================
% Negative evapotranspiration (ET) is dew
%==========================================================================
% post-compute sublimation from ESN
SSN = ESN.*(Ts<0);
ET(ET < 0) = 0; 

%% Output manager
%==========================================================================
%{
As outputs:
    Hourly outputs
    Daily outputs
    Parameters
%}
%==========================================================================

% Make date usable for R
Date_R = char(Date); 

%% Parameter output
Param_t = table(Lat,Lon,Zbas,dbThick,'VariableNames',{'Lat','Lon','Zbas','dbThick'});
Param_t = [Param_t, struct2table(SnowIce_Param), struct2table(Deb_Par)];
Param_t = rows2vars(Param_t);
Param_t = renamevars(Param_t,{'OriginalVariableNames','Var1'},{'Parameter','Value'});

% Exporting as .txt
writetable(Param_t, strcat(outlocation, '3_Params/',id_location,'_param.txt') )

% Here I manually choose the T&C outputs I want to save at each point.
Outputs_t = table(Date, ...
    EICE,ESN,SND,SWE,...
    Ta,Ws,U,N, ...
    SAD1+SAD2+SAB1+SAB2,Pre,Pr,Pr_sno, ...
    ALB,Smelt,Imelt,SSN, ...
    ICE,ET, ETen ,QE, ...
    ros, NDVI, T_H, T_L, ...
    O, FROCK, ...
'VariableNames',{'Date', ...
    'EICE','ESN','SND','SWE',...
    'Ta','Ws','U','N', ...
    'Rsw', 'Pre','Pr','Pr_sno', ...
    'Albedo','Smelt','Imelt','SSN', ...
    'ICE','ET','ETen','QE', ...
    'ros', 'NDVI', 'T_H', 'T_L', ...
    'O_mois','FROCK'});

%% Hourly output
%writetable(Outputs_t, strcat(outlocation, '1_Hourly/',id_location,'_hourly_results.txt'))

%% Daily outputs
% If daily outputs are activated

if output_daily == 1

Outputs_tt = table2timetable(Outputs_t); 

%% Mean of variables
Outputs_dm = retime(Outputs_tt,'daily',@nanmean); % For average variables
% Initial matrix set based on averages
Outputs_d = Outputs_dm; 

%% Some variables must be summed
Outputs_ds = retime(Outputs_tt,'daily',@nansum); % For sum variables

% Changing variables that must be summed in Outputs_d
Outputs_d.Pr = Outputs_ds.Pr;
Outputs_d.Pr_sno = Outputs_ds.Pr_sno;
Outputs_d.ET = Outputs_ds.ET;
Outputs_d.ETen = Outputs_ds.ETen;
Outputs_d.T_H = Outputs_ds.T_H;
Outputs_d.T_L = Outputs_ds.T_L;

%% Adding daily outputs
Outputs_d.LAI_H = sum(LAI_H.*Ccrown,2); % Leaf area index
Outputs_d.LAI_L = sum(LAI_L.*Ccrown,2);
Outputs_d.LAI = Outputs_d.LAI_H+Outputs_d.LAI_L;

Outputs_d.NPP_H = sum(NPP_H.*Ccrown,2); % Net primary production
Outputs_d.NPP_L = sum(NPP_L.*Ccrown,2);
Outputs_d.NPP = Outputs_d.NPP_H + Outputs_d.NPP_L;

Outputs_t = timetable2table(Outputs_d);

% Exporting as .txt
writetable(Outputs_t, strcat(outlocation, '2_Daily/',id_location,'_daily_results.txt'))

end 

%% Save of the environment
% Here only for one point - For the snow station. More points can be
% defined
if strcmp(POI{POI.Name == Point,'Feature'}{1}, 'Snow_station')
save(strcat(outlocation, '5_Env/Env_',id_location,'.mat'))
end

%% OTHER CALCULATIONS OUT OF THE MODEL
%==========================================================================
%
%==========================================================================
%{
% Resample the DEM to a new resolution (e.g., 30 meters)
DEM_resampled_out = resample(DEM, 500); 

% Display the original and resampled DEMs
imagesc(DEM);
figure;
imagesc(DEM_resampled);
%}

end