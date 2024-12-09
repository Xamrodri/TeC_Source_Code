%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation initiated on Dec 03, 2024
% Author: Maximiliano Rodriguez
% Code originally from: Achille Jouberton
% Area of Study: Pyrenees
% Region: Cinca Subcatchment
% Code explanation: This codes prepare the matrices for forcing for the
% Multi Point model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clear all;
close all;
clc;
clear;

%% INITIAL PARAMETERS FOR MODEL INITIALIZATION 
Bdrive = 0; 
ISTA = 1; % 1:local computer, 2: Cluster
s = 2; % catchment selection
%% DIRECTORIES

if ISTA == 1
root = 'M:/19_ISTA/1_TC/3_Model_Source/2_MegaWat/';
elseif ISTA == 2 %%% << HYPERION >>
root = '/nfs/scistore18/pelligrp/mrodrigu/1_TC/2_MegaWat';
else %%% << Incorrect input >>
error("Incorrect parameter for directory")
end

%% Other variables
if ~exist('fnum ','var')
     fnum = 1;  % If not be defined to the bash script
end

SITEs = ["Langtang", "Cinca_Mid"];
SITE = char(SITEs(s));

FORCINGs = ["ERA5Land","NHM-200m"]; 
FORCING = char(FORCINGs(fnum));

%% Period of modelling
x1s = ["01-Oct-2016 00:00:00","01-Sep-2022 00:00:00"];
x2s = ["30-Sep-2019 23:00:00","30-Sep-2023 23:00:00"];

x1 = char(x1s(s));
x2 = char(x2s(s));


%'FF'; 'LWIN'; 'PP'; 'RH'; 'TA'; 'SWPART'; 'PRESS'
vars_DSs = [
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,0,0,0,0,0;
    0,0,1,0,0,0,0;
    0,0,1,0,0,0,0;
    0,0,0,0,0,0,0];

vars_DS = vars_DSs(s,:);

sim_nm = [SITE '_' FORCING '_' x1(8:11) '_' x2(8:11)];

Lats = [28.2108, 39.0969];
Lons = [85.5695, 71.4176];
DeltaGMTs=[5.75,5,8,5.45,8,8];
Zbass=[3862,3579,3800,5449,4600,5850]; %in this setup: elevation of the reference pressure logger

if strcmp('ERA5Land',FORCING)
    t_befs = [1,0.5,0.75,0.75,1,1];
    t_afts = [0,0.5,0.25,0.25,0,0];
elseif strcmp('NHM-200m',FORCING)
    t_befs = [1,0.5,0.75,1.5];
    t_afts = [0,0.5,0.25,-0.5];
end

%% DTM
dtm_files = ["dtm_Langtang_100m.mat", "dtm_Cinca_Mid_250m.mat"];

Lat = Lats(s); Lon = Lons(s);
DeltaGMT=DeltaGMTs(s); t_bef=t_befs(s); t_aft=t_afts(s); 
dtm_file = char(dtm_files(s)); Zbas=Zbass(s); 
res = str2num(extractBetween(string(dtm_file), [SITE '_'], 'm.mat')); % model resolution

%% Path to distributed forcing
addpath(genpath([root '5_Common_inputs']));
addpath(genpath([root '1_Functions']));
outlocation = [root '3_Pyrenees_PointScale/2_Forcing/'];

load(dtm_file)

%% EXTRACT TOPOGRAPHIC/SHADING DATA FROM DISTRIBUTED set-up
%==========================================================================
% 
%==========================================================================

[Slo_top,Aspect]=Slope_Aspect_indexes(DTM_orig,cellsize,'mste');
Aspect(isnan(Aspect))=0;
Slo_top(Slo_top<0.001)=0.001;
[HZ,Zasp] = Horizon_Angle(DTM_orig,cellsize);
[SvF,Ct] = Sky_View_Factor(DTM_orig,atan(Slo_top)*180/pi,Aspect,HZ,Zasp);

x1 = datetime(x1);
x2 = datetime(x2);
Top_Date = (x1:hours(1):x2)';
numcell = numel(DTM);

[YE,MO,DA,HO,MI,SE] = datevec(Top_Date);
Datam(:,1) = YE; Datam(:,2)= MO; Datam(:,3)= DA; Datam(:,4)= HO;

%if from Achille code.  
%if (exist('poi_id','var')&& poi_id==1) %||  ~exist('poi_id','var') %Opt.1 if run on a cluster as an array task, use cluster array task ID as loc

[Slo_top_S,Aspect_S]=Slope_Aspect_indexes(DTM_orig,cellsize,'mste');
Aspect_S(isnan(Aspect_S))=0;
Slo_top_S(Slo_top_S<0.001)=0.001;
[HZ,Zasp] = Horizon_Angle(DTM_orig,cellsize);
[SvF_S,Ct_S] = Sky_View_Factor(DTM_orig,atan(Slo_top)*180/pi,Aspect_S,HZ,Zasp);

ShF_S= zeros(numcell,length(Top_Date));
h_Sts = zeros(length(Top_Date),1);
zeta_Sts = zeros(length(Top_Date),1);

parfor t = 1:length(Top_Date)
    t
    Datam_S=Datam(t,:);
    [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day]= SetSunVariables(Datam_S,DeltaGMT,Lon,Lat,t_bef,t_aft);
    h_Sts(t)=h_S;
    zeta_Sts(t) = zeta_S;
    ShF_S(:,t) = reshape(Shadow_Effect(DTM,h_S,zeta_S,HZ,Zasp),numcell,1);
end
Slo_top_S=reshape(Slo_top_S,numcell,1);
Aspect_S = reshape(Aspect_S,numcell,1);
SvF_S = reshape(SvF_S, numcell,1);

save([outlocation SITE '_ShF_' FORCING '.mat'], 'Top_Date','ShF_S','Slo_top_S','SvF_S','Ct_S','h_Sts','Aspect_S','zeta_Sts','cellsize', '-v7.3')

%end

%==========================================================================
% EXTRACT FORCING DATA FROM DOWNSCALED/BIASCORRECTED ERA5-Land or NHM-200m
%==========================================================================

%{
% Read POI table
tic
T=readtable([path_poi '/' SITE '_MultiPoints.txt']);

T.Name_orig = T.Name;
T.Name = strrep(T.Name,'.','_');

%define which variables to load as monthly spatial input     
vars_sto = string({'FF'; 'LWIN'; 'PP'; 'RH'; 'TA'; 'PRESS';'SAD1';'SAD2';'SAB1';'SAB2';'PARD';'PARB' });

DateHR = x1:hours(1):x2;
Date = x1:calmonths(1):x2;
YR = year(Date);
MO = month(Date);

%load forcing_POIs.mat
%Create basic strucutre to load data into

tab  = array2table(NaN(length(DateHR),length(vars_sto)), 'VariableNames', vars_sto);
tab_s = table2timetable(tab,'RowTimes',DateHR);

% if exist('poi_id','var') %Opt.1 if run on a cluster as an array task, use cluster array task ID as loc
%     xx_all = poi_id;
% else
xx_all = 1:size(T,1);
% end 

for xx=xx_all
    BASE.(T.Name{xx}) = tab_s;
end

%%
if strcmp(FORCING,'ERA5Land')

vars = string({'FF'; 'LWIN'; 'PP'; 'RH'; 'TA'; 'SWPART'; 'PRESS'});
if strcmp(SITE,"Rolwaling"); SITE = "Trambau"; end

for m = 1:length(Date)
%     tic
    yr = YR(m);
    mo = MO(m);
    idx = DateHR.Month == mo & DateHR.Year == yr;
        for v = 1:length(vars)
            if vars_DS(v)==1 && v ~= 6
                fn = char(strcat('Downscaled_', vars(v,:), '_',upper(SITE),'_ERA5Land_YY_', string(yr), '_MM_', string(mo), '.mat'));
            elseif vars_DS(v)==1 && v == 6
                fn = char(strcat('Downscaled_', vars(v,:), '_',upper(SITE),'_ERA5Land_', string(yr), '_', string(mo), '.mat'));
            else
                fn = char(strcat('BiasCorrected_', vars(v,:), '_',upper(SITE),'_ERA5Land_', string(yr), '_', string(mo), '.mat'));
            end
            load(fn)   
        end 

        if exist('FFds',"var");FF_BC = FFds; clear FFds; end
        if exist('LWINds',"var");LWIN_BC = LWINds; clear LWINds; end
        if exist('PPds',"var");PP_BC = PPds; clear PPds; end
        if exist('RHds',"var");RH_BC = RHds; clear RHds; end
        if exist('TAds',"var");TA_BC = TAds; clear TAds;end
        if exist('SWPARTds',"var");SWPART_BC = SWPARTds; clear SWPARTds; end
        if exist('PRESSds',"var")
          if strcmp(SITE,"Rolwaling")
              PRESS_BC = PRESSds.*0.01; clear PRESSds; % To remove once the pressure unit problem is solved from Trambau
          else 
              PRESS_BC = PRESSds; clear PRESSds; 
          end
        end 

        FF_BC = flipud(FF_BC);
        LWIN_BC = flipud(LWIN_BC);
        RH_BC = flipud(RH_BC);
        TA_BC = flipud(TA_BC);
        PRESS_BC = flipud(PRESS_BC);
        PP_BC = flipud(PP_BC);
        SAD1 = flipud(SAD1);
        SAD2 = flipud(SAD2);
        SAB1 = flipud(SAB1);
        SAB2 = flipud(SAB2);
        PARB = flipud(PARB);
        PARD = flipud(PARD); 
        
  for s = xx_all
            row = T.ROW(s);
            col = T.COL(s);
            BASE.(T.Name{s}).FF(idx) = squeeze(FF_BC(row,col,:));
            BASE.(T.Name{s}).LWIN(idx) = squeeze(LWIN_BC(row,col,:));
            BASE.(T.Name{s}).RH(idx) = squeeze(RH_BC(row,col,:));
            BASE.(T.Name{s}).TA(idx) = squeeze(TA_BC(row,col,:));
            BASE.(T.Name{s}).PRESS(idx) = squeeze(PRESS_BC(row,col,:));
            BASE.(T.Name{s}).PP(idx) = squeeze(PP_BC(row,col,:));
            BASE.(T.Name{s}).SAD1(idx) = squeeze(SAD1(row,col,:));
            BASE.(T.Name{s}).SAD2(idx) = squeeze(SAD2(row,col,:));
            BASE.(T.Name{s}).SAB1(idx) = squeeze(SAB1(row,col,:));
            BASE.(T.Name{s}).SAB2(idx) = squeeze(SAB2(row,col,:));
            BASE.(T.Name{s}).PARB(idx) = squeeze(PARB(row,col,:));
            BASE.(T.Name{s}).PARD(idx) = squeeze(PARD(row,col,:));
            
  end
%   toc
end

%store into individual files for each point, round to 4 digits, to save space
%saving yearly by overwriting old files, in case it crashes somewhere through the process
    for s = xx_all 
        BASE1 = BASE;
        forcing_all = BASE.(T.Name{s});
        forcing_all.Variables = round(forcing_all.Variables,4);
        save(strcat(outlocation,T.Name_orig{s},'.mat'),'forcing_all')
    end
end
%%
if strcmp(FORCING,'NHM-200m')

        vars = string({'FF'; 'LWIN'; 'PP'; 'RH'; 'TA'; 'PRESS'; 'SWPART'});
        vars_nhm = ["WS","LWin","Ptot","RH","Ta","Pres","SWPART"];
 
      for v = 1:length(vars)-1
             load(strcat("NHM_", vars(v), "_200m_2018_2019"),vars_nhm(v))
      end 
             load('NHM_SWPART_200m_2018_2019')


        FF_BC = flipud(WS);
        LWIN_BC = flipud(LWin);
        RH_BC = flipud(RH);
        TA_BC = flipud(Ta);
        PRESS_BC = flipud(Pres.*0.01);
        PP_BC = flipud(Ptot);
        SAD1 = flipud(SAD1);
        SAD2 = flipud(SAD2);
        SAB1 = flipud(SAB1);
        SAB2 = flipud(SAB2);
        PARB = flipud(PARB);
        PARD = flipud(PARD);

        for s = 1:size(T,1)
            row = T.ROW(s);
            col = T.COL(s);
            BASE.(T.Name{s}).FF = squeeze(FF_BC(row,col,:));
            BASE.(T.Name{s}).LWIN = squeeze(LWIN_BC(row,col,:));
            BASE.(T.Name{s}).RH = squeeze(RH_BC(row,col,:));
            BASE.(T.Name{s}).TA = squeeze(TA_BC(row,col,:));
            BASE.(T.Name{s}).PRESS = squeeze(PRESS_BC(row,col,:));
            BASE.(T.Name{s}).PP = squeeze(PP_BC(row,col,:));
            BASE.(T.Name{s}).SAD1 = squeeze(SAD1(row,col,:));
            BASE.(T.Name{s}).SAD2 = squeeze(SAD2(row,col,:));
            BASE.(T.Name{s}).SAB1 = squeeze(SAB1(row,col,:));
            BASE.(T.Name{s}).SAB2 = squeeze(SAB2(row,col,:));
            BASE.(T.Name{s}).PARB = squeeze(PARB(row,col,:));
            BASE.(T.Name{s}).PARD = squeeze(PARD(row,col,:));
            
        end

    for s =xx_all
      BASE1 = BASE;
      forcing_all = BASE.(T.Name{s});
      forcing_all.Variables = round(forcing_all.Variables,4);
      save(strcat(outlocation,T.Name_orig{s},'.mat'),'forcing_all')
    end

end
%%    
%store into individual files for each point, round to 4 digits, to save space
%saving yearly by overwriting old files, in case it crashes somewhere through the process
toc

% p = profile('info');
% % save([outlocation 'profile_outputs'],'p')
% profsave(p, [outlocation 'profile_outputs'])

%}