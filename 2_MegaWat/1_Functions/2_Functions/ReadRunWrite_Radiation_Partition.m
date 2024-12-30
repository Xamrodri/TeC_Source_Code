%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
CODE: RADIATION PARTITION 
AUTHOR: MAXIMILIANO RODRIGUEZ
INSTITUTION: ISTA
FROM: MIKE McCARTHY
DATE: DEC 30, 2024

PURPOSE: It uses Automatic_Radiation_Partition_ERA.m to perform the radiation
partition.

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CLEAN
close all
clear all

%% SITE AND DATE
%==========================================================================
% DEFINE SITE & PERIOD
%==========================================================================

% site name
site='CH';

% years of interest
yrs=2008:2008;

%% PATHS
%==========================================================================
% DEFINE PATHS
%==========================================================================
path_func='M:/19_ISTA/1_TC/3_Model_Source/2_MegaWat/1_Functions/2_Functions';
addpath(path_func)
path_forcings='M:/19_ISTA/18_Forcing/2_Apennines/';
% addpath(path)

%% Read initial input variables
forc = readtable(strcat(path_forcings,'10_Output/Forcing_ERA5_Land_2008_sta1.csv'));
GRAPH_VAR=0;

%% Date manipulation
date_old = num2str(forc.Date);
time_old = forc.Time;
forc.Date = datetime(str2num(date_old(:, 1:4)), str2num(date_old(:, 5:6)), str2num(date_old(:, 7:8)), time_old, 0, 0);
forc = removevars(forc, 'Time');

% prepare empty arrays
SAD1_df=NaN(size(forc,1),1);
SAD2_df=NaN(size(forc,1),1);
SAB1_df=NaN(size(forc,1),1);
SAB2_df=NaN(size(forc,1),1);
PARB_df=NaN(size(forc,1),1);
PARD_df=NaN(size(forc,1),1);
    
%use "parfor" instead?
%for j=1:n_sta

% Example
Date = forc.Date;
LAT=forc.Lat(1);
LON=forc.Lat(2);
ELE=500; %%%%%%%%%%%%%% CHECK
PR= forc.Total_Precipitation;
TDP= forc.Dew_Point_Temp;
SW=forc.SW_rad_downward; %%%%%%%%%%%%%%% CHECK WHY UNTIL HOUR 22, 23, 24
DeltaGMT = 1; %%%%%%%%%%%% CHECK
%% MAIN FUNCTION
%==================================================================
%{
"AUTOMATIC COMPUTATION OF RADIATION PARTITION.M"

Outputs:
    SAD1 [W m-2]: First band diffuse radiation
    SAD2 [W m-2]: Second band diffuse radiation 
    SAB1 [W m-2]: First band direct radiation 
    SAB2 [W m-2]: Second band direct radiation 
    PARB [W m-2]: PAR radiation direct 
    PARD [W m-2]: PAR radiation diffuse 
    t_bef [hours]: Optimized intergration interval for solar radiation variables 
    t_aft [hours or fraction after]: Optimized intergration interval for solar radiation variables
    N=Latm: cloudcover
%}
%==================================================================

% Achille function
[SAD1,SAD2,SAB1,SAB2,PARB,PARD]=Automatic_Radiation_Partition_ERA(Date,LAT,LON,ELE,DeltaGMT,PR,TDP,SW,GRAPH_VAR);
% Mike function
%[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Automatic_Radiation_Partition_ERA(Date,LAT,LON,ELE,DeltaGMT,PR,TDP,SW,GRAPH_VAR);

%% Function outputs
%==========================================================================
% WRITE FUNCTION OUTPUT TO ARRAY
%==========================================================================
        
forc.SAD1 = SAD1';
forc.SAD2 = SAD2';
forc.SAB1 = SAB1';
forc.SAB2 = SAB2';
forc.PARB = PARB';
forc.PARD = PARD';
        
%clearvars SD SB SAD1 SAD2 SAB1 SAB2 PARB PARD N Rsws t_bef t_aft
%clearvars LAT LON ELE PR TDP SW

%end

%==========================================================================
% WRITE FUNCTION OUTPUT TO FILE
%==========================================================================
fn=strcat(path_forcings,'10_Output/Forcing_ERA5_Land_2008_sta1_all.mat');

save(fn, 'forc');

%save(fn,'Date', 'SD_df','SB_df','SAD1_df','SAD2_df','SAB1_df','SAB2_df',...
%    'PARB_df','PARD_df'|,'N_df','Rsws_df','t_bef_df','t_aft_df');



