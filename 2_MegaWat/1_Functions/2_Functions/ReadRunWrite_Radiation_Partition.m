close all
clear all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE SITE & PERIOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%site name
site='CH';

%%%years of interest
yrs=2015:2022;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_func='C:\Users\Buri\switchdrive\WSL\T&C\Code\Functions';
addpath(path_func)
path='E:\Meteodata\CH\RadiationPartition';
% addpath(path)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  READ INPUT FOR FUNCTION PER YEAR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_yrs=length(yrs); %number of years to partition
for i=1:n_yrs
    YR=yrs(i);

    %%%loading input variables
    load(strcat(path,'\Automatic_Radiation_Partition_I_INPUT___',site,'_',num2str(YR),'.mat'))
    GRAPH_VAR=0;

    %%%extract dimensions
    sta_nms=fieldnames(Pr);
    n_sta=length(sta_nms); %number of stations to partition
    n_ts=length(Date);

    %replace german umlauts & french accents
    sta_nms=rep_specialcharacters(sta_nms);
    nms=rep_specialcharacters(fieldnames(Pr)); Pr=RenameField(Pr,fieldnames(Pr),nms); clear nms
    nms=rep_specialcharacters(fieldnames(Rsw)); Rsw=RenameField(Rsw,fieldnames(Rsw),nms); clear nms
    nms=rep_specialcharacters(fieldnames(Tdew)); Tdew=RenameField(Tdew,fieldnames(Tdew),nms); clear nms
    
    superdisp('+++ Applying "Radiation Partition" to ',n_sta,' stations in ',YR,' (year ',i,' of ',n_yrs,') +++')

    %%%prepare empty arrays
    SD_df=NaN(n_ts,n_sta);
    SB_df=SD_df;
    SAD1_df=SD_df;
    SAD2_df=SD_df;
    SAB1_df=SD_df;
    SAB2_df=SD_df;
    PARB_df=SD_df;
    PARD_df=SD_df;
    N_df=SD_df;
    Rsws_df=SD_df;
    t_bef_df=NaN(1,n_sta);
    t_aft_df=t_bef_df;
    %use "parfor" instead?
    for j=1:n_sta
        nm=sta_nms{j};

        LAT=Lat(j);
        LON=Lon(j);
        ELE=double(Zbas(j));
        PR=Pr.(nm);
        TDP=Tdew.(nm);
        SW=double(Rsw.(nm));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % RUN FUNCTION
        % "AUTOMATIC COMPUTATION OF RADIATION PARTITION.M"
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Automatic_Radiation_Partition_I(Date,LAT,LON,ELE,DeltaGMT,PR,TDP,SW,GRAPH_VAR);
        %This outputs:
        %SAD1 %First band diffuse radiation (W m-2)
        %SAD2 %Second band diffuse radiation (W m-2)
        %SAB1 %First band direct radiation (W m-2)
        %SAB2 %Second band direct radiation (W m-2)
        %PARB %PAR radiation direct (W m-2)
        %PARD %PAR radiation diffuse (W m-2)
        %t_bef %Optimized intergration interval for solar radiation variables (hours)
        %t_aft %Optimized intergration interval for solar radiation variables (hours or fraction after)
        %N=Latm; cloudcover

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % WRITE FUNCTION OUTPUT TO ARRAY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SD_df(:,j)=SD;
        SB_df(:,j)=SB;
        SAD1_df(:,j)=SAD1;
        SAD2_df(:,j)=SAD2;
        SAB1_df(:,j)=SAB1;
        SAB2_df(:,j)=SAB2;
        PARB_df(:,j)=PARB;
        PARD_df(:,j)=PARD;
        N_df(:,j)=N;
        Rsws_df(:,j)=Rsws;
        t_bef_df(:,j)=t_bef;
        t_aft_df(:,j)=t_aft;

        clearvars SD SB SAD1 SAD2 SAB1 SAB2 PARB PARD N Rsws t_bef t_aft
        clearvars LAT LON ELE PR TDP SW

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % WRITE FUNCTION OUTPUT TO FILE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fn=strcat(path,'\Automatic_Radiation_Partition_I_OUTPUT___',site,'_',num2str(YR),'.mat');
    save(fn,'Date', 'SD_df','SB_df','SAD1_df','SAD2_df','SAB1_df','SAB2_df',...
        'PARB_df','PARD_df','N_df','Rsws_df','t_bef_df','t_aft_df');

    clearvars Date DeltaGMT Lat Lon Zbas Pr Tdew Rsw sta_nms n_sta n_ts ...
        SD_df SB_df SAD1_df SAD2_df SAB1_df SAB2_df ...
        PARB_df PARD_df N_df Rsws_df t_bef_df t_aft_df
    superdisp('+++ ',fn,' +++')
    superdisp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
end
