clc; clear all;

%Folder where files are stored from the distributed model
path_output = 'M:/19_ISTA/1_TC/3_Model_Source/3_Output_files/1_Model_outputs/';
result = 'Region_(Distributed_Ebro)_SubBasin_(Cinca_Mid)_Forcing_(Chelsa)_RunningDate_(251124_1200)_RunningPeriod_(2021_2021)/';
%Directory=path_output;
%cd(Directory)

%% Creating table tb from results
tb = readtable([path_output, result,'Cell_data_final/OUTPUT_Cinca_Mid_PIXEL_Cinca_Mid_OUT.dat']);
NN =  tb.Date;
                
switch_plots = [1 ... % ICE
                1 ... % SNOW
                1 ... % PRECIPITATION
                1 ... % RUNOFF
                1 ... % VEGETATION
                0 ... %
                ];

%% ICE PLOTS
if switch_plots(1) == 1
    
    %======================================================================
    % Ice Depth, ICE water
    %======================================================================
    figure(1)
    set(gca,'FontSize',11);
    subplot(2,1,1);
    plot(NN,tb.ICE_D,'r','LineWidth', 1.5);
    %hold on; grid on;
    title('Ice Depth');
    ylabel('[m]')
    subplot(2,1,2);
    plot(NN,tb.ICE,'g','LineWidth', 1.5);
    %hold on; grid on;
    title('ICE water');
    ylabel('[mm]')
    
    %======================================================================
    % Water inside Icepack, Water released from Icepack, Ice Evaporation
    %======================================================================
    figure(2)
    subplot(3,1,1);
    plot(NN,tb.IP_wc,'g','LineWidth', 1.5);
    %hold on; 
    grid on;
    title('Water inside Icepack');
    ylabel('[mm]')
    set(gca,'FontSize',11);
    subplot(3,1,2);
    plot(NN,tb.WR_IP,'r','LineWidth', 1.5);
    %hold on; 
    grid on;
    title('Water released from Icepack');
    ylabel('[mm]')
    subplot(3,1,3);
    plot(NN,tb.EICE,'g','LineWidth', 1.5);
    %hold on; 
    grid on;
    title('Ice Evaporation');
    ylabel('[mm/h]')

end

%% PLOTS FOR ATMOSPHERIC VARIABLES
if switch_plots(3) == 1

    %======================================================================
    % Snow Depth, Snow Density
    %======================================================================
    figure(3)
    set(gca,'FontSize',11);
    subplot(3,1,1);
    plot(NN,tb.Pre_S,'r','LineWidth', 1.5);
    hold on;
    grid on;
    title('Solid Precipitacion');
    ylabel('[m]')
    subplot(3,1,2);
    plot(NN,tb.SWE,'g','LineWidth', 1.5);
    hold on; grid on;
    title('SWE');
    ylabel('[mm]')    
    xlabel('Hour'); ylabel('[kg/m^3]')
    
    
end

%% PLOTS FOR RUNOFF
if switch_plots(4) == 1

    %======================================================================
    % Channel runoff, Runon
    %======================================================================
    figure(4)
    set(gca,'FontSize',11);
    subplot(2,1,1);
    plot(NN,tb.Q_channel,'r','LineWidth', 1.5);
    grid on;
    title('Channel discharge'); %%% Water in channels
    ylabel('[mm]')

    set(gca,'FontSize',11);
    subplot(2,1,2);
    plot(NN,tb.q_runon,'g','LineWidth', 1.5);
    grid on;
    title('Runon');
    ylabel('[mm/h]')    
    xlabel('Hour');
    
end

%% PLOTS FOR VEGETATION
if switch_plots(5) == 1

    %======================================================================
    % LAI
    %======================================================================
    figure(5)
    set(gca,'FontSize',11);
    subplot(2,1,1);
    plot(NN,tb.LAI_H,'r','LineWidth', 1.5);
    grid on;
    title('Leaf area index - High'); %%% Water in channels
    ylabel('[mm2 leaf area/m2 ground area]')

    set(gca,'FontSize',11);
    subplot(2,1,2);
    plot(NN,tb.LAI_L,'g','LineWidth', 1.5);
    grid on;
    title('Leaf area index - Low');
    ylabel('[mm2 leaf area/m2 ground area]')    
    xlabel('Hour');
    
end

%{
%% SUMMARY VARIABLES
switch_summary = 1;
if switch_summary == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Retrieving information from table
    Pr = tb.Pr_S;
    Date = tb.Date;

    Yrs=year(tb.Date);
    r=0;
    for i=min(Yrs):max(Yrs)
        if length(find(Yrs==i))>350*24
            r=r+1;
            Pr_yr(r)=sum(Pr(find(Yrs==i)));
            ET_yr(r)=sum(sum(T_H(find(Yrs==i),:),2)+sum(T_L(find(Yrs==i),:),2) + EG(find(Yrs==i)) +  ...
                sum(EIn_H(find(Yrs==i),:),2)+sum(EIn_L(find(Yrs==i),:),2)+EIn_urb(find(Yrs==i))+EIn_rock(find(Yrs==i)) + ...
                ESN(find(Yrs==i)) + ESN_In(find(Yrs==i)) );
            Lk_yr(r)=sum(Lk(find(Yrs==i)));
            T_yr(r)=sum(sum(T_H(find(Yrs==i),:),2)+sum(T_L(find(Yrs==i),:),2));
            VPD_yr(r)=mean(Ds(find(Yrs==i)));
        end
    end
    Yrs=year(Date(1):1:Date(end)+1);

    NPP_H = tb.NPP_H; 
    NPP_H = tb.NPP_H; 
    RA_H = tb.RA_H;
    RA_L = tb.RA_L;

    r=0;
    GPP_H=(NPP_H+RA_H);
    GPP_L=(NPP_L+RA_L);
    Yrs=Yrs(1:length(GPP_H));
    for i=min(Yrs):max(Yrs)
        if length(find(Yrs==i))>350
            r=r+1;
            GPP_yr(r)=  sum((GPP_H(find(Yrs==i),:)+GPP_L(find(Yrs==i),:))*Ccrown');
            NPP_yr(r)= sum((NPP_H(find(Yrs==i),:)+NPP_L(find(Yrs==i),:))*Ccrown');
            ANPP_yr(r)=  sum((ANPP_H(find(Yrs==i),:)+ANPP_L(find(Yrs==i),:))*Ccrown');
            Yrs_yr(r)= i;
        end
    end

    figure(20)
    subplot(1,2,1)
    plot(Pr_yr,ET_yr,'xk','LineWidth', 1.5);
    %hold on; grid on;
    xlabel('Pr [mm]'); ylabel('ET [mm]')
    title('Annual Value')
    subplot(1,2,2)
    plot(Pr_yr,Lk_yr,'xk','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('Recharge [mm]')
    title('Annual Value')
    figure(18)
    subplot(1,3,1)
    plot(Pr_yr,GPP_yr,'xm','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('GPP [gC/m^2]')
    title('Annual Value')
    subplot(1,3,2)
    plot(Pr_yr,NPP_yr,'xm','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('NPP [gC/m^2]')
    title('Annual Value')
    subplot(1,3,3)
    plot(Pr_yr,ANPP_yr,'xm','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('ANPP [gC/m^2]')
    title('Annual Value')
    figure(21)
    subplot(2,1,1)
    plot(Yrs_yr,GPP_yr,'-k','LineWidth', 1.5);
    hold on; grid on;
    plot(Yrs_yr,NPP_yr,'-r','LineWidth', 1.5);
    plot(Yrs_yr,ANPP_yr,'-g','LineWidth', 1.5);
    xlabel('Year'); ylabel('[gC/m^2]')
    title('Vegetation Productivities')
    legend('GPP','NPP','ANPP')
    subplot(2,1,2)
    plot(Yrs_yr,Pr_yr,'-k','LineWidth', 1.5);
    hold on; grid on;
    plot(Yrs_yr,ET_yr,'-g','LineWidth', 1.5);
    plot(Yrs_yr,T_yr,'-r','LineWidth', 1.5);
    plot(Yrs_yr,Lk_yr,'-m','LineWidth', 1.5);
    xlabel('Year'); ylabel('[mm]')
    title('Hydrologic Components')
    legend('Pr','ET','Transp.','Rech.')
    figure(22)
    %%%% Huang et al 2015 GCB
    subplot(3,2,1)
    plot(Yrs_yr,GPP_yr./ET_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC/m^2 mm]')
    title('EWUE = GPP/ET')
    subplot(3,2,2)
    plot(Yrs_yr,GPP_yr./T_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC/m^2 mm]')
    title('WUE_T = GPP/T')
    subplot(3,2,3)
    plot(Yrs_yr,0.001*GPP_yr.*VPD_yr./ET_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC kPa/m^2 mm]')
    title('IWUE = GPP*VPD/ET')
    subplot(3,2,4)
    plot(Yrs_yr,GPP_yr./T_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC/m^2 mm]')
    title('WUE_{leaf} = Ag/T')
    subplot(3,2,5)
    plot(Yrs_yr,ANPP_yr./Pr_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC/m^2 mm]')
    title('RAIN USE EFFICIENCY (Huxman et al., 2004)')
    subplot(3,2,6)
    plot(Yrs_yr,ET_yr./(ET_yr+Lk_yr),'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[-]')
    title('HORTON INDEX')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EG T_H T_L ESN ESN_In Rd Rh EIn_H+EIn_L+EIn_urb+EIn_rock Qsub Qi
    Nyr_plot = length(Date)/8766;
    clear Tm ESNm EGm Einm Qim Rhm Rdm Prm Lkm TT
    PRECIP = (Pr_liq+Pr_sno)*dth;
    TRASP = (T_H + T_L)*dth;
    EINT = (ELitter+EIn_H+EIn_L+EIn_urb+EIn_rock)*dth;
    QPER = sum(Qi_out-Qi_in,2) ;
    ESNOW = (EICE+ESN + ESN_In)*dth;
    %ET=ET*dth;
    for j=1:12
        %%%%%%%%%%%%%%%%%%%%%%%%%
        Tm(j)=sum(TRASP(Datam(:,2)==j));
        ESNm(j)=sum(ESNOW(Datam(:,2)==j));
        EGm(j)=sum(EG(Datam(:,2)==j));
        Einm(j)=sum(EINT(Datam(:,2)==j));
        Qim(j)=sum(QPER(Datam(:,2)==j));
        Rhm(j)=sum(Rh(Datam(:,2)==j));
        Rdm(j)=sum(Rd(Datam(:,2)==j));
        Prm(j)=sum(PRECIP(Datam(:,2)==j));
        Lkm(j)=sum(Lk(Datam(:,2)==j));
        %ETm(j)=sum(ET(Datam(:,2)==j));
        %Qsubm(j)= Prm(j)-Tm(j)-EGm(j) - ESNm(j)- Einm(j) -Qim(j) -Rhm(j) - Rdm(j); %%%
        TT(j)=Tm(j)+ESNm(j)+EGm(j)+Einm(j)+Qim(j)+Lkm(j)+Rhm(j)+Rdm(j);
        
        %%%%%%%%%%%%%%%%%%%%%%
    end
    figure(1001)
    set(gca,'FontSize',9);
    mmm=[0,1:12,13];
    fill(mmm,[0 Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm+EGm 0]/Nyr_plot,'y')
    hold on ;
    fill(mmm,[0 Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm 0]/Nyr_plot,'r')
    fill(mmm,[0 Rhm+Rdm+ESNm+Einm+Lkm+Qim 0]/Nyr_plot,'g')
    fill(mmm,[0 Rhm+Rdm+ESNm+Einm 0]/Nyr_plot,'c')
    fill(mmm,[0 Rhm+Rdm 0]/Nyr_plot,'b')
    plot(1:12,Prm/Nyr_plot,'o--b','LineWidth', 1.5);
    %plot(1:12,ETm/Nyr_plot,'o--r','LineWidth', 1.5);
    grid on;
    xlim([1 12])
    %axis([1 12 0 max(Prm/Nyr_plot)+12])
    legend('Soil Evap.','Transp.','Rec. + Lat. Flow','Interc. Evap.','Runoff','Precip.')
    ylabel('[mm]'); xlabel('Month')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1002)
    set(gca,'FontSize',9);
    mmm=[0,1:12,13];
    fill(mmm,[0 Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm+EGm 0]/Nyr_plot,'y')
    hold on ;
    fill(mmm,[0 Rhm+Rdm+Lkm+Qim 0]/Nyr_plot,'r')
    fill(mmm,[0 Rhm+Rdm+Qim 0]/Nyr_plot,'b')
    fill(mmm,[0 Rhm+Rdm 0]/Nyr_plot,'g')
    fill(mmm,[0 Rhm 0]/Nyr_plot,'c')
    plot(1:12,Prm/Nyr_plot,'o--b','LineWidth', 1.5);
    grid on;
    xlim([1 12])
    %axis([1 12 0 max(Prm/Nyr_plot)+12])
    legend('Evapotransp.','Rec.','Lat. flow','Sat. Excess','Infilt. Excess','Precip.')
    ylabel('[mm]'); xlabel('Month')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1003)
    set(gca,'FontSize',9);
    mmm=[0,1:12,13];
    fill(mmm,[0 (Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm+EGm)./TT 0],'y')
    hold on ;
    fill(mmm,[0 (Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm)./TT 0],'r')
    fill(mmm,[0 (Rhm+Rdm+ESNm+Einm+Lkm+Qim)./TT 0],'g')
    fill(mmm,[0 (Rhm+Rdm+ESNm+Einm)./TT 0],'c')
    fill(mmm,[0 (Rhm+Rdm)./TT 0],'b')
    grid on;
    axis([1 12 0 1])
    legend('Soil Evap.','Transp.','Rec. + Lat. Flow','Interc. Evap.','Runoff')
    ylabel('Fraction'); xlabel('Month')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1004)
    set(gca,'FontSize',9);
    mmm=[0,1:12,13];
    fill(mmm,[0 (Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm+EGm)./TT 0],'y')
    hold on ;
    fill(mmm,[0 (Rhm+Rdm+Lkm+Qim)./TT 0],'r')
    fill(mmm,[0 (Rhm+Rdm+Qim)./TT 0],'b')
    fill(mmm,[0 (Rhm+Rdm)./TT 0],'g')
    fill(mmm,[0 (Rhm)./TT 0],'c')
    grid on;
    axis([1 12 0 1])
    legend('Evapotransp.','Rec.','Lat. flow','Sat. Excess','Infilt. Excess')
    ylabel('Fraction'); xlabel('Month')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%}

disp("Plots done")