clc; clear all;

%Folder where files are stored from the distributed model
path_output = 'C:/Users/mxtvc/Desktop/TC/1_Output/ERA5Land_211024_2021_2021_Distributed_casestudy_test';
Directory=path_output;
cd(Directory)

load("Osat_1.mat")
load("OUTPUT_Kyzylsu_SNOWMAP_13.mat")

%Surface map of the catchment
figure(1)
surf(DTM)

%Spatial SWE
figure(2)
SWE_MAT = reshape(Smelt_spatial_daily, length(x_cell), length(y_cell))
surf(SWE_MAT)

%Spatial SWE
figure(3)
SSN_MAT = reshape(SSN_spatial_daily, length(x_cell), length(y_cell))
surf(SSN_MAT)

%Osat
load("Osat_1.mat")
figure(4)
surf(Osat)

%Osat
load("Ohy_1.mat")
figure(5)
surf(Ohy)

%Matrix
figure(6)
[M, c]=contour(DTM, [1000 1500 2000 2500 3000 3500 4000], "ShowText",true, ...
    "LabelFormat","%d m")
c.LineWidth = 2

%Contour and heatmap overlapped
figure(7)
Ohy_flipped = flipud(Ohy)
imagesc(Ohy_flipped);
caxis([min(Ohy_flipped(:)) max(Ohy_flipped(:))]);
colormap(sky);
colorbar;
set(gca, 'XTick', [], 'YTick', []);

hold on;
DTM_flipped = flipud(DTM)
[M2, c2]=contour(DTM_flipped, [1000 1500 2000 2500 3000 3500 4000 5000], "ShowText",true, ...
    "LabelFormat","%d m", 'LineColor', 'k');
hold off;

%Labels
xlabel({'EAST'})
ylabel({'NORTH'})

%xl = xlim; yl = ylim;

%X ticks
xcord_x = 1:length(x_cell);
xcord_y = 240*ones(1, length(x_cell));
xcord_text = strsplit(num2str(round(x_cell,0)));
xc_selected= 1:20:length(x_cell)

xticks(xc_selected);
xticklabels(xcord_text(xc_selected));

%Y ticks
new_y_cell = flip(y_cell)
ycord_x = 1:length(new_y_cell);
ycord_y = 240*ones(1, length(new_y_cell));
ycord_text = strsplit(num2str(round(new_y_cell,0)));
yc_selected= 1:20:length(new_y_cell)

yticks(yc_selected);
yticklabels(ycord_text(yc_selected));
%{
text(xcord_x(xc_selected), xcord_y(xc_selected),xcord_text(xc_selected) , ...
    'Rotation',45, 'VerticalAlignment','bottom', ...
    'HorizontalAlignment','center');
%}


switch_ice = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_ice == 1
    figure(17)
    set(gca,'FontSize',11);
    subplot(2,1,1);
    plot(NN,ICE_D,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Ice Depth');
    ylabel('[m]')
    subplot(2,1,2);
    plot(NN,ICE,'g','LineWidth', 1.5);
    hold on; grid on;
    title('ICE water');
    ylabel('[mm]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    figure(18)
    subplot(3,1,1);
    plot(NN,IP_wc,'g','LineWidth', 1.5);
    hold on; grid on;
    title('Water inside Icepack');
    ylabel('[mm]')
    set(gca,'FontSize',11);
    subplot(3,1,2);
    plot(NN,WR_IP,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Water released from Icepack');
    ylabel('[mm]')
    subplot(3,1,3);
    plot(NN,EICE,'g','LineWidth', 1.5);
    hold on; grid on;
    title('Ice Evaporation');
    ylabel('[mm/h]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

switch_snow=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_snow == 1
    figure(14)
    set(gca,'FontSize',11);
    subplot(3,1,1);
    plot(NN,SND,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Snow Depth');
    ylabel('[m]')
    subplot(3,1,2);
    plot(NN,SWE,'g','LineWidth', 1.5);
    hold on; grid on;
    title('SWE');
    ylabel('[mm]')
    subplot(3,1,3);
    plot(NN,ros,'g','LineWidth', 1.5);
    hold on; grid on;
    title('Snow Density');
    xlabel('Hour'); ylabel('[kg/m^3]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    figure(15)
    subplot(2,1,1);
    plot(NN,SP_wc,'g','LineWidth', 1.5);
    hold on; grid on;
    title('Water inside Snowpack');
    ylabel('[mm]')
    subplot(2,1,2);
    plot(NN,Qfm,'b','LineWidth', 1.5);
    hold on; grid on;
    title('Freezing/Melting Water in the Snowpack Heat');
    xlabel('Hour'); ylabel('[W/m^2]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    figure(16)
    set(gca,'FontSize',11);
    subplot(3,1,1);
    plot(NN,WR_SP,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Water from Snowpack');
    ylabel('[mm]')
    subplot(3,1,2);
    plot(NN,ESN,'g','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,ESN_In,'b','LineWidth', 1.5);
    title('Snow Evaporation');
    legend(' Snowpack Evap.','Intercepted Snow Evap.')
    ylabel('[mm/h]')
    subplot(3,1,3);
    plot(NN,U_SWE,'k','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,In_SWE,'r','LineWidth', 1.5);
    title('Intercepted Snow');
    xlabel('Hour'); ylabel('[mm]')
    legend(' Unload Snow','Intercepted Snow')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

switch_summary = 1;

if switch_summary == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear Pr_yr
    Yrs=year(Date);
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
    hold on; grid on;
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



