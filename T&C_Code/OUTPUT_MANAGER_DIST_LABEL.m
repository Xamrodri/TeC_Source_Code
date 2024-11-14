%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% OUTPUT WRITING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ADAPTED: PASCAL BURI, 28 FEBRUARY 2022
% Modified by Achille Jouberton since 2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% MASS BALANCE VARIABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    Qlat_in_tgtm1 = 0;
    q_runon_tgtm1 = 0;
    Q_channel_tgtm1 = 0;
else
    V_tgtm1 = V_tg;
    Vice_tgtm1 = Vice_tg;
    SWE_tgtm1 = SWE_tg;
    In_tgtm1 =In_tg;
    ICE_tgtm1 =  ICE_tg; %%
    WAT_tgtm1 =  WAT_tg; %%
    FROCK_tgtm1 = FROCK_tg; %%
    %%%%---
    Qlat_in_tgtm1 = Qlat_in_tg;
    q_runon_tgtm1 = q_runon_tg;
    Q_channel_tgtm1 = Q_channel_tg;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% SET NON-SOIL PIXELS IN 
%%%%%%%  SOIL MOISTURE VARIABLES TO ZERO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O(ksv>4,:)=0;
OF(ksv>4,:)=0;
OS(ksv>4,:)=0;
V(ksv>4,:)=0;
V_ice(ksv>4,:)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL AVERAGE OVER THE WATERSHED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alp_soil_tg	=	mean(alp_soil(Kinde));	%%[-]
Ca_tg       =	mean(Ca_S(Kinde));	%%[ppm]
Cice_tg     =	mean(Cice(Kinde));	%%[]
Cicew_tg	=	mean(Cicew(Kinde));	%%[]
CK1_tg      =	mean(CK1(Kinde));	%%[mm]
Csno_tg     =	mean(Csno(Kinde));	%%[]
Csnow_tg	=	mean(Csnow(Kinde));	%%[]
DQ_S_tg     =	mean(DQ_S(Kinde));	%%[]
dQ_S_tg     =	mean(dQ_S(Kinde));	%%[]
Dr_H_space	=	sum(Dr_H,2);	
Dr_H_tg     =	mean(Dr_H_space(Kinde));	%%[mm]
Dr_L_space	=	sum(Dr_L,2);	
Dr_L_tg     =	mean(Dr_L_space(Kinde));	%%[mm]
Ds_tg       =	mean(Ds_S(Kinde));	%%[°C]
DT_S_tg     =	mean(DT_S(Kinde));	%%[]
dw_SNO_tg	=	mean(dw_SNO(Kinde));	%%[]
ea_tg       =	mean(ea_S(Kinde));	%%[Pa]
EG_tg       =	mean(EG(Kinde));	%%[mm/h]
EICE_tg     =	mean(EICE(Kinde));	%%[mm/h]
EIn_H_space	=	sum(EIn_H,2);	
EIn_H_tg	=	mean(EIn_H_space(Kinde));	%%[mm/h]
EIn_L_space	=	sum(EIn_L,2);	
EIn_L_tg	=	mean(EIn_L_space(Kinde));	%%[mm/h]
EIn_rock_tg	=	mean(EIn_rock(Kinde));	%%[mm/h]
EIn_urb_tg	=	mean(EIn_urb(Kinde));	%%[mm/h]
er_tg       =	mean(er(Kinde));	%%[kg/s m^2]
ESN_tg      =	mean(ESN(Kinde)+ESN_In(Kinde));	%%[mm/h]
SSN_tg      =	mean(SSN(Kinde)+SSN_In(Kinde));	%%[mm/h]
EWAT_tg     =	mean(EWAT(Kinde));	%%[mm/h]
f_tg        =	mean(f(Kinde)*dth);	%%[mm]
Fract_sat_tg =	sum(ZWT(Kinde)==0)/num_cell;	%%[-]
FROCK_tg	=	mean(FROCK(Kinde));	%%[mm]
G_tg        =	mean(G(Kinde));	%%[W/m^2]
Gfin_tg     =	mean(Gfin(Kinde));	%%[W/m^2]
H_tg        =	mean(H(Kinde));	%%[W/m^2]
ICE_D_tg	=	mean(ICE_D(Kinde));	%%[m]
ICE_tg      =	mean(ICE(Kinde));	%%[mm]
Imelt_tg	=	mean(Imelt(Kinde));	%%[mm]
In_H_space	=	sum(In_H,2);	
In_L_space	=	sum(In_L,2);	
In_rock_tg	=	mean(In_rock(Kinde));	%%[mm]
In_SWE_tg	=	mean(In_SWE(Kinde));	%%[mm]
In_urb_tg	=	mean(In_urb(Kinde));	%%[mm]	
IP_wc_tg	=	mean(IP_wc(Kinde));	%%[mm]
Lk_rock_tg	=	mean(Lk_rock(Kinde)*dth);	%%[mm]
Lk_tg       =	mean(Lk(Kinde)*dth);	%%[mm]
Lk_wat_tg	=	mean(Lk_wat(Kinde)*dth);	%%[mm]
N_tg        =	mean(N_S(Kinde));	%%[-]
NDVI_tg     =   mean(NDVI(Kinde));
NIce_tg     =	mean(NIce(Kinde));	%%[mm]
NIn_SWE_tg	=	mean(NIn_SWE(Kinde));	%%[mm]
OF_tg       =	mean(OF(Kinde));	%%[]
OS_tg       =	mean(OS(Kinde));	%%[]	
Pr_liq_tg	=	mean(Pr_liq(Kinde));	%%[mm]
Pr_sno_tg	=	mean(Pr_sno(Kinde));	%%[mm]
Pr_tg       =	mean(Pr_S(Kinde));	%%[mm]
Pre_tg      =	mean(Pre_S(Kinde));	%%[mbar]
Q_channel_tg =	mean(Q_channel(Kinde));	%%[mm]
q_runon_tg	=	mean(q_runon(Kinde)*dth);	%%[mm]
QE_tg       =	mean(QE(Kinde));	%%[W/m^2]
Qfm_tg      =	mean(Qfm(Kinde));	%%[W/m^2]
Qi_in_space	=	sum(Qi_in*dth,2);	
Qi_out_space =	sum(Qi_out*dth,2);	
Qlat_in_tg	=	mean(Qi_in_space(Kinde));	%%[mm]
Qlat_out_tg	=	mean(Qi_out_space(Kinde));	%%[mm]
Qv_tg       =	mean(Qv(Kinde));	%%[W/m^2]
r_soil_tg	=	mean(r_soil(Kinde));	%%[s/m]
ra_tg       =	mean(ra(Kinde));	%%[s/m]
Rd_tg       =	mean(Rd(Kinde));	%%[mm]
Rh_tg       =	mean(Rh(Kinde));	%%[mm]
Rn_tg       =	mean(Rn(Kinde));	%%[W/m^2]  %%% NET RADIATION [W/m^2] = Rnet_ground + Rnet_vegH + Rnet_vegL + Rnet_snow + Rnet_rock +Rnet_wat + Rnet_urb + Rnet_ice + Rnet_deb;
ros_tg      =	mean(ros(Kinde));	%%[kg/m^3]	
SE_rock_tg	=	mean(SE_rock(Kinde));	%%[]
SE_urb_tg	=	mean(SE_urb(Kinde));	%%[]
Smelt_tg	=	mean(Smelt(Kinde));	%%[mm]
SND_tg      =	mean(SND(Kinde));	%%[m]
snow_albedo_tg = mean(snow_albedo(Kinde)); %%[-]
SP_wc_tg	=	mean(SP_wc(Kinde));	%%[mm]
SWE_avalanched_tg =   mean(SWE_avalanched(Kinde)); %%[mm]
SWE_tg      =	mean(SWE(Kinde));	%%[mm]
T_H_space	=	sum(T_H,2);	
T_H_tg      =	mean(T_H_space(Kinde));	%%[mm/h]
T_L_space	=	sum(T_L,2);	
T_L_tg      =	mean(T_L_space(Kinde));	%%[mm/h]
Ta_tg       =	mean(Ta_S(Kinde));	%%[°C]
Tdamp_tg	=	mean(Tdamp(Kinde));	%%[°C]
Tdew_tg     =	mean(Tdew_S(Kinde));	%%[°C]
Tdp_space	=	mean(Tdp,2);	
Tdp_tg      =	mean(Tdp_space(Kinde));	%%[°C]
Tice_tg     =	mean(Tice(Kinde));	%%[C]
Ts_tg       =	mean(Ts(Kinde));	%%[°C]
TsVEG_tg	=	mean(TsVEG(Kinde));	%%[°C]
U_SWE_tg	=	mean(U_SWE(Kinde));	%%[mm]
V_space     =	sum(V,2);	
Vice_space	=	sum(Vice,2);	
WAT_tg      =	mean(WAT(Kinde));	%%[mm]
WIS_tg      =	mean(WIS(Kinde));	%%[mm]rface_albedo
WR_IP_tg	=	mean(WR_IP(Kinde));	%%[mm]
WR_SP_tg	=	mean(WR_SP(Kinde));	%%[mm]
Ws_tg       =	mean(Ws_S(Kinde));	%%[m/s]
ZWT_tg      =	mean(ZWT(Kinde));	%%[mm]
%%%
EIn_tg      =	EIn_H_tg+ EIn_L_tg;	%%[mm/h]
Inveg_tg	=	mean(In_H_space(Kinde) +  In_L_space(Kinde));
O_space     =	V_space./Zs_OUT + Ohy_OUT;	%%[-]
O_tg        =	mean(O_space(Kinde));	%%[-]
PAR_space	=	PARB_S + PARD_S;	%%[W/m^2]
PAR_tg      =	mean(PAR_space(Kinde));%%[W/m^2]
Rsw_space	=	SAB1_S+ SAB2_S + SAD1_S+ SAD2_S;
Rsw_tg      =	mean(Rsw_space(Kinde));	%%[W/m^2]
T_tg        =	T_L_tg +T_H_tg;	%%[mm/h]
V_tg        =	mean(Asur(Kinde).*V_space(Kinde));	%%[mm]
Vice_tg     =	mean(Asur(Kinde).*Vice_space(Kinde));	%%[mm]
In_tg       =	mean(In_H_space(Kinde) +  In_L_space(Kinde) +  SP_wc(Kinde) + In_SWE(Kinde) + In_urb(Kinde) + In_rock(Kinde)+ IP_wc(Kinde) );	%%[mm]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SIMPLE MASS CHECK CONTROL
if sum(SNn)>=1
    CKt = (V_tgtm1 - V_tg) + (Vice_tgtm1 - Vice_tg) + Pr_tg - EG_tg*dth - T_tg*dth - EIn_tg*dth - Lk_tg...
        - ESN_tg*dth - EIn_urb_tg*dth - EWAT_tg*dth - EIn_rock_tg*dth - EICE_tg*dth ...
        + (SWE_tgtm1 -SWE_tg) + (In_tgtm1 -In_tg) ...
        +  (ICE_tgtm1 -ICE_tg) +  (WAT_tgtm1 -WAT_tg)  +  (FROCK_tgtm1 -FROCK_tg) ...
        + Qlat_in_tgtm1 + q_runon_tgtm1 + Q_channel_tgtm1  ...
        - Qlat_in_tg - Q_exit - Qsub_exit - q_runon_tg - Q_channel_tg -Swe_exit;
else
    CKt = (V_tgtm1 - V_tg) + (Vice_tgtm1 - Vice_tg) + Pr_tg - EG_tg*dth - T_tg*dth - EIn_tg*dth - Lk_tg...
        - ESN_tg*dth - EIn_urb_tg*dth - EWAT_tg*dth - EIn_rock_tg*dth - EICE_tg*dth ...
        + (SWE_tgtm1 -SWE_tg) + (In_tgtm1 -In_tg) ...
        +  (ICE_tgtm1 -ICE_tg) +  (WAT_tgtm1 -WAT_tg)  +  (FROCK_tgtm1 -FROCK_tg) ...
        + Qlat_in_tgtm1 + q_runon_tgtm1   ...
        - Qlat_in_tg - Q_exit - Qsub_exit - q_runon_tg  -Swe_exit;
end

%%%%%%%%%%%%%%%%%%% Variable names %%%%%%%%%%%%%%%%%%%%%%%
vars_avg = {'Date','alp_soil_tg','Ca_tg','Cicew_tg','Cice_tg','CK1_tg','Csnow_tg','Csno_tg','DQ_S_tg','dQ_S_tg','Dr_H_tg',...
        'Dr_L_tg','Ds_tg','DT_S_tg','dw_SNO_tg','ea_tg','EG_tg','EICE_tg','EIn_H_tg','EIn_L_tg','EIn_rock_tg','EIn_tg',...
        'EIn_urb_tg','er_tg','ESN_tg','SSN_tg','EWAT_tg','Fract_sat_tg','FROCK_tg','f_tg','Gfin_tg','G_tg','H_tg','ICE_D_tg',...
        'ICE_tg','Imelt_tg','Inveg_tg','In_rock_tg','In_SWE_tg','In_tg','In_urb_tg','IP_wc_tg','Lk_rock_tg','Lk_tg',...
        'Lk_wat_tg','NIce_tg','NIn_SWE_tg','N_tg','NDVI_tg','OF_tg','OS_tg','O_tg','PAR_tg','Pre_tg','Pr_liq_tg','Pr_sno_tg','Pr_tg',...
        'QE_tg','Qfm_tg','Qlat_in_tg','Qlat_out_tg','Qsub_exit','Qv_tg','Q_channel_tg','Q_exit','q_runon_tg','ra_tg',...
        'Rd_tg','Rh_tg','Rn_tg','ros_tg','Rsw_tg','r_soil_tg','SE_rock_tg','SE_urb_tg','Smelt_tg','SND_tg','snow_albedo_tg',...
        'SP_wc_tg','SWE_avalanched_tg','SWE_tg','t','Ta_tg','Tdamp_tg','Tdew_tg','Tdp_tg','Tice_tg','TsVEG_tg','Ts_tg','T_H_tg',...
        'T_L_tg','T_tg','U_SWE_tg','Vice_tg','V_tg','WAT_tg','WIS_tg','WR_IP_tg','WR_SP_tg','Ws_tg','ZWT_tg'};

%%%%%%%%%%%%%%%%% Current date %%%%%%%%%%%%%%%%%%%%%%%%%%%

Date_str = datestr(datetime(Datam_S(1),Datam_S(2), Datam_S(3), Datam_S(4),0,0),'dd-mmm-yyyy HH:MM:ss');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output_manag(1) == 1

%%% [Output: columns in alphabetical order] %%%
if t==2
    tit{1}=strcat(outlocation, '/OUTPUT_',SITE,'_AVG.dat'); 
    fid(1)=fopen(tit{1},'a');

    % Add labels to column list
    for ii = 1:length(vars_avg)-1
        fprintf(fid(1),'%s\t',vars_avg{ii});
    end 
    fprintf(fid(1),'%s\t\n',vars_avg{length(vars_avg)});
end

if t==t1_reinit % if re-starting a T&C run at timestep t1_reinit

   tit{1}=strcat(outlocation, '/OUTPUT_',SITE,'_AVG.dat'); % file name

   % Check the line number after which delete everything (based on datestamp)
   fid(1)=fopen(tit{1},'r+');  % Open the file
   tline = fgetl(fid(1)); %Get the first line
   lineCounter = 1;
   Date_str_m1 = datestr(datetime(Date_str) - hours(1)); %Get the datestamp of the last line to keep
    while ischar(tline)
       if length(tline) > 19 && strcmp(tline(1:20), Date_str_m1)
         break;
       end

      % Read next line
      tline = fgetl(fid(1));
      lineCounter = lineCounter + 1;
    end

   % Copy all the lines until the line number found in the previous section
     fid(1) = fopen(tit{1},'r');
     your_text = cell(lineCounter,1);
     for ii = 1:lineCounter
       your_text(ii) = {fgetl(fid(1))}; 
     end
     fclose(fid(1));

    % Write all the lines until the line number found in the previous section
      fid(1) = fopen(tit{1},'w');
        for ii = 1:lineCounter
          fprintf(fid(1),'%s\n',your_text{ii});
        end
      clear your_text 
end 

%%% START <<OUTPUT_ZZZ_AVG.dat>> %%%
fprintf(fid(1),'%s\t',Date_str);
fprintf(fid(1),'%g\t',alp_soil_tg);
fprintf(fid(1),'%g\t',Ca_tg);      
fprintf(fid(1),'%g\t',Cicew_tg);	
fprintf(fid(1),'%g\t',Cice_tg);    
fprintf(fid(1),'%g\t',CK1_tg);     
fprintf(fid(1),'%g\t',Csnow_tg);	
fprintf(fid(1),'%g\t',Csno_tg);    
fprintf(fid(1),'%g\t',DQ_S_tg);    
fprintf(fid(1),'%g\t',dQ_S_tg);    
fprintf(fid(1),'%g\t',Dr_H_tg);    
fprintf(fid(1),'%g\t',Dr_L_tg);    
fprintf(fid(1),'%g\t',Ds_tg);      
fprintf(fid(1),'%g\t',DT_S_tg);    
fprintf(fid(1),'%g\t',dw_SNO_tg);	
fprintf(fid(1),'%g\t',ea_tg);      
fprintf(fid(1),'%g\t',EG_tg);      
fprintf(fid(1),'%g\t',EICE_tg);    
fprintf(fid(1),'%g\t',EIn_H_tg);	
fprintf(fid(1),'%g\t',EIn_L_tg);	
fprintf(fid(1),'%g\t',EIn_rock_tg);
fprintf(fid(1),'%g\t',EIn_tg);     
fprintf(fid(1),'%g\t',EIn_urb_tg);	
fprintf(fid(1),'%g\t',er_tg);      
fprintf(fid(1),'%g\t',ESN_tg);     
fprintf(fid(1),'%g\t',SSN_tg);     
fprintf(fid(1),'%g\t',EWAT_tg);    
fprintf(fid(1),'%g\t',Fract_sat_tg);
fprintf(fid(1),'%g\t',FROCK_tg);	
fprintf(fid(1),'%g\t',f_tg);       
fprintf(fid(1),'%g\t',Gfin_tg);     
fprintf(fid(1),'%g\t',G_tg);        
fprintf(fid(1),'%g\t',H_tg);        
fprintf(fid(1),'%g\t',ICE_D_tg);	
fprintf(fid(1),'%g\t',ICE_tg);      
fprintf(fid(1),'%g\t',Imelt_tg);	
fprintf(fid(1),'%g\t',Inveg_tg);	
fprintf(fid(1),'%g\t',In_rock_tg);	
fprintf(fid(1),'%g\t',In_SWE_tg);	
fprintf(fid(1),'%g\t',In_tg);       
fprintf(fid(1),'%g\t',In_urb_tg);	
fprintf(fid(1),'%g\t',IP_wc_tg);	   
fprintf(fid(1),'%g\t',Lk_rock_tg);	
fprintf(fid(1),'%g\t',Lk_tg);       
fprintf(fid(1),'%g\t',Lk_wat_tg);
fprintf(fid(1),'%g\t',NIce_tg);     
fprintf(fid(1),'%g\t',NIn_SWE_tg);	
fprintf(fid(1),'%g\t',N_tg);    
fprintf(fid(1),'%g\t',NDVI_tg);   
fprintf(fid(1),'%g\t',OF_tg);       
fprintf(fid(1),'%g\t',OS_tg);       
fprintf(fid(1),'%g\t',O_tg);        
fprintf(fid(1),'%g\t',PAR_tg);      
fprintf(fid(1),'%g\t',Pre_tg);      
fprintf(fid(1),'%g\t',Pr_liq_tg);	
fprintf(fid(1),'%g\t',Pr_sno_tg);	
fprintf(fid(1),'%g\t',Pr_tg);       
fprintf(fid(1),'%g\t',QE_tg);       
fprintf(fid(1),'%g\t',Qfm_tg);      
fprintf(fid(1),'%g\t',Qlat_in_tg);	
fprintf(fid(1),'%g\t',Qlat_out_tg);	
fprintf(fid(1),'%g\t',Qsub_exit);	
fprintf(fid(1),'%g\t',Qv_tg);       
fprintf(fid(1),'%g\t',Q_channel_tg);
fprintf(fid(1),'%g\t',Q_exit);      
fprintf(fid(1),'%g\t',q_runon_tg);	
fprintf(fid(1),'%g\t',ra_tg);       
fprintf(fid(1),'%g\t',Rd_tg);       
fprintf(fid(1),'%g\t',Rh_tg);       
fprintf(fid(1),'%g\t',Rn_tg);       
fprintf(fid(1),'%g\t',ros_tg);      
fprintf(fid(1),'%g\t',Rsw_tg);      
fprintf(fid(1),'%g\t',r_soil_tg);	
fprintf(fid(1),'%g\t',SE_rock_tg);	
fprintf(fid(1),'%g\t',SE_urb_tg);	
fprintf(fid(1),'%g\t',Smelt_tg);	
fprintf(fid(1),'%g\t',SND_tg);      
fprintf(fid(1),'%g\t',snow_albedo_tg);
fprintf(fid(1),'%g\t',SP_wc_tg);
fprintf(fid(1),'%g\t',SWE_avalanched_tg);	
fprintf(fid(1),'%g\t',SWE_tg);      
fprintf(fid(1),'%g\t',t);           
fprintf(fid(1),'%g\t',Ta_tg);       
fprintf(fid(1),'%g\t',Tdamp_tg);	
fprintf(fid(1),'%g\t',Tdew_tg);     
fprintf(fid(1),'%g\t',Tdp_tg);      
fprintf(fid(1),'%g\t',Tice_tg);     
fprintf(fid(1),'%g\t',TsVEG_tg);	
fprintf(fid(1),'%g\t',Ts_tg);       
fprintf(fid(1),'%g\t',T_H_tg);      
fprintf(fid(1),'%g\t',T_L_tg);      
fprintf(fid(1),'%g\t',T_tg);        
fprintf(fid(1),'%g\t',U_SWE_tg);	
fprintf(fid(1),'%g\t',Vice_tg);     
fprintf(fid(1),'%g\t',V_tg);        
fprintf(fid(1),'%g\t',WAT_tg);      
fprintf(fid(1),'%g\t',WIS_tg);      
fprintf(fid(1),'%g\t',WR_IP_tg);	
fprintf(fid(1),'%g\t',WR_SP_tg);	
fprintf(fid(1),'%g\t',Ws_tg);       
fprintf(fid(1),'%g\t\n',ZWT_tg);

%%% END <<OUTPUT_ZZZ_AVG.dat>> %%%
if t==N_time_step
    fclose(fid(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL STANDARD DEVIATION OVER THE WATERSHED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output_manag(2) == 1

std_alp_soil_tg	=	std(alp_soil(Kinde));	%%[]
std_Ca_tg       =	std(Ca_S(Kinde));	%%[W/m^2]
std_Ds_tg       =	std(Ds_S(Kinde));	%%[°C]
std_dw_SNO_tg	=	std(dw_SNO(Kinde));	%%[]
std_ea_tg       =	std(ea_S(Kinde));	%%[Pa]
std_EG_tg       =	std(EG(Kinde));	%%[mm/h]
std_EICE_tg     =	std(EICE(Kinde));	%%[mm/h]
std_EIn_H_tg	=	std(EIn_H_space(Kinde));	%%[mm/h]
std_EIn_L_tg	=	std(EIn_L_space(Kinde));	%%[mm/h]
std_EIn_rock_tg	=	std(EIn_rock(Kinde));	%%[mm/h]
std_EIn_urb_tg	=	std(EIn_urb(Kinde));	%%[mm/h]
std_er_tg       =	std(er(Kinde));	%%[kg/s m^2]
std_ESN_tg      =	std(ESN(Kinde)+ESN_In(Kinde));	%%[mm/h]
std_SSN_tg      =	std(SSN(Kinde)+SSN_In(Kinde));	%%[mm/h]
std_EWAT_tg     =	std(EWAT(Kinde));	%%[mm/h]
std_f_tg        =	std(f(Kinde)*dth);	%%[mm]
std_FROCK_tg	=	std(FROCK(Kinde));	%%[mm]
std_G_tg        =	std(G(Kinde));	%%[W/m^2]
std_Gfin_tg     =	std(Gfin(Kinde));	%%[W/m^2]
std_H_tg        =	std(H(Kinde));	%%[W/m^2]
std_ICE_D_tg	=	std(ICE_D(Kinde));	%%[m]
std_ICE_tg      =	std(ICE(Kinde));	%%[mm]
std_In_SWE_tg	=	std(In_SWE(Kinde));	%%[mm]
std_In_tg       =	std(In_H_space(Kinde) +  In_L_space(Kinde) +  SP_wc(Kinde) + In_SWE(Kinde) +  In_urb(Kinde) + In_rock(Kinde)+ IP_wc(Kinde) );	%%[mm]
std_Inveg_tg	=	std(In_H_space(Kinde) +  In_L_space(Kinde));	
std_IP_wc_tg	=	std(IP_wc(Kinde));	%%[mm]
std_Lk_rock_tg	=	std(Lk_rock(Kinde)*dth);	%%[mm]
std_Lk_tg       =	std(Lk(Kinde)*dth);	%%[mm]
std_Lk_wat_tg	=	std(Lk_wat(Kinde)*dth);	%%[mm]
std_N_tg        =	std(N_S(Kinde));	%%[-]
std_NDVI_tg     =   std(NDVI(Kinde));
std_NIce_tg     =	std(NIce(Kinde));	%%
std_NIn_SWE_tg	=	std(NIn_SWE(Kinde));	%%[m/s]
std_O_tg        =	std(O_space(Kinde));	%%
std_OF_tg       =	std(OF(Kinde));	%%[]
std_OS_tg       =	std(OS(Kinde));	%%[]
std_PAR_tg      =	std(PAR_space(Kinde));	%%[W/m^2]
std_Pr_liq_tg	=	std(Pr_liq(Kinde));	%%[mm]
std_Pr_sno_tg	=	std(Pr_sno(Kinde));	%%[mm]
std_Pr_tg       =	std(Pr_S(Kinde));	%%[mm]
std_Pre_tg      =	std(Pre_S(Kinde));	%%[mbar]
std_Q_channel_tg =	std(Q_channel(Kinde));	%%[mm]
std_q_runon_tg	=	std(q_runon(Kinde)*dth);	%%[mm]
std_QE_tg       =	std(QE(Kinde));	%%[W/m^2]
std_Qfm_tg      =	std(Qfm(Kinde));	%%[W/m^2]
std_Qlat_in_tg	=	std(Qi_in_space(Kinde));	%%[mm]
std_Qlat_out_tg	=	std(Qi_out_space(Kinde));	%%[mm]
std_Qv_tg       =	std(Qv(Kinde));	%%[W/m^2]
std_r_soil_tg	=	std(r_soil(Kinde));	%%[°C]
std_ra_tg       =	std(ra(Kinde));	%%[m/s]
std_Rd_tg       =	std(Rd(Kinde));	%%[mm]
std_Rh_tg       =	std(Rh(Kinde));	%%[mm]
std_Rn_tg       =	std(Rn(Kinde));	%%[W/m^2]
std_ros_tg      =	std(ros(Kinde));	%%[kg/m^3]
std_Rsw_tg      =	std(Rsw_space(Kinde));	%%[W/m^2]
std_SND_tg      =	std(SND(Kinde));	%%[m]
std_snow_albedo_tg = std(snow_albedo(Kinde)); %%[-]
std_SP_wc_tg	=	std(SP_wc(Kinde));	%%[mm]
std_SWE_tg      =	std(SWE(Kinde));	%%[mm]
std_SWE_avalanched_tg =	std(SWE_avalanched(Kinde));	%%[mm]
std_T_H_tg      =	std(T_H_space(Kinde));	%%[mm/h]
std_T_L_tg      =	std(T_L_space(Kinde));	%%[mm/h]
std_Ta_tg       =	std(Ta_S(Kinde));	%%[°C]
std_Tdamp_tg	=	std(Tdamp(Kinde));	%%[°C]
std_Tdew_tg     =	std(Tdew_S(Kinde));	%%[°C]
std_Tdp_tg      =	std(Tdp_space(Kinde));	
std_Ts_tg       =	std(Ts(Kinde));	%%[°C]
std_TsVEG_tg	=	std(TsVEG(Kinde));	%%[°C]
std_U_SWE_tg	=	std(U_SWE(Kinde));	%%[]
std_V_tg        =	std(V_space(Kinde));	%%
std_WAT_tg      =	std(WAT(Kinde));	%%[mm]
std_WIS_tg      =	std(WIS(Kinde));	%%[mm]
std_WR_IP_tg	=	std(WR_IP(Kinde));	%%[]
std_WR_SP_tg	=	std(WR_SP(Kinde));	%%[]
std_Ws_tg       =	std(Ws_S(Kinde));	%%[m/s]
std_ZWT_tg      =	std(ZWT(Kinde));	%%[mm]

%%%%%%%%%%%%%% STD variable names %%%%%%%%%%%%

vars_std = {'Date','std_alp_soil_tg','std_Ca_tg','std_Ds_tg','std_dw_SNO_tg','std_ea_tg','std_EG_tg','std_EICE_tg',...
        'std_EIn_H_tg','std_EIn_L_tg','std_EIn_rock_tg','std_EIn_urb_tg','std_er_tg','std_ESN_tg','std_SSN_tg','std_EWAT_tg','std_FROCK_tg',...
        'std_f_tg','std_Gfin_tg','std_G_tg','std_H_tg','std_ICE_D_tg','std_ICE_tg','std_Inveg_tg','std_In_SWE_tg','std_In_tg',...
        'std_IP_wc_tg','std_Lk_rock_tg','std_Lk_tg','std_Lk_wat_tg','std_NIce_tg','std_NIn_SWE_tg','std_N_tg','std_NDVI_tg','std_OF_tg',...
        'std_OS_tg','std_O_tg','std_PAR_tg','std_Pre_tg','std_Pr_liq_tg','std_Pr_sno_tg','std_Pr_tg','std_QE_tg','std_Qfm_tg',...
        'std_Qlat_in_tg','std_Qlat_out_tg','std_Qv_tg','std_Q_channel_tg','std_q_runon_tg','std_ra_tg','std_Rd_tg',...
        'std_Rh_tg','std_Rn_tg','std_ros_tg','std_Rsw_tg','std_r_soil_tg','std_SND_tg','std_snow_albedo_tg','std_SP_wc_tg',...
        'std_SWE_avalanched_tg','std_SWE_tg','std_Ta_tg','std_Tdamp_tg','std_Tdew_tg','std_Tdp_tg','std_TsVEG_tg',...
        'std_Ts_tg','std_T_H_tg','std_T_L_tg','std_U_SWE_tg','std_V_tg','std_WAT_tg','std_WIS_tg','std_WR_IP_tg',...
        'std_WR_SP_tg','std_Ws_tg','std_ZWT_tg'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [Output: columns in alphabetical order] %%%
if t==2
    tit4{1}=strcat(outlocation,'OUTPUT_',SITE,'_STD.dat');
    fid4(1)=fopen(tit4{1},'a');

        % Add labels to column list
    for ii = 1:length(vars_std)-1
        fprintf(fid4(1),'%s\t',vars_std{ii});
    end 
    fprintf(fid4(1),'%s\t\n',vars_std{length(vars_std)});
end

if t==t1_reinit

   tit4{1}=strcat(outlocation,'OUTPUT_',SITE,'_STD.dat');

   % Check the line number after which delete everything (based on datestamp)
   fid4(1)=fopen(tit4{1},'r+');  % Open the file
   tline = fgetl(fid4(1)); %Get the first line
   lineCounter = 1;
   Date_str_m1 = datestr(datetime(Date_str) - hours(1)); %Get the datestamp of the last line to keep
    while ischar(tline)
       if length(tline) > 19 && strcmp(tline(1:20), Date_str_m1)
         break;
       end
      % Read next line
      tline = fgetl(fid4(1));
      lineCounter = lineCounter + 1;
    end

   % Copy all the lines until the line number found in the previous section
     fid4(1) = fopen(tit4{1},'r');
     your_text = cell(lineCounter,1);
     for ii = 1:lineCounter
       your_text(ii) = {fgetl(fid4(1))}; 
     end
     fclose(fid4(1));

    % Write all the lines until the line number found in the previous section
      fid4(1) = fopen(tit4{1},'w');
        for ii = 1:lineCounter
          fprintf(fid4(1),'%s\n',your_text{ii});
        end
      clear your_text 
end 


%%% START <<OUTPUT_ZZZ_STD.dat>> %%%
fprintf(fid4(1),'%s\t',Date_str);
fprintf(fid4(1),'%g\t',std_alp_soil_tg);	
fprintf(fid4(1),'%g\t',std_Ca_tg);          
fprintf(fid4(1),'%g\t',std_Ds_tg);          
fprintf(fid4(1),'%g\t',std_dw_SNO_tg);		
fprintf(fid4(1),'%g\t',std_ea_tg);          
fprintf(fid4(1),'%g\t',std_EG_tg);          
fprintf(fid4(1),'%g\t',std_EICE_tg);		
fprintf(fid4(1),'%g\t',std_EIn_H_tg);		
fprintf(fid4(1),'%g\t',std_EIn_L_tg);		
fprintf(fid4(1),'%g\t',std_EIn_rock_tg);	
fprintf(fid4(1),'%g\t',std_EIn_urb_tg);		
fprintf(fid4(1),'%g\t',std_er_tg);          
fprintf(fid4(1),'%g\t',std_ESN_tg);         
fprintf(fid4(1),'%g\t',std_SSN_tg);         
fprintf(fid4(1),'%g\t',std_EWAT_tg);		
fprintf(fid4(1),'%g\t',std_FROCK_tg);		
fprintf(fid4(1),'%g\t',std_f_tg);           
fprintf(fid4(1),'%g\t',std_Gfin_tg);		
fprintf(fid4(1),'%g\t',std_G_tg);           
fprintf(fid4(1),'%g\t',std_H_tg);           
fprintf(fid4(1),'%g\t',std_ICE_D_tg);		
fprintf(fid4(1),'%g\t',std_ICE_tg);         
fprintf(fid4(1),'%g\t',std_Inveg_tg);		
fprintf(fid4(1),'%g\t',std_In_SWE_tg);		
fprintf(fid4(1),'%g\t',std_In_tg);          
fprintf(fid4(1),'%g\t',std_IP_wc_tg);		
fprintf(fid4(1),'%g\t',std_Lk_rock_tg);		
fprintf(fid4(1),'%g\t',std_Lk_tg);          
fprintf(fid4(1),'%g\t',std_Lk_wat_tg);		
fprintf(fid4(1),'%g\t',std_NIce_tg);		
fprintf(fid4(1),'%g\t',std_NIn_SWE_tg);		
fprintf(fid4(1),'%g\t',std_N_tg);  
fprintf(fid4(1),'%g\t',std_NDVI_tg);  
fprintf(fid4(1),'%g\t',std_OF_tg);          
fprintf(fid4(1),'%g\t',std_OS_tg);          
fprintf(fid4(1),'%g\t',std_O_tg);           
fprintf(fid4(1),'%g\t',std_PAR_tg);         
fprintf(fid4(1),'%g\t',std_Pre_tg);         
fprintf(fid4(1),'%g\t',std_Pr_liq_tg);		
fprintf(fid4(1),'%g\t',std_Pr_sno_tg);		
fprintf(fid4(1),'%g\t',std_Pr_tg);          
fprintf(fid4(1),'%g\t',std_QE_tg);          
fprintf(fid4(1),'%g\t',std_Qfm_tg);         
fprintf(fid4(1),'%g\t',std_Qlat_in_tg);		
fprintf(fid4(1),'%g\t',std_Qlat_out_tg);	
fprintf(fid4(1),'%g\t',std_Qv_tg);          
fprintf(fid4(1),'%g\t',std_Q_channel_tg);	
fprintf(fid4(1),'%g\t',std_q_runon_tg);		
fprintf(fid4(1),'%g\t',std_ra_tg);          
fprintf(fid4(1),'%g\t',std_Rd_tg);          
fprintf(fid4(1),'%g\t',std_Rh_tg);          
fprintf(fid4(1),'%g\t',std_Rn_tg);          
fprintf(fid4(1),'%g\t',std_ros_tg);         
fprintf(fid4(1),'%g\t',std_Rsw_tg);         
fprintf(fid4(1),'%g\t',std_r_soil_tg);		
fprintf(fid4(1),'%g\t',std_SND_tg);         
fprintf(fid4(1),'%g\t',std_snow_albedo_tg); 
fprintf(fid4(1),'%g\t',std_SP_wc_tg);		
fprintf(fid4(1),'%g\t',std_SWE_avalanched_tg);    
fprintf(fid4(1),'%g\t',std_SWE_tg);         
fprintf(fid4(1),'%g\t',std_Ta_tg);          
fprintf(fid4(1),'%g\t',std_Tdamp_tg);		
fprintf(fid4(1),'%g\t',std_Tdew_tg);		
fprintf(fid4(1),'%g\t',std_Tdp_tg);         
fprintf(fid4(1),'%g\t',std_TsVEG_tg);		
fprintf(fid4(1),'%g\t',std_Ts_tg);          
fprintf(fid4(1),'%g\t',std_T_H_tg);         
fprintf(fid4(1),'%g\t',std_T_L_tg);         
fprintf(fid4(1),'%g\t',std_U_SWE_tg);		
fprintf(fid4(1),'%g\t',std_V_tg);           
fprintf(fid4(1),'%g\t',std_WAT_tg);         
fprintf(fid4(1),'%g\t',std_WIS_tg);         
fprintf(fid4(1),'%g\t',std_WR_IP_tg);		
fprintf(fid4(1),'%g\t',std_WR_SP_tg);		
fprintf(fid4(1),'%g\t',std_Ws_tg);          
fprintf(fid4(1),'%g\t\n',std_ZWT_tg);
%%% END <<OUTPUT_ZZZ_STD.dat>> %%%
if t==N_time_step
    fclose(fid4(1));
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL AVERAGE OVER THE VEGETATION TYPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% START VEGETATION TYPES %%%
%%% 1 = Fir
%%% 2 = Larch
%%% 3 = Grass
%%% 4 = Shrub
%%% 5 = Broadleaf-evergreen
%%% 6 = Broadleaf-deciduous
%%% 7 = Rock/Ice
if output_manag(3) == 1

Veg_names = {'Fir','Larch','Grass','Shrub','Broadleaf-evergreen','Broadleaf-deciduous','Rock/Ice'};
%%% END VEGETATION TYPES %%%

for ijki=1:cc_max
    for ievc=1:length(EVcode)
        riev= find(ksv==EVcode(ievc));
        %%%%%%
        OH_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*OH(riev,ijki)); %%[]
        std_OH_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*OH(riev,ijki)); %%[]
        OL_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*OL(riev,ijki)); %%[]
        std_OL_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*OL(riev,ijki)); %%[]
        An_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*An_H(riev,ijki)); %%[]
        std_An_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*An_H(riev,ijki)); %%[]
        An_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*An_L(riev,ijki)); %%[]
        std_An_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*An_L(riev,ijki)); %%[]
        Tdp_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Tdp_H(riev,ijki)); %%[]
        std_Tdp_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Tdp_H(riev,ijki)); %%[]
        Tdp_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Tdp_L(riev,ijki)); %%[]
        std_Tdp_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Tdp_L(riev,ijki)); %%[]
        Rdark_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Rdark_H(riev,ijki)); %%[]
        std_Rdark_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Rdark_H(riev,ijki)); %%[]
        Rdark_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Rdark_L(riev,ijki)); %%[]
        std_Rdark_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Rdark_L(riev,ijki)); %%[]
        LAI_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*LAI_H(riev,ijki)); %%[]
        std_LAI_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*LAI_H(riev,ijki)); %%[]
        LAI_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*LAI_L(riev,ijki)); %%[]
        std_LAI_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*LAI_L(riev,ijki)); %%[]
        NDVI_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*NDVI(riev,ijki)); %%[]
        std_NDVI_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*NDVI(riev,ijki)); %%[]
        NPP_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*NPP_H(riev,ijki)); %%[]
        std_NPP_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*NPP_H(riev,ijki)); %%[]
        NPP_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*NPP_L(riev,ijki)); %%[]
        std_NPP_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*NPP_L(riev,ijki)); %%[]
        ANPP_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*ANPP_H(riev,ijki)); %%[]
        std_ANPP_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*ANPP_H(riev,ijki)); %%[]
        ANPP_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*ANPP_L(riev,ijki)); %%[]
        std_ANPP_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*ANPP_L(riev,ijki)); %%[]
        RA_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*RA_H(riev,ijki)); %%[]
        std_RA_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*RA_H(riev,ijki)); %%[]
        RA_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*RA_L(riev,ijki)); %%[]
        std_RA_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*RA_L(riev,ijki)); %%[]
        Rg_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Rg_H(riev,ijki)); %%[]
        std_Rg_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Rg_H(riev,ijki)); %%[]
        Rg_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Rg_L(riev,ijki)); %%[]
        std_Rg_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Rg_L(riev,ijki)); %%[]
        LAIdead_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*LAIdead_H(riev,ijki)); %%[]
        std_LAIdead_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*LAIdead_H(riev,ijki)); %%[]
        LAIdead_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*LAIdead_L(riev,ijki)); %%[]
        std_LAIdead_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*LAIdead_L(riev,ijki)); %%[]
        hc_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*hc_H(riev,ijki)); %%[]
        std_hc_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*hc_H(riev,ijki)); %%[]
        hc_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*hc_L(riev,ijki)); %%[]
        std_hc_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*hc_L(riev,ijki)); %%[]
        AgeL_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*AgeL_H(riev,ijki)); %%[]
        std_AgeL_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*AgeL_H(riev,ijki)); %%[]
        AgeL_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*AgeL_L(riev,ijki)); %%[]
        std_AgeL_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*AgeL_L(riev,ijki)); %%[]
        SAI_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*SAI_H(riev,ijki)); %%[]
        std_SAI_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*SAI_H(riev,ijki)); %%[]
        SAI_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*SAI_L(riev,ijki)); %%[]
        std_SAI_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*SAI_L(riev,ijki)); %%[]
        PHE_S_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*PHE_S_H(riev,ijki)); %%[]
        std_PHE_S_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*PHE_S_H(riev,ijki)); %%[]
        PHE_S_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*PHE_S_L(riev,ijki)); %%[]
        std_PHE_S_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*PHE_S_L(riev,ijki)); %%[]
        %Llitter_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Llitter(riev,ijki)); %%[]
        %std_Llitter_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Llitter(riev,ijki)); %%[]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %,B_H(ij,:,:),,,
        %,B_L(ij,:,:),
    end
end

%%%%%% Vegetation vars name

vars_veg = {'Date','AgeL_H_tg','AgeL_L_tg','ANPP_H_tg','ANPP_L_tg','An_H_tg','An_L_tg','hc_H_tg','LAIdead_H_tg','LAIdead_L_tg',...
        'LAI_H_tg','LAI_L_tg','NDVI_tg','NPP_H_tg','NPP_L_tg','OH_tg','OL_tg','PHE_S_H_tg','PHE_S_L_tg','RA_H_tg','RA_L_tg','Rdark_H_tg',...
        'Rdark_L_tg','Rg_H_tg','Rg_L_tg','SAI_H_tg','SAI_L_tg','Tdp_H_tg','Tdp_L_tg','hc_L_tg'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [Output: columns in alphabetical order] %%%
if t==2
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0
                tit2{ievc,ijki}=strcat(outlocation,'OUTPUT_',SITE,'_AVG_PFT_',Veg_names{ievc},'.dat');
                fid2(ievc,ijki)=fopen(tit2{ievc,ijki},'a');

                % Add labels to column list
                for ii = 1:length(vars_veg)-1
                   fprintf(fid2(ievc,ijki),'%s\t',vars_veg{ii});
                end 
                fprintf(fid2(ievc,ijki),'%s\t\n',vars_veg{length(vars_veg)});
            end
        end
    end
end


if t==t1_reinit
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0

                  tit2{ievc,ijki}=strcat(outlocation,'OUTPUT_',SITE,'_AVG_PFT_',Veg_names{ievc},'.dat');

                  % Check the line number after which delete everything (based on datestamp)
                  fid2(ievc,ijki)=fopen(tit2{ievc,ijki},'r+');  % Open the file
                  tline = fgetl(fid2(ievc,ijki)); %Get the first line
                  lineCounter = 1;
                  Date_str_m1 = datestr(datetime(Date_str) - hours(1)); %Get the datestamp of the last line to keep
                  while ischar(tline)
                     if length(tline) > 19 && strcmp(tline(1:20), Date_str_m1)
                         break;
                     end
                     % Read next line
                  tline = fgetl(fid2(ievc,ijki));
                  lineCounter = lineCounter + 1;
                  end

                % Copy all the lines until the line number found in the previous section
                 fid2(ievc,ijki) = fopen(tit2{ievc,ijki},'r');
                 your_text = cell(lineCounter,1);
                 for ii = 1:lineCounter
                   your_text(ii) = {fgetl(fid2(ievc,ijki))}; 
                 end
                 fclose(fid2(ievc,ijki));

                % Write all the lines until the line number found in the previous section
                 fid2(ievc,ijki) = fopen(tit2{ievc,ijki},'w');
                      for ii = 1:lineCounter
                            fprintf(fid2(ievc,ijki),'%s\n',your_text{ii});
                      end
                 clear your_text 
            end
       end
    end
end

for ijki=1:cc_max
    for ievc=1:length(EVcode)
        if Ccrown_OUT(ievc,ijki)>0
		    %%% START <<OUTPUT_ZZZ_AVG_PFT_YYY_code_XXX.dat>> %%%
            fprintf(fid2(ievc,ijki),'%s\t',Date_str);
            fprintf(fid2(ievc,ijki),'%g\t',AgeL_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',AgeL_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',ANPP_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',ANPP_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',An_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',An_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',hc_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',LAIdead_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',LAIdead_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',LAI_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',LAI_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',NDVI_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',NPP_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',NPP_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',OH_tg(ievc,ijki)); 
            fprintf(fid2(ievc,ijki),'%g\t',OL_tg(ievc,ijki)); 
            fprintf(fid2(ievc,ijki),'%g\t',PHE_S_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',PHE_S_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',RA_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',RA_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',Rdark_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',Rdark_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',Rg_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',Rg_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',SAI_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',SAI_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',Tdp_H_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t',Tdp_L_tg(ievc,ijki));
            fprintf(fid2(ievc,ijki),'%g\t\n',hc_L_tg(ievc,ijki));
			%%% END <<OUTPUT_ZZZ_AVG_PFT_YYY_code_XXX.dat>> %%%
        end
    end
end

%%%
if t==N_time_step
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0
                fclose(fid2(ievc,ijki));
            end
        end
    end
end

end


%==========================================================================
%==========================================================================
%%%%%%%%% SPATIAL STD OVER THE VEGETATION TYPE
%==========================================================================
%==========================================================================

%%%%%% Vegetation vars name
if output_manag(4) == 1

vars_veg_std = {'Date','std_AgeL_H_tg','std_AgeL_L_tg','std_ANPP_H_tg','std_ANPP_L_tg','std_An_H_tg','std_An_L_tg',...
    'std_hc_H_tg','std_LAIdead_H_tg','std_LAIdead_L_tg','std_LAI_H_tg','std_LAI_L_tg','std_NDVI_tg','std_NPP_H_tg',...
    'std_NPP_L_tg','std_OH_tg','std_OL_tg','std_PHE_S_H_tg','std_PHE_S_L_tg','std_RA_H_tg','std_RA_L_tg','std_Rdark_H_tg',...
        'std_Rdark_L_tg','std_Rg_H_tg','std_Rg_L_tg','std_SAI_H_tg','std_SAI_L_tg','std_Tdp_H_tg','std_Tdp_L_tg','std_hc_L_tg'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [Output: columns in alphabetical order] %%%
if t==2
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0
                tit3{ievc,ijki}=strcat(outlocation,'OUTPUT_',SITE,'_STD_PFT_',Veg_names{ievc},'.dat');
                fid3(ievc,ijki)=fopen(tit3{ievc,ijki},'a');

                % Add labels to column list
                for ii = 1:length(vars_veg_std)-1
                   fprintf(fid3(ievc,ijki),'%s\t',vars_veg_std{ii});
                end 
                fprintf(fid3(ievc,ijki),'%s\t\n',vars_veg_std{length(vars_veg_std)});
            end
        end
    end
end

if t==t1_reinit
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0

                  tit3{ievc,ijki}=strcat(outlocation,'OUTPUT_',SITE,'_STD_PFT_',Veg_names{ievc},'.dat');

                  % Check the line number after which delete everything (based on datestamp)
                 fid3(ievc,ijki)=fopen(tit3{ievc,ijki},'r+');  % Open the file
                  tline = fgetl(fid3(ievc,ijki)); %Get the first line
                  lineCounter = 1;
                  Date_str_m1 = datestr(datetime(Date_str) - hours(1)); %Get the datestamp of the last line to keep
                  while ischar(tline)
                     if length(tline) > 19 && strcmp(tline(1:20), Date_str_m1)
                         break;
                     end
                     % Read next line
                  tline = fgetl(fid3(ievc,ijki));
                  lineCounter = lineCounter + 1;
                  end

                % Copy all the lines until the line number found in the previous section
                 fid3(ievc,ijki) = fopen(tit3{ievc,ijki},'r');
                 your_text = cell(lineCounter,1);
                 for ii = 1:lineCounter
                   your_text(ii) = {fgetl(fid3(ievc,ijki))}; 
                 end
                 fclose(fid3(ievc,ijki));

                % Write all the lines until the line number found in the previous section
                 fid3(ievc,ijki) = fopen(tit3{ievc,ijki},'w');
                      for ii = 1:lineCounter
                            fprintf(fid3(ievc,ijki),'%s\n',your_text{ii});
                      end
                 clear your_text    
            end
        end
    end
end

for ijki=1:cc_max
    for ievc=1:length(EVcode)
        if Ccrown_OUT(ievc,ijki)>0
		    %%% START <<OUTPUT_ZZZ_STD_PFT_YYY_code_XXX.dat>> %%%
            fprintf(fid3(ievc,ijki),'%s\t',Date_str);
            fprintf(fid3(ievc,ijki),'%g\t',std_AgeL_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_AgeL_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_ANPP_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_ANPP_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_An_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_An_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_hc_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_LAIdead_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_LAIdead_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_LAI_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_LAI_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_LAI_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_NDVI_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_NPP_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_OH_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_OL_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_PHE_S_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_PHE_S_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_RA_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_RA_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_Rdark_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_Rdark_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_Rg_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_Rg_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_SAI_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_SAI_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_Tdp_H_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t',std_Tdp_L_tg(ievc,ijki));
            fprintf(fid3(ievc,ijki),'%g\t\n',std_hc_L_tg(ievc,ijki));
			%%% END <<OUTPUT_ZZZ_STD_PFT_YYY_code_XXX.dat>> %%%
        end
    end
end
%%%
if t==N_time_step
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0
                fclose(fid3(ievc,ijki));
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL AVERAGE OVER EACH LAND COVER CLASS
%%%%%%%%%  (vegetation types, rock, ice, clean-ice, debris-covered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Land cover class index "LCinde"
% 1: idx_Veg1 (veg 1 index (Fir))
% 2: idx_Veg2 (veg 2 index (Larch))
% 3: idx_Veg3 (veg 3 index (Grass))
% 4: idx_Veg4 (veg 4 index (Shrub))
% 5: idx_Veg5 (veg 5 index (Broadleaf-evergreen)
% 6: idx_Veg6 (veg 6 index (Broadleaf-deciduous)
% 7: idx_Rock (rock index)
% 8: idx_Ice (ice index)
% 9: idx_Cleanice (clean-ice index)
% 10: idx_Debice (debris-covered ice index)
if output_manag(5) == 1

LC_names = {'Fir','Larch','Grass','Shrub','Broadleaf-evergreen','Broadleaf-deciduous',...
            'Rock','Ice','Cleanice','Debice'};

%%% START LAND COVER CLASSES %%%
%%% 1 = Fir
%%% 2 = Larch
%%% 3 = Grass
%%% 4 = Shrub
%%% 5 = Broadleaf-evergreen
%%% 6 = Broadleaf-deciduous
%%% 7 = Rock
%%% 8 = Ice
%%% 9 = Clean-Ice
%%% 10 = Debris-covered-ice
%%% END LAND COVER CLASSES %%%

for ilc=1:length(LCinde) %%%loop over classes
    idx= LCinde{ilc}; %index for current class
    %%%%%%
    alp_soil_tg_LC(ilc)	=	mean(alp_soil(idx));	%%[-]
    Ca_tg_LC(ilc)       =	mean(Ca_S(idx));	%%[ppm]
    Cice_tg_LC(ilc)     =	mean(Cice(idx));	%%[]
    Cicew_tg_LC(ilc)	=	mean(Cicew(idx));	%%[]
    CK1_tg_LC(ilc)      =	mean(CK1(idx));	%%[mm]
    Csno_tg_LC(ilc)     =	mean(Csno(idx));	%%[]
    Csnow_tg_LC(ilc)	=	mean(Csnow(idx));	%%[]
    DQ_S_tg_LC(ilc)     =	mean(DQ_S(idx));	%%[]
    dQ_S_tg_LC(ilc)     =	mean(dQ_S(idx));	%%[]
    Dr_H_tg_LC(ilc)     =	mean(sum(Dr_H(idx),2));	%%[mm]
    Dr_L_tg_LC(ilc)     =	mean(sum(Dr_L(idx),2));	%%[mm]
    Ds_tg_LC(ilc)       =	mean(Ds_S(idx));	%%[°C]
    DT_S_tg_LC(ilc)     =	mean(DT_S(idx));	%%[]
    dw_SNO_tg_LC(ilc)	=	mean(dw_SNO(idx));	%%[]
    ea_tg_LC(ilc)       =	mean(ea_S(idx));	%%[Pa]
    EG_tg_LC(ilc)       =	mean(EG(idx));	%%[mm/h]
    EICE_tg_LC(ilc)     =	mean(EICE(idx));	%%[mm/h]
    EIn_H_tg_LC(ilc)	=	mean(sum(EIn_H(idx),2));	%%[mm/h]
    EIn_L_tg_LC(ilc)	=	mean(sum(EIn_L(idx),2));	%%[mm/h]
    EIn_rock_tg_LC(ilc)	=	mean(EIn_rock(idx));	%%[mm/h]
    EIn_urb_tg_LC(ilc)	=	mean(EIn_urb(idx));	%%[mm/h]
    er_tg_LC(ilc)       =	mean(er(idx));	%%[kg/s m^2]
    ESN_tg_LC(ilc)      =	mean(ESN(idx)+ESN_In(idx));	%%[mm/h]
    SSN_tg_LC(ilc)      =	mean(SSN(idx)+SSN_In(idx));	%%[mm/h]
    EWAT_tg_LC(ilc)     =	mean(EWAT(idx));	%%[mm/h]
    f_tg_LC(ilc)        =	mean(f(idx)*dth);	%%[mm]
    Fract_sat_tg_LC(ilc) =	sum(ZWT(idx)==0)/numel(ZWT(idx));	%%[-]
    FROCK_tg_LC(ilc)	=	mean(FROCK(idx));	%%[mm]
    G_tg_LC(ilc)        =	mean(G(idx));	%%[W/m^2]
    Gfin_tg_LC(ilc)     =	mean(Gfin(idx));	%%[W/m^2]
    H_tg_LC(ilc)        =	mean(H(idx));	%%[W/m^2]
    ICE_D_tg_LC(ilc)	=	mean(ICE_D(idx));	%%[m]
    ICE_tg_LC(ilc)      =	mean(ICE(idx));	%%[mm]
    Imelt_tg_LC(ilc)	=	mean(Imelt(idx));	%%[mm]
    In_rock_tg_LC(ilc)	=	mean(In_rock(idx));	%%[mm]
    In_SWE_tg_LC(ilc)	=	mean(In_SWE(idx));	%%[mm]
    In_urb_tg_LC(ilc)	=	mean(In_urb(idx));	%%[mm]
    IP_wc_tg_LC(ilc)	=	mean(IP_wc(idx));	%%[mm]
    Lk_rock_tg_LC(ilc)	=	mean(Lk_rock(idx)*dth);	%%[mm]
    Lk_tg_LC(ilc)       =	mean(Lk(idx)*dth);	%%[mm]
    Lk_wat_tg_LC(ilc)	=	mean(Lk_wat(idx)*dth);	%%[mm]
    N_tg_LC(ilc)        =	mean(N_S(idx));	%%[-]
    NDVI_tg_LC(ilc)     =   mean(NDVI(idx));
    NIce_tg_LC(ilc)     =	mean(NIce(idx));	%%[mm]
    NIn_SWE_tg_LC(ilc)	=	mean(NIn_SWE(idx));	%%[mm]
    OF_tg_LC(ilc)       =	mean(OF(idx));	%%[]
    OS_tg_LC(ilc)       =	mean(OS(idx));	%%[]
    Pr_liq_tg_LC(ilc)	=	mean(Pr_liq(idx));	%%[mm]
    Pr_sno_tg_LC(ilc)	=	mean(Pr_sno(idx));	%%[mm]
    Pr_tg_LC(ilc)       =	mean(Pr_S(idx));	%%[mm]
    Pre_tg_LC(ilc)      =	mean(Pre_S(idx));	%%[mbar]
    Q_channel_tg_LC(ilc) =	mean(Q_channel(idx));	%%[mm]
    q_runon_tg_LC(ilc)	=	mean(q_runon(idx)*dth);	%%[mm]
    QE_tg_LC(ilc)       =	mean(QE(idx));	%%[W/m^2]
    Qfm_tg_LC(ilc)      =	mean(Qfm(idx));	%%[W/m^2]
    Qlat_in_tg_LC(ilc)	=	mean(sum(Qi_in(idx)*dth,2));	%%[mm]
    Qlat_out_tg_LC(ilc)	=	mean(sum(Qi_out(idx)*dth,2));	%%[mm]
    Qv_tg_LC(ilc)       =	mean(Qv(idx));	%%[W/m^2]
    r_soil_tg_LC(ilc)	=	mean(r_soil(idx));	%%[s/m]
    ra_tg_LC(ilc)       =	mean(ra(idx));	%%[s/m]
    Rd_tg_LC(ilc)       =	mean(Rd(idx));	%%[mm]
    Rh_tg_LC(ilc)       =	mean(Rh(idx));	%%[mm]
    Rn_tg_LC(ilc)       =	mean(Rn(idx));	%%[W/m^2]
    ros_tg_LC(ilc)      =	mean(ros(idx));	%%[kg/m^3]
    SE_rock_tg_LC(ilc)	=	mean(SE_rock(idx));	%%[]
    SE_urb_tg_LC(ilc)	=	mean(SE_urb(idx));	%%[]
    Smelt_tg_LC(ilc)	=	mean(Smelt(idx));	%%[mm]
    SND_tg_LC(ilc)      =	mean(SND(idx));	%%[m]
    snow_albedo_tg_LC(ilc)= mean(snow_albedo(idx)); %%[-]
    SP_wc_tg_LC(ilc)	=	mean(SP_wc(idx));	%%[mm]
    SWE_tg_LC(ilc)      =	mean(SWE(idx));	%%[mm]
    SWE_avalanched_tg_LC(ilc) =	mean(SWE_avalanched(idx));	%%[mm]
    T_H_tg_LC(ilc)      =	mean(sum(T_H(idx),2));	%%[mm/h]
    T_L_tg_LC(ilc)      =	mean(sum(T_L(idx),2));	%%[mm/h]
    Ta_tg_LC(ilc)       =	mean(Ta_S(idx));	%%[°C]
    Tdamp_tg_LC(ilc)	=	mean(Tdamp(idx));	%%[°C]
    Tdew_tg_LC(ilc)     =	mean(Tdew_S(idx));	%%[°C]
    Tdp_tg_LC(ilc)      =	mean(Tdp_space(idx));	%%[°C]
    Tice_tg_LC(ilc)     =	mean(Tice(idx));	%%[C]
    Ts_tg_LC(ilc)       =	mean(Ts(idx));	%%[°C]
    TsVEG_tg_LC(ilc)	=	mean(TsVEG(idx));	%%[°C]
    U_SWE_tg_LC(ilc)	=	mean(U_SWE(idx));	%%[mm]
    WAT_tg_LC(ilc)      =	mean(WAT(idx));	%%[mm]
    WIS_tg_LC(ilc)      =	mean(WIS(idx));	%%[mm]
    WR_IP_tg_LC(ilc)	=	mean(WR_IP(idx));	%%[mm]
    WR_SP_tg_LC(ilc)	=	mean(WR_SP(idx));	%%[mm]
    Ws_tg_LC(ilc)       =	mean(Ws_S(idx));	%%[m/s]
    ZWT_tg_LC(ilc)      =	mean(ZWT(idx));	%%[mm]
    %%%
    EIn_tg_LC(ilc)      =	EIn_H_tg_LC(ilc)+ EIn_L_tg_LC(ilc);	%%[mm/h]
    Inveg_tg_LC(ilc)	=	mean(sum(In_H(idx),2) +  sum(In_L(idx),2));
    O_tg_LC(ilc)        =	mean(O_space(idx));	%%[-]
    PAR_tg_LC(ilc)      =	mean(PAR_space(idx));%%[W/m^2]
    Rsw_tg_LC(ilc)      =	mean(Rsw_space(idx));	%%[W/m^2]
    T_tg_LC(ilc)        =	T_L_tg_LC(ilc) +T_H_tg_LC(ilc);	%%[mm/h]
    V_tg_LC(ilc)        =	mean(Asur(idx).*sum(V(idx),2));	%%[mm]
    Vice_tg_LC(ilc)     =	mean(Asur(idx).*sum(Vice(idx),2));	%%[mm]
    In_tg_LC(ilc)       =	mean(In_H_space(ilc) +  In_L_space(idx) +  SP_wc(idx) + In_SWE(idx) + In_urb(idx) + In_rock(idx)+ IP_wc(idx) );	%%[mm]
end

%%%%%%%%%%%%%% Land cover class variable names %%%%%%%%%%%

vars_LC = {'Date','alp_soil_tg_LC','Ca_tg_LC','Cicew_tg_LC','Cice_tg_LC','CK1_tg_LC','Csnow_tg_LC','Csno_tg_LC','DQ_S_tg_LC',...
    'dQ_S_tg_LC','Dr_H_tg_LC','Dr_L_tg_LC','Ds_tg_LC','DT_S_tg_LC','dw_SNO_tg_LC','ea_tg_LC','EG_tg_LC','EICE_tg_LC','EIn_H_tg_LC',...
    'EIn_L_tg_LC','EIn_rock_tg_LC','EIn_tg_LC','EIn_urb_tg_LC','er_tg_LC','ESN_tg_LC','SSN_tg_LC','EWAT_tg_LC','Fract_sat_tg_LC',...
    'FROCK_tg_LC','f_tg_LC','Gfin_tg_LC','G_tg_LC','H_tg_LC','ICE_D_tg_LC','ICE_tg_LC','Imelt_tg_LC','Inveg_tg_LC','In_rock_tg_LC',...
    'In_SWE_tg_LC','In_urb_tg_LC','IP_wc_tg_LC','Lk_rock_tg_LC','Lk_tg_LC','Lk_wat_tg_LC','NIce_tg_LC','NIn_SWE_tg_LC','N_tg_LC','NDVI_tg_LC',...
    'OF_tg_LC','OS_tg_LC','O_space','O_tg_LC','PAR_space','PAR_tg_LC','Pre_tg_LC','Pr_liq_tg_LC','Pr_sno_tg_LC','Pr_tg_LC',...
    'QE_tg_LC','Qfm_tg_LC','Qlat_in_tg_LC','Qlat_out_tg_LC','Qv_tg_LC','Q_channel_tg_LC','q_runon_tg_LC','ra_tg_LC','Rd_tg_LC',...
    'Rh_tg_LC','Rn_tg_LC','ros_tg_LC','Rsw_space','Rsw_tg_LC','r_soil_tg_LC','SE_rock_tg_LC','SE_urb_tg_LC','Smelt_tg_LC',...
    'SND_tg_LC','snow_albedo_tg_LC','SP_wc_tg_LC','SWE_avalanched_tg_LC','SWE_tg_LC','Ta_tg_LC','Tdamp_tg_LC','Tdew_tg_LC',...
    'Tdp_space','Tdp_tg_LC','Tice_tg_LC','TsVEG_tg_LC','Ts_tg_LC','T_H_tg_LC','T_L_tg_LC','T_tg_LC','U_SWE_tg_LC','Vice_tg_LC',...
    'V_tg_LC','WAT_tg_LC','WIS_tg_LC','WR_IP_tg_LC','WR_SP_tg_LC','Ws_tg_LC','ZWT_tg_LC','In_tg_LC'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [Output: columns in alphabetical order] %%%
if t==2
    for ilc=1:length(LCinde) %%%loop over classes
        tit8{ilc}=strcat(outlocation,'OUTPUT_',SITE,'_AVG_LC_',LC_names{ilc},'.dat');
        fid8(ilc)=fopen(tit8{ilc},'a');

        % Add labels to column list
           for ii = 1:length(vars_LC)-1
               fprintf(fid8(ilc),'%s\t',vars_LC{ii});
           end 
           fprintf(fid8(ilc),'%s\t\n',vars_LC{length(vars_LC)});        
    end
end
if t==t1_reinit
    for ilc=1:length(LCinde) %%%loop over classes

   tit8{ilc}=strcat(outlocation,'OUTPUT_',SITE,'_AVG_LC_',LC_names{ilc},'.dat');

   % Check the line number after which delete everything (based on datestamp)
   fid8(ilc)=fopen(tit8{ilc},'r+');  % Open the file
   tline = fgetl(fid8(ilc)); %Get the first line
   lineCounter = 1;
   Date_str_m1 = datestr(datetime(Date_str) - hours(1)); %Get the datestamp of the last line to keep
    while ischar(tline)
       if length(tline) > 19 && strcmp(tline(1:20), Date_str_m1)
         break;
       end
      % Read next line
      tline = fgetl(fid8(ilc));
      lineCounter = lineCounter + 1;
    end

   % Copy all the lines until the line number found in the previous section
     fid8(ilc) = fopen(tit8{ilc},'r');
     your_text = cell(lineCounter,1);
     for ii = 1:lineCounter
       your_text(ii) = {fgetl(fid8(ilc))}; 
     end
     fclose(fid8(ilc));

    % Write all the lines until the line number found in the previous section
      fid8(ilc) = fopen(tit8{ilc},'w');
        for ii = 1:lineCounter
          fprintf(fid8(ilc),'%s\n',your_text{ii});
        end
      clear your_text               
    end
end

for ilc=1:length(LCinde) %%%loop over classes
    %%% START <<OUTPUT_ZZZ_AVG_LC_code_YYY.dat>> %%%
    fprintf(fid8(ilc),'%s\t',Date_str);
    fprintf(fid8(ilc),'%g\t',alp_soil_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',Ca_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',Cicew_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',Cice_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',CK1_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',Csnow_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',Csno_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',DQ_S_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',dQ_S_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',Dr_H_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',Dr_L_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',Ds_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',DT_S_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',dw_SNO_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',ea_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',EG_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',EICE_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',EIn_H_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',EIn_L_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',EIn_rock_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',EIn_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',EIn_urb_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',er_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',ESN_tg_LC(ilc));  
    fprintf(fid8(ilc),'%g\t',SSN_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',EWAT_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',Fract_sat_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',FROCK_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',f_tg_LC(ilc));        	
    fprintf(fid8(ilc),'%g\t',Gfin_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',G_tg_LC(ilc));        	
    fprintf(fid8(ilc),'%g\t',H_tg_LC(ilc));        	
    fprintf(fid8(ilc),'%g\t',ICE_D_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',ICE_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',Imelt_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',Inveg_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',In_rock_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',In_SWE_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',In_urb_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',IP_wc_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',Lk_rock_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',Lk_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',Lk_wat_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',NIce_tg_LC(ilc));     	
    fprintf(fid8(ilc),'%g\t',NIn_SWE_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',N_tg_LC(ilc));        	
    fprintf(fid8(ilc),'%g\t',NDVI_tg_LC(ilc));        	
    fprintf(fid8(ilc),'%g\t',OF_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',OS_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',O_space(ilc));     	
    fprintf(fid8(ilc),'%g\t',O_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',PAR_space(ilc));	
    fprintf(fid8(ilc),'%g\t',PAR_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',Pre_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',Pr_liq_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',Pr_sno_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',Pr_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',QE_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',Qfm_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',Qlat_in_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',Qlat_out_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',Qv_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',Q_channel_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',q_runon_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',ra_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',Rd_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',Rh_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',Rn_tg_LC(ilc));       	
    fprintf(fid8(ilc),'%g\t',ros_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',Rsw_space(ilc));	
    fprintf(fid8(ilc),'%g\t',Rsw_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',r_soil_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',SE_rock_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',SE_urb_tg_LC(ilc));	
    fprintf(fid8(ilc),'%g\t',Smelt_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',SND_tg_LC(ilc));      	
    fprintf(fid8(ilc),'%g\t',snow_albedo_tg_LC(ilc));
    fprintf(fid8(ilc),'%g\t',SP_wc_tg_LC(ilc));
    fprintf(fid8(ilc),'%g\t',SWE_avalanched_tg_LC(ilc));
    fprintf(fid8(ilc),'%g\t',SWE_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',Ta_tg_LC(ilc));      
    fprintf(fid8(ilc),'%g\t',Tdamp_tg_LC(ilc));   
    fprintf(fid8(ilc),'%g\t',Tdew_tg_LC(ilc));    
    fprintf(fid8(ilc),'%g\t',Tdp_space(ilc));
    fprintf(fid8(ilc),'%g\t',Tdp_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',Tice_tg_LC(ilc));    
    fprintf(fid8(ilc),'%g\t',TsVEG_tg_LC(ilc)); 
    fprintf(fid8(ilc),'%g\t',Ts_tg_LC(ilc));      
    fprintf(fid8(ilc),'%g\t',T_H_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',T_L_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',T_tg_LC(ilc));       
    fprintf(fid8(ilc),'%g\t',U_SWE_tg_LC(ilc));   
    fprintf(fid8(ilc),'%g\t',Vice_tg_LC(ilc));
    fprintf(fid8(ilc),'%g\t',V_tg_LC(ilc));       
    fprintf(fid8(ilc),'%g\t',WAT_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',WIS_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t',WR_IP_tg_LC(ilc));   
    fprintf(fid8(ilc),'%g\t',WR_SP_tg_LC(ilc));   
    fprintf(fid8(ilc),'%g\t',Ws_tg_LC(ilc));      
    fprintf(fid8(ilc),'%g\t',ZWT_tg_LC(ilc));     
    fprintf(fid8(ilc),'%g\t\n',In_tg_LC(ilc));
	%%% END <<OUTPUT_ZZZ_AVG_LC_code_YYY.dat>> %%%
end

%%%
if t==N_time_step
    for ilc=1:length(LCinde)
        fclose(fid8(ilc));
    end
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL STD OVER EACH LAND COVER CLASS
%%%%%%%%%  (vegetation types, rock, ice, clean-ice, debris-covered)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Land cover class index "LCinde"
% 1: idx_Veg1 (veg 1 index (Fir))
% 2: idx_Veg2 (veg 2 index (Larch))
% 3: idx_Veg3 (veg 3 index (Grass))
% 4: idx_Veg4 (veg 4 index (Shrub))
% 5: idx_Rock (rock index)
% 6: idx_Ice (ice index)
% 7: idx_Cleanice (clean-ice index)
% 8: idx_Debice (debris-covered ice index)
if output_manag(6) == 1

for ilc=1:length(LCinde) %%%loop over classes
    idx= LCinde{ilc}; %index for current class
    %%%%%%
    std_alp_soil_tg_LC(ilc)	=	std(alp_soil(idx));	%%[-]
    std_Ca_tg_LC(ilc)       =	std(Ca_S(idx));	%%[ppm]
    std_Cice_tg_LC(ilc)     =	std(Cice(idx));	%%[]
    std_Cicew_tg_LC(ilc)	=	std(Cicew(idx));	%%[]
    std_CK1_tg_LC(ilc)      =	std(CK1(idx));	%%[mm]
    std_Csno_tg_LC(ilc)     =	std(Csno(idx));	%%[]
    std_Csnow_tg_LC(ilc)	=	std(Csnow(idx));	%%[]
    std_DQ_S_tg_LC(ilc)     =	std(DQ_S(idx));	%%[]
    std_dQ_S_tg_LC(ilc)     =	std(dQ_S(idx));	%%[]
    std_Dr_H_tg_LC(ilc)     =	std(sum(Dr_H(idx),2));	%%[mm]
    std_Dr_L_tg_LC(ilc)     =	std(sum(Dr_L(idx),2));	%%[mm]
    std_Ds_tg_LC(ilc)       =	std(Ds_S(idx));	%%[°C]
    std_DT_S_tg_LC(ilc)     =	std(DT_S(idx));	%%[]
    std_dw_SNO_tg_LC(ilc)	=	std(dw_SNO(idx));	%%[]
    std_ea_tg_LC(ilc)       =	std(ea_S(idx));	%%[Pa]
    std_EG_tg_LC(ilc)       =	std(EG(idx));	%%[mm/h]
    std_EICE_tg_LC(ilc)     =	std(EICE(idx));	%%[mm/h]
    std_EIn_H_tg_LC(ilc)	=	std(sum(EIn_H(idx),2));	%%[mm/h]
    std_EIn_L_tg_LC(ilc)	=	std(sum(EIn_L(idx),2));	%%[mm/h]
    std_EIn_rock_tg_LC(ilc)	=	std(EIn_rock(idx));	%%[mm/h]
    std_EIn_urb_tg_LC(ilc)	=	std(EIn_urb(idx));	%%[mm/h]
    std_er_tg_LC(ilc)       =	std(er(idx));	%%[kg/s m^2]
    std_ESN_tg_LC(ilc)      =	std(ESN(idx)+ESN_In(idx));	%%[mm/h]
    std_SSN_tg_LC(ilc)      =	std(SSN(idx)+SSN_In(idx));	%%[mm/h]
    std_EWAT_tg_LC(ilc)     =	std(EWAT(idx));	%%[mm/h]
    std_f_tg_LC(ilc)        =	std(f(idx)*dth);	%%[mm]
    std_Fract_sat_tg_LC(ilc) =	sum(ZWT(idx)==0)/numel(ZWT(idx));	%%[-]
    std_FROCK_tg_LC(ilc)	=	std(FROCK(idx));	%%[mm]
    std_G_tg_LC(ilc)        =	std(G(idx));	%%[W/m^2]
    std_Gfin_tg_LC(ilc)     =	std(Gfin(idx));	%%[W/m^2]
    std_H_tg_LC(ilc)        =	std(H(idx));	%%[W/m^2]
    std_ICE_D_tg_LC(ilc)	=	std(ICE_D(idx));	%%[m]
    std_ICE_tg_LC(ilc)      =	std(ICE(idx));	%%[mm]
    std_Imelt_tg_LC(ilc)	=	std(Imelt(idx));	%%[mm]
    std_In_rock_tg_LC(ilc)	=	std(In_rock(idx));	%%[mm]
    std_In_SWE_tg_LC(ilc)	=	std(In_SWE(idx));	%%[mm]
    std_In_urb_tg_LC(ilc)	=	std(In_urb(idx));	%%[mm]
    std_IP_wc_tg_LC(ilc)	=	std(IP_wc(idx));	%%[mm]
    std_Lk_rock_tg_LC(ilc)	=	std(Lk_rock(idx)*dth);	%%[mm]
    std_Lk_tg_LC(ilc)       =	std(Lk(idx)*dth);	%%[mm]
    std_Lk_wat_tg_LC(ilc)	=	std(Lk_wat(idx)*dth);	%%[mm]
    std_N_tg_LC(ilc)        =	std(N_S(idx));	%%[-]
    std_NDVI_tg_LC(ilc)     =	std(NDVI(idx));	%%[-]
    std_NIce_tg_LC(ilc)     =	std(NIce(idx));	%%[mm]
    std_NIn_SWE_tg_LC(ilc)	=	std(NIn_SWE(idx));	%%[mm]
    std_O_tg_LC(ilc)        =	std(O_space(idx));	%%[-]
    std_OF_tg_LC(ilc)       =	std(OF(idx));	%%[]
    std_OS_tg_LC(ilc)       =	std(OS(idx));	%%[]
    std_PAR_tg_LC(ilc)      =	std(PAR_space(idx));%%[W/m^2]
    std_Pr_liq_tg_LC(ilc)	=	std(Pr_liq(idx));	%%[mm]
    std_Pr_sno_tg_LC(ilc)	=	std(Pr_sno(idx));	%%[mm]
    std_Pr_tg_LC(ilc)       =	std(Pr_S(idx));	%%[mm]
    std_Pre_tg_LC(ilc)      =	std(Pre_S(idx));	%%[mbar]
    std_Q_channel_tg_LC(ilc) =	std(Q_channel(idx));	%%[mm]
    std_q_runon_tg_LC(ilc)	=	std(q_runon(idx)*dth);	%%[mm]
    std_QE_tg_LC(ilc)       =	std(QE(idx));	%%[W/m^2]
    std_Qfm_tg_LC(ilc)      =	std(Qfm(idx));	%%[W/m^2]
    std_Qlat_in_tg_LC(ilc)	=	std(sum(Qi_in(idx)*dth,2));	%%[mm]
    std_Qlat_out_tg_LC(ilc)	=	std(sum(Qi_out(idx)*dth,2));	%%[mm]
    std_Qv_tg_LC(ilc)       =	std(Qv(idx));	%%[W/m^2]
    std_r_soil_tg_LC(ilc)	=	std(r_soil(idx));	%%[s/m]
    std_ra_tg_LC(ilc)       =	std(ra(idx));	%%[s/m]
    std_Rd_tg_LC(ilc)       =	std(Rd(idx));	%%[mm]
    std_Rh_tg_LC(ilc)       =	std(Rh(idx));	%%[mm]
    std_Rn_tg_LC(ilc)       =	std(Rn(idx));	%%[W/m^2]
    std_ros_tg_LC(ilc)      =	std(ros(idx));	%%[kg/m^3]
    std_Rsw_tg_LC(ilc)      =	std(Rsw_space(idx));	%%[W/m^2]
    std_SE_rock_tg_LC(ilc)	=	std(SE_rock(idx));	%%[]
    std_SE_urb_tg_LC(ilc)	=	std(SE_urb(idx));	%%[]
    std_Smelt_tg_LC(ilc)	=	std(Smelt(idx));	%%[mm]
    std_SND_tg_LC(ilc)      =	std(SND(idx));	%%[m]
    std_snow_albedo_tg_LC(ilc) = std(snow_albedo(idx)); %%[-]
    std_SP_wc_tg_LC(ilc)	=	std(SP_wc(idx));	%%[mm]
    std_SWE_tg_LC(ilc)      =	std(SWE(idx));	%%[mm]
    std_SWE_avalanched_tg_LC(ilc) =	std(SWE_avalanched(idx));	%%[mm]
    std_T_H_tg_LC(ilc)      =	std(sum(T_H(idx),2));	%%[mm/h]
    std_T_L_tg_LC(ilc)      =	std(sum(T_L(idx),2));	%%[mm/h]
    std_Ta_tg_LC(ilc)       =	std(Ta_S(idx));	%%[°C]
    std_Tdamp_tg_LC(ilc)	=	std(Tdamp(idx));	%%[°C]
    std_Tdew_tg_LC(ilc)     =	std(Tdew_S(idx));	%%[°C]
    std_Tdp_space_LC        =	std(Tdp,[],2);
    std_Tdp_tg_LC(ilc)      =	std(std_Tdp_space_LC(idx));	%%[°C]
    std_Tice_tg_LC(ilc)     =	std(Tice(idx));	%%[C]
    std_Ts_tg_LC(ilc)       =	std(Ts(idx));	%%[°C]
    std_TsVEG_tg_LC(ilc)	=	std(TsVEG(idx));	%%[°C]
    std_U_SWE_tg_LC(ilc)	=	std(U_SWE(idx));	%%[mm]
    std_WAT_tg_LC(ilc)      =	std(WAT(idx));	%%[mm]
    std_WIS_tg_LC(ilc)      =	std(WIS(idx));	%%[mm]
    std_WR_IP_tg_LC(ilc)	=	std(WR_IP(idx));	%%[mm]
    std_WR_SP_tg_LC(ilc)	=	std(WR_SP(idx));	%%[mm]
    std_Ws_tg_LC(ilc)       =	std(Ws_S(idx));	%%[m/s]
    std_ZWT_tg_LC(ilc)      =	std(ZWT(idx));	%%[mm]
    %%%
    std_EIn_tg_LC(ilc)      =	EIn_H_tg_LC(ilc)+ EIn_L_tg_LC(ilc);	%%[mm/h]
    std_Inveg_tg_LC(ilc)	=	std(sum(In_H(idx),2) +  sum(In_L(idx),2));
    std_O_space_LC          =	sum(V,2)./Zs_OUT + Ohy_OUT;	%%[-]
    std_PAR_space_LC        =	PARB_S + PARD_S;	%%[W/m^2]
    std_Rsw_space_LC        =	SAB1_S+ SAB2_S + SAD1_S+ SAD2_S;
    std_T_tg_LC(ilc)        =	T_L_tg_LC(ilc) +T_H_tg_LC(ilc);	%%[mm/h]
    std_V_tg_LC(ilc)        =	std(Asur(idx).*sum(V(idx),2));	%%[mm]
    std_Vice_tg_LC(ilc)     =	std(Asur(idx).*sum(Vice(idx),2));	%%[mm]
    std_In_tg_LC(ilc)       =	std(In_H_space(ilc) +  In_L_space(idx) +  SP_wc(idx) + In_SWE(idx) + In_urb(idx) + In_rock(idx)+ IP_wc(idx) );	%%[mm]
end

%%%%%%%%%%%%%% Land cover class std variable names %%%%%%%%%%%

vars_LC_std = {'Date','std_alp_soil_tg_LC','std_Ca_tg_LC','std_Cicew_tg_LC','std_Cice_tg_LC','std_CK1_tg_LC','std_Csnow_tg_LC',...
    'std_Csno_tg_LC','std_DQ_S_tg_LC','std_dQ_S_tg_LC','std_Dr_H_tg_LC','std_Dr_L_tg_LC','std_Ds_tg_LC','std_DT_S_tg_LC',...
    'std_dw_SNO_tg_LC','std_ea_tg_LC','std_EG_tg_LC','std_EICE_tg_LC','std_EIn_H_tg_LC','std_EIn_L_tg_LC','std_EIn_rock_tg_LC',...
    'std_EIn_tg_LC','std_EIn_urb_tg_LC','std_er_tg_LC','std_ESN_tg_LC','std_SSN_tg_LC','std_EWAT_tg_LC','std_Fract_sat_tg_LC',...
    'std_FROCK_tg_LC','std_f_tg_LC','std_Gfin_tg_LC','std_G_tg_LC','std_H_tg_LC','std_ICE_D_tg_LC','std_ICE_tg_LC','std_Imelt_tg_LC',...
    'std_Inveg_tg_LC','std_In_rock_tg_LC','std_In_SWE_tg_LC','std_In_urb_tg_LC','std_IP_wc_tg_LC','std_Lk_rock_tg_LC','std_Lk_tg_LC',...
    'std_Lk_wat_tg_LC','std_NIce_tg_LC','std_NIn_SWE_tg_LC','std_N_tg_LC','std_NDVI_tg_LC',...
    'std_OF_tg_LC','std_OS_tg_LC','std_O_space','std_O_tg_LC','std_PAR_space','std_PAR_tg_LC','std_Pre_tg_LC','std_Pr_liq_tg_LC',...
    'std_Pr_sno_tg_LC','std_Pr_tg_LC','std_QE_tg_LC','std_Qfm_tg_LC','std_Qlat_in_tg_LC','std_Qlat_out_tg_LC','std_Qv_tg_LC','std_Q_channel_tg_LC',...
    'std_q_runon_tg_LC','std_ra_tg_LC','std_Rd_tg_LC','std_Rh_tg_LC','std_Rn_tg_LC','std_ros_tg_LC','std_Rsw_space','std_Rsw_tg_LC',...
    'std_r_soil_tg_LC','std_SE_rock_tg_LC','std_SE_urb_tg_LC','std_Smelt_tg_LC','std_std_SND_tg_LC','std_snow_albedo_tg_LC','std_SP_wc_tg_LC',...
    'std_SWE_avalanched_tg_LC','std_SWE_tg_LC','std_Ta_tg_LC','std_Tdamp_tg_LC','std_Tdew_tg_LC','std_Tdp_space','std_Tdp_tg_LC','std_Tice_tg_LC',...
    'std_TsVEG_tg_LC','std_Ts_tg_LC','std_T_H_tg_LC','std_T_L_tg_LC','std_T_tg_LC','std_U_SWE_tg_LC','std_Vice_tg_LC','std_V_tg_LC','std_WAT_tg_LC',...
    'std_WIS_tg_LC','std_WR_IP_tg_LC','std_WR_SP_tg_LC','std_Ws_tg_LC','std_ZWT_tg_LC','std_In_tg_LC'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [Output: columns in alphabetical order] %%%
if t==2
    for ilc=1:length(LCinde) %%%loop over classes
        tit9{ilc}=strcat(outlocation,'OUTPUT_',SITE,'_STD_LC_',LC_names{ilc},'.dat');
        fid9(ilc)=fopen(tit9{ilc},'a');
        % Add labels to column list
           for ii = 1:length(vars_LC_std)-1
               fprintf(fid9(ilc),'%s\t',vars_LC_std{ii});
           end 
           fprintf(fid9(ilc),'%s\t\n',vars_LC_std{length(vars_LC_std)}); 
    end
end

if t==t1_reinit
    for ilc=1:length(LCinde) %%%loop over classes

   tit9{ilc}=strcat(outlocation,'OUTPUT_',SITE,'_STD_LC_',LC_names{ilc},'.dat');

   % Check the line number after which delete everything (based on datestamp)
   fid9(ilc)=fopen(tit9{ilc},'r+');  % Open the file
   tline = fgetl(fid9(ilc)); %Get the first line
   lineCounter = 1;
   Date_str_m1 = datestr(datetime(Date_str) - hours(1)); %Get the datestamp of the last line to keep
    while ischar(tline)
       if length(tline) > 19 && strcmp(tline(1:20), Date_str_m1)
         break;
       end
      % Read next line
      tline = fgetl(fid9(ilc));
      lineCounter = lineCounter + 1;
    end

   % Copy all the lines until the line number found in the previous section
     fid9(ilc) = fopen(tit9{ilc},'r');
     your_text = cell(lineCounter,1);
     for ii = 1:lineCounter
       your_text(ii) = {fgetl(fid9(ilc))}; 
     end
     fclose(fid9(ilc));

    % Write all the lines until the line number found in the previous section
      fid9(ilc) = fopen(tit9{ilc},'w');
        for ii = 1:lineCounter
          fprintf(fid9(ilc),'%s\n',your_text{ii});
        end
      clear your_text               

    end
end
for ilc=1:length(LCinde) %%%loop over classes
    %%% START <<OUTPUT_ZZZ_STD_LC_code_YYY.dat>> %%%
    fprintf(fid9(ilc),'%s\t',Date_str);
    fprintf(fid9(ilc),'%g\t',std_alp_soil_tg_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_Ca_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_Cicew_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_Cice_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_CK1_tg_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_Csnow_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_Csno_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_DQ_S_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_dQ_S_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_Dr_H_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_Dr_L_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_Ds_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_DT_S_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_dw_SNO_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_ea_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_EG_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_EICE_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_EIn_H_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_EIn_L_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_EIn_rock_tg_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_EIn_tg_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_EIn_urb_tg_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_er_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_ESN_tg_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_SSN_tg_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_EWAT_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_Fract_sat_tg_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_FROCK_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_f_tg_LC(ilc));        
    fprintf(fid9(ilc),'%g\t',std_Gfin_tg_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_G_tg_LC(ilc));        
    fprintf(fid9(ilc),'%g\t',std_H_tg_LC(ilc));        
    fprintf(fid9(ilc),'%g\t',std_ICE_D_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_ICE_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_Imelt_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_Inveg_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_In_rock_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_In_SWE_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_In_urb_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_IP_wc_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_Lk_rock_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_Lk_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_Lk_wat_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_NIce_tg_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_NIn_SWE_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_N_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_NDVI_tg_LC(ilc));        
    fprintf(fid9(ilc),'%g\t',std_OF_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_OS_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_O_space_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_O_tg_LC(ilc));       	
    fprintf(fid9(ilc),'%g\t',std_PAR_space_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_PAR_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_Pre_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_Pr_liq_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_Pr_sno_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_Pr_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_QE_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_Qfm_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_Qlat_in_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_Qlat_out_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_Qv_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_Q_channel_tg_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_q_runon_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_ra_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_Rd_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_Rh_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_Rn_tg_LC(ilc));       
    fprintf(fid9(ilc),'%g\t',std_ros_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_Rsw_space_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_Rsw_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_r_soil_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_SE_rock_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_SE_urb_tg_LC(ilc));	
    fprintf(fid9(ilc),'%g\t',std_Smelt_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_SND_tg_LC(ilc));      
    fprintf(fid9(ilc),'%g\t',std_snow_albedo_tg_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_SP_wc_tg_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_SWE_avalanched_tg_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_SWE_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_Ta_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_Tdamp_tg_LC(ilc)); 
    fprintf(fid9(ilc),'%g\t',std_Tdew_tg_LC(ilc));  
    fprintf(fid9(ilc),'%g\t',std_Tdp_space_LC(ilc));
    fprintf(fid9(ilc),'%g\t',std_Tdp_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_Tice_tg_LC(ilc));  
    fprintf(fid9(ilc),'%g\t',std_TsVEG_tg_LC(ilc)); 
    fprintf(fid9(ilc),'%g\t',std_Ts_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_T_H_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_T_L_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_T_tg_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_U_SWE_tg_LC(ilc)); 
    fprintf(fid9(ilc),'%g\t',std_Vice_tg_LC(ilc));  
    fprintf(fid9(ilc),'%g\t',std_V_tg_LC(ilc));     
    fprintf(fid9(ilc),'%g\t',std_WAT_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_WIS_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t',std_WR_IP_tg_LC(ilc)); 
    fprintf(fid9(ilc),'%g\t',std_WR_SP_tg_LC(ilc)); 
    fprintf(fid9(ilc),'%g\t',std_Ws_tg_LC(ilc));    
    fprintf(fid9(ilc),'%g\t',std_ZWT_tg_LC(ilc));   
    fprintf(fid9(ilc),'%g\t\n',std_In_tg_LC(ilc)); 
    %%% END <<OUTPUT_ZZZ_STD_LC_code_YYY.dat>> %%%	
end

%%%
if t==N_time_step
    for ilc=1:length(LCinde)
        fclose(fid9(ilc));
    end
end

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TEMPORAL AVERAGE MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if output_manag(7) == 1

if t==2
    toutp=2;
else
    toutp=toutp+1;
end
%%%
if toutp==2
    %%%%%
    An_H_spatial	=	An_H;
    An_L_spatial	=	An_L;
    ANPP_H_spatial	=	ANPP_H;
    ANPP_L_spatial	=	ANPP_L;
    B_H_spatial     =	B_H;
    B_L_spatial     =	B_L;
    Ca_spatial      =	Ca_S;
    Ci_shdH_spatial	=	Ci_shdH;
    Ci_shdL_spatial	=	Ci_shdL;
    Ci_sunH_spatial	=	Ci_sunH;
    Ci_sunL_spatial	=	Ci_sunL;
    Cice_spatial	=	Cice;
    Cicew_spatial	=	Cicew;
    CK1_spatial     =	CK1;
    Csno_spatial	=	Csno;
    Csnow_spatial	=	Csnow;
    DQ_S_spatial	=	DQ_S;
    dQ_S_spatial	=	dQ_S;
    Dr_H_spatial	=	Dr_H_space;
    Dr_L_spatial	=	Dr_L_space;
    Ds_spatial      =	Ds_S;
    DT_S_spatial	=	DT_S;
    dw_SNO_spatial	=	dw_SNO;
    ea_spatial      =	ea_S;
    EG_spatial      =	EG;
    EICE_spatial	=	EICE;
    EIn_H_spatial	=	EIn_H_space;
    EIn_L_spatial	=	EIn_L_space;
    EIn_rock_spatial =	EIn_rock;
    er_spatial      =	er;
    ESN_spatial     =	ESN + ESN_In;
    SSN_spatial     =	SSN + SSN_In;
    EWAT_spatial	=	EWAT;
    f_spatial       =	f;
    FROCK_spatial	=	FROCK;
    G_spatial       =	G;
    Gfin_spatial	=	Gfin;
    GPP_H_spatial	=	NPP_H+RA_H;
    GPP_L_spatial	=	NPP_L +RA_L;
    H_spatial       =	H;
    hc_H_spatial	=	hc_H;
    hc_L_spatial	=	hc_L;
    ICE_D_spatial	=	ICE_D;
    ICE_spatial     =	ICE;
    Imelt_spatial	=	Imelt;
    In_rock_spatial	=	In_rock;
    In_spatial      =	In_H_space + In_L_space + SP_wc + In_SWE  + In_urb + In_rock + IP_wc;
    In_SWE_spatial	=	In_SWE;
    Inveg_spatial	=	In_H_space + In_L_space;
    IP_wc_spatial	=	IP_wc;
    LAI_H_spatial	=	LAI_H;
    LAI_L_spatial	=	LAI_L;
    LAIdead_H_spatial =	LAIdead_H;
    LAIdead_L_spatial =	LAIdead_L;
    Lk_rock_spatial	=	Lk_rock;
    Lk_spatial      =	Lk;
    Lk_wat_spatial	=	Lk_wat;
    N_spatial       =	N_S;
    NDVI_spatial       =	NDVI;
    NICe_spatial	=	NIce;
    NPP_H_spatial	=	NPP_H;
    NPP_L_spatial	=	NPP_L;
    O_spatial       =	O_space;
    OF_spatial      =	OF;
    OH_spatial      =	OH;
    OL_spatial      =	OL;
    OS_spatial      =	OS;
    PAR_spatial     =	PAR_space;
    Pr_liq_spatial	=	Pr_liq;
    Pr_sno_spatial	=	Pr_sno;
    Pr_spatial      =	Pr_S;
    Pre_spatial     =	Pre_S;
    Q_channel_spatial =	Q_channel;
    q_runon_spatial	=	q_runon;
    QE_spatial      =	QE;
    Qfm_spatial     =	Qfm;
    Qlat_in_spatial	=	Qi_in_space;
    Qlat_out_spatial =	Qi_out_space;
    Qv_spatial      =	Qv;
    RA_H_spatial	=	RA_H;
    RA_L_spatial	=	RA_L;
    Rd_spatial      =	Rd;
    Rdark_H_spatial	=	Rdark_H;
    Rdark_L_spatial	=	Rdark_L;
    Rg_H_spatial	=	Rg_H;
    Rg_L_spatial	=	Rg_L;
    Rh_spatial      =	Rh;
    Rmc_H_spatial	=	Rmc_H;
    Rmc_L_spatial	=	Rmc_L;
    Rmr_H_spatial	=	Rmr_H;
    Rmr_L_spatial	=	Rmr_L;
    Rms_H_spatial	=	Rms_H;
    Rms_L_spatial	=	Rms_L;
    Rn_spatial      =	Rn;
    ros_spatial     =	ros;
    Rsw_spatial     =	Rsw_space;
    SAI_H_spatial	=	SAI_H;
    SAI_L_spatial	=	SAI_L;
    SAT_spatial     =	(ZWT==0);
    Smelt_spatial	=	Smelt;
    SND_spatial     =	SND;
    snow_albedo_spatial = snow_albedo;
    SP_wc_spatial	=	SP_wc;
    SWE_spatial     =	SWE;
    SWE_avalanched_spatial =	SWE_avalanched;
    T_H_spatial     =	T_H_space;
    T_L_spatial     =	T_L_space;
    Ta_spatial      =	Ta_S;
    Tdamp_spatial	=	Tdamp;
    Tdew_spatial	=	Tdew_S;
    Tdp_H_spatial	=	Tdp_H;
    Tdp_L_spatial	=	Tdp_L;
    Tdp_spatial     =	Tdp_space;
    Tice_spatial	=	Tice;
    Ts_spatial      =	Ts;
    TsVEG_spatial	=	TsVEG;
    U_SWE_spatial	=	U_SWE;
    V_spatial       =	V_space;
    Vice_spatial	=	Vice_space;
    WAT_spatial     =	WAT;
    WIS_spatial     =	WIS;
    WR_IP_spatial	=	WR_IP;
    WR_SP_spatial	=	WR_SP;
    Ws_spatial      =	Ws_S;
    ZWT_spatial     =	ZWT;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %%%%%
    An_H_spatial        = ((toutp-2)*An_H_spatial + An_H)/(toutp-1) ;
    An_L_spatial        = ((toutp-2)*An_L_spatial + An_L)/(toutp-1) ;
    ANPP_H_spatial      = ((toutp-2)*ANPP_H_spatial + ANPP_H)/(toutp-1) ;
    ANPP_L_spatial      = ((toutp-2)*ANPP_L_spatial + ANPP_L)/(toutp-1) ;
    B_H_spatial         = ((toutp-2)*B_H_spatial + B_H)/(toutp-1) ;
    B_L_spatial         = ((toutp-2)*B_L_spatial + B_L)/(toutp-1) ;
    Ca_spatial          = ((toutp-2)*Ca_spatial + Ca_S)/(toutp-1) ;
    Ci_shdH_spatial     = ((toutp-2)*Ci_shdH_spatial + Ci_shdH)/(toutp-1) ;
    Ci_shdL_spatial     = ((toutp-2)*Ci_shdL_spatial + Ci_shdL)/(toutp-1) ;
    Ci_sunH_spatial     = ((toutp-2)*Ci_sunH_spatial + Ci_sunH)/(toutp-1) ;
    Ci_sunL_spatial     = ((toutp-2)*Ci_sunL_spatial + Ci_sunL)/(toutp-1) ;
    Cice_spatial        = Cice_spatial + Cice ;
    Cicew_spatial       = Cicew_spatial + Cicew ;
    CK1_spatial         = CK1_spatial + CK1;
    Csno_spatial        = Csno_spatial + Csno ;
    Csnow_spatial       = Csnow_spatial + Csnow ;
    DQ_S_spatial        = DQ_S_spatial + DQ_S;
    dQ_S_spatial        = dQ_S_spatial + dQ_S;
    Dr_H_spatial        = ((toutp-2)*Dr_H_spatial + Dr_H_space)/(toutp-1) ;
    Dr_L_spatial        = ((toutp-2)*Dr_L_spatial + Dr_L_space)/(toutp-1) ;
    Ds_spatial          = ((toutp-2)*Ds_spatial + Ds_S)/(toutp-1) ;
    DT_S_spatial        = DT_S_spatial + DT_S;
    dw_SNO_spatial      = ((toutp-2)*dw_SNO_spatial + dw_SNO)/(toutp-1) ;
    ea_spatial          = ((toutp-2)*ea_spatial + ea_S)/(toutp-1) ;
    EG_spatial          = ((toutp-2)*EG_spatial + EG)/(toutp-1) ;
    EICE_spatial        = ((toutp-2)*EICE_spatial + EICE)/(toutp-1) ;
    EIn_H_spatial       = ((toutp-2)*EIn_H_spatial + EIn_H_space)/(toutp-1) ;
    EIn_L_spatial       = ((toutp-2)*EIn_L_spatial + EIn_L_space)/(toutp-1) ;
    EIn_rock_spatial    = ((toutp-2)*EIn_rock_spatial + EIn_rock)/(toutp-1) ;
    er_spatial          = ((toutp-2)*er_spatial + er)/(toutp-1) ;
    ESN_spatial         = ((toutp-2)*ESN_spatial + ESN + ESN_In)/(toutp-1) ;
    SSN_spatial         = ((toutp-2)*SSN_spatial + SSN + SSN_In)/(toutp-1) ;
    EWAT_spatial        = ((toutp-2)*EWAT_spatial + EWAT)/(toutp-1) ;
    f_spatial           = ((toutp-2)*f_spatial +  f)/(toutp-1) ;
    FROCK_spatial       = ((toutp-2)*FROCK_spatial + FROCK)/(toutp-1) ;
    G_spatial           = ((toutp-2)*G_spatial + G)/(toutp-1) ;
    Gfin_spatial        = ((toutp-2)*Gfin_spatial + Gfin)/(toutp-1) ;
    GPP_H_spatial       = ((toutp-2)*GPP_H_spatial + NPP_H+RA_H)/(toutp-1) ;
    GPP_L_spatial       = ((toutp-2)*GPP_L_spatial + NPP_L+RA_L)/(toutp-1) ;
    H_spatial           = ((toutp-2)*H_spatial + H)/(toutp-1) ;
    hc_H_spatial        = ((toutp-2)*hc_H_spatial +  hc_H)/(toutp-1) ;
    hc_L_spatial        = ((toutp-2)*hc_L_spatial +  hc_L)/(toutp-1) ;
    ICE_D_spatial       = ((toutp-2)*ICE_D_spatial + ICE_D)/(toutp-1) ;
    ICE_spatial         = ((toutp-2)*ICE_spatial + ICE)/(toutp-1) ;
    Imelt_spatial       = ((toutp-2)*Imelt_spatial + Imelt)/(toutp-1) ;
    In_rock_spatial     = ((toutp-2)*In_rock_spatial + In_rock)/(toutp-1) ;
    In_spatial          = ((toutp-2)*In_spatial + In_H_space + In_L_space + SP_wc + In_SWE + In_urb + In_rock + IP_wc)/(toutp-1) ;
    In_SWE_spatial      = ((toutp-2)*In_SWE_spatial + In_SWE)/(toutp-1) ;
    Inveg_spatial       = ((toutp-2)*Inveg_spatial + In_H_space + In_L_space)/(toutp-1) ;
    IP_wc_spatial       = ((toutp-2)*IP_wc_spatial + IP_wc)/(toutp-1) ;
    LAI_H_spatial       = ((toutp-2)*LAI_H_spatial + LAI_H)/(toutp-1) ;
    LAI_L_spatial       = ((toutp-2)*LAI_L_spatial + LAI_L)/(toutp-1) ;
    LAIdead_H_spatial   = ((toutp-2)*LAIdead_H_spatial +  LAIdead_H)/(toutp-1) ;
    LAIdead_L_spatial   = ((toutp-2)*LAIdead_L_spatial +  LAIdead_L)/(toutp-1) ;
    Lk_rock_spatial     = ((toutp-2)*Lk_rock_spatial +  Lk_rock)/(toutp-1) ;
    Lk_spatial          = ((toutp-2)*Lk_spatial +  Lk)/(toutp-1) ;
    Lk_wat_spatial      = ((toutp-2)*Lk_wat_spatial +  Lk_wat)/(toutp-1) ;
    N_spatial           = ((toutp-2)*N_spatial + N_S)/(toutp-1) ;
    NDVI_spatial        = ((toutp-2)*NDVI_spatial + NDVI)/(toutp-1) ;
    NICe_spatial        = ((toutp-2)*NICe_spatial + NIce)/(toutp-1) ;
    NPP_H_spatial       = ((toutp-2)*NPP_H_spatial + NPP_H)/(toutp-1) ;
    NPP_L_spatial       = ((toutp-2)*NPP_L_spatial + NPP_L)/(toutp-1) ;
    O_spatial           = ((toutp-2)*O_spatial + O_space)/(toutp-1) ;
    OF_spatial          = ((toutp-2)*OF_spatial +  OF)/(toutp-1) ;
    OH_spatial          = ((toutp-2)*OH_spatial + OH)/(toutp-1) ;
    OL_spatial          = ((toutp-2)*OL_spatial + OL)/(toutp-1) ;
    OS_spatial          = ((toutp-2)*OS_spatial +  OS)/(toutp-1) ;
    PAR_spatial         = ((toutp-2)*PAR_spatial + PAR_space)/(toutp-1) ;
    Pr_liq_spatial      = Pr_liq_spatial + Pr_liq;
    Pr_sno_spatial      = Pr_sno_spatial + Pr_sno;
    Pr_spatial          = Pr_spatial + Pr_S;
    Pre_spatial         = ((toutp-2)*Pre_spatial + Pre_S)/(toutp-1) ;
    Q_channel_spatial   = ((toutp-2)*Q_channel_spatial +  Q_channel)/(toutp-1) ;
    q_runon_spatial     = ((toutp-2)*q_runon_spatial +  q_runon)/(toutp-1) ;
    QE_spatial          = ((toutp-2)*QE_spatial + QE)/(toutp-1) ;
    Qfm_spatial         = ((toutp-2)*Qfm_spatial + Qfm)/(toutp-1) ;
    Qlat_in_spatial     = ((toutp-2)*Qlat_in_spatial +  Qi_in_space)/(toutp-1) ;
    Qlat_out_spatial    = ((toutp-2)*Qlat_out_spatial +  Qi_out_space)/(toutp-1) ;
    Qv_spatial          = ((toutp-2)*Qv_spatial + Qv)/(toutp-1) ;
    RA_H_spatial        = ((toutp-2)*RA_H_spatial + RA_H)/(toutp-1) ;
    RA_L_spatial        = ((toutp-2)*RA_L_spatial + RA_L)/(toutp-1) ;
    Rd_spatial          = ((toutp-2)*Rd_spatial +  Rd)/(toutp-1) ;
    Rdark_H_spatial     = ((toutp-2)*Rdark_H_spatial + Rdark_H)/(toutp-1) ;
    Rdark_L_spatial     = ((toutp-2)*Rdark_L_spatial + Rdark_L)/(toutp-1) ;
    Rg_H_spatial        = ((toutp-2)*Rg_H_spatial + Rg_H)/(toutp-1) ;
    Rg_L_spatial        = ((toutp-2)*Rg_L_spatial + Rg_L)/(toutp-1) ;
    Rh_spatial          = ((toutp-2)*Rh_spatial +  Rh)/(toutp-1) ;
    Rmc_H_spatial       = ((toutp-2)*Rmc_H_spatial + Rmc_H)/(toutp-1) ;
    Rmc_L_spatial       = ((toutp-2)*Rmc_L_spatial + Rmc_L)/(toutp-1) ;
    Rmr_H_spatial       = ((toutp-2)*Rmr_H_spatial + Rmr_H)/(toutp-1) ;
    Rmr_L_spatial       = ((toutp-2)*Rmr_L_spatial + Rmr_L)/(toutp-1) ;
    Rms_H_spatial       = ((toutp-2)*Rms_H_spatial + Rms_H)/(toutp-1) ;
    Rms_L_spatial       = ((toutp-2)*Rms_L_spatial + Rms_L)/(toutp-1) ;
    Rn_spatial          = ((toutp-2)*Rn_spatial + Rn)/(toutp-1) ;
    ros_spatial         = ((toutp-2)*ros_spatial + ros)/(toutp-1) ;
    Rsw_spatial         = ((toutp-2)*Rsw_spatial + Rsw_space)/(toutp-1) ;
    SAI_H_spatial       = ((toutp-2)*SAI_H_spatial +  SAI_H)/(toutp-1) ;
    SAI_L_spatial       = ((toutp-2)*SAI_L_spatial +  SAI_L)/(toutp-1) ;
    SAT_spatial         = SAT_spatial + (ZWT==0) ;
    Smelt_spatial       = ((toutp-2)*Smelt_spatial + Smelt)/(toutp-1) ;
    SND_spatial         = ((toutp-2)*SND_spatial + SND)/(toutp-1) ;
    snow_albedo_spatial = ((toutp-2)*snow_albedo_spatial + snow_albedo_out)/(toutp-1) ;
    SP_wc_spatial       = ((toutp-2)*SP_wc_spatial + SP_wc)/(toutp-1) ;
    SWE_spatial         = ((toutp-2)*SWE_spatial + SWE)/(toutp-1) ;
    SWE_avalanched_spatial    = ((toutp-2)*SWE_avalanched_spatial + SWE_avalanched)/(toutp-1) ;
    T_H_spatial         = ((toutp-2)*T_H_spatial + T_H_space)/(toutp-1) ;
    T_L_spatial         = ((toutp-2)*T_L_spatial + T_L_space)/(toutp-1) ;
    Ta_spatial          = ((toutp-2)*Ta_spatial + Ta_S)/(toutp-1) ;
    Tdamp_spatial       = ((toutp-2)*Tdamp_spatial + Tdamp)/(toutp-1) ;
    Tdew_spatial        = ((toutp-2)*Tdew_spatial + Tdew_S)/(toutp-1) ;
    Tdp_spatial         = ((toutp-2)*Tdp_spatial + Tdp_space)/(toutp-1) ;
    Tice_spatial        = ((toutp-2)*Tice_spatial + Tice)/(toutp-1) ;
    Ts_spatial          = ((toutp-2)*Ts_spatial + Ts)/(toutp-1) ;
    TsVEG_spatial       = ((toutp-2)*TsVEG_spatial + TsVEG)/(toutp-1) ;
    U_SWE_spatial       = ((toutp-2)*U_SWE_spatial + U_SWE)/(toutp-1) ;
    V_spatial           = ((toutp-2)*V_spatial + V_space)/(toutp-1) ;
    Vice_spatial        = ((toutp-2)*Vice_spatial + Vice_space)/(toutp-1) ;
    WAT_spatial         = ((toutp-2)*WAT_spatial + WAT)/(toutp-1) ;
    WIS_spatial         = ((toutp-2)*WIS_spatial +  WIS)/(toutp-1) ;
    WR_IP_spatial       = ((toutp-2)*WR_IP_spatial + WR_IP)/(toutp-1) ;
    WR_SP_spatial       = ((toutp-2)*WR_SP_spatial + WR_SP)/(toutp-1) ;
    Ws_spatial          = ((toutp-2)*Ws_spatial + Ws_S)/(toutp-1) ;
    ZWT_spatial         = ((toutp-2)*ZWT_spatial +  ZWT)/(toutp-1) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [Output: columns in alphabetical order] %%%
if  length(intersect(t,tstore))==1 ||  t==N_time_step
    Title_save = strcat(outlocation,'OUTPUT_',SITE,'_SPATIAL_',num2str(t));
	%%% START <<OUTPUT_ZZZ_SPATIAL_YYY.dat>> %%%
    save(Title_save,...
        'ANPP_H_spatial',...
        'ANPP_L_spatial',...
        'An_H_spatial',...
        'An_L_spatial',...
        'B_H_spatial',...
        'B_L_spatial',...
        'Ca_spatial',...
        'cellsize',...
        'Cicew_spatial',...
        'Cice_spatial',...
        'Ci_shdH_spatial',...
        'Ci_shdL_spatial',...
        'Ci_sunH_spatial',...
        'Ci_sunL_spatial',...
        'CK1_spatial',...
        'Csnow_spatial',...
        'Csno_spatial',...
        'dQ_S_spatial',...
        'DQ_S_spatial',...
        'Dr_H_spatial',...
        'Dr_L_spatial',...
        'Ds_spatial',...
        'DTM',...
        'DT_S_spatial',...
        'dw_SNO_spatial',...
        'ea_spatial',...
        'EG_spatial',...
        'EICE_spatial',...
        'EIn_H_spatial',...
        'EIn_L_spatial',...
        'EIn_rock_spatial',...
        'er_spatial',...
        'ESN_spatial',...
        'SSN_spatial',...
        'EWAT_spatial',...
        'FROCK_spatial',...
        'f_spatial',...
        'Gfin_spatial',...
        'GPP_H_spatial',...
        'GPP_L_spatial',...
        'G_spatial',...
        'hc_H_spatial',...
        'hc_L_spatial',...
        'H_spatial',...
        'ICE_D_spatial',...
        'ICE_spatial',...
        'Imelt_spatial',...
        'Inveg_spatial',...
        'In_rock_spatial',...
        'In_spatial',...
        'In_SWE_spatial',...
        'IP_wc_spatial',...
        'LAIdead_H_spatial',...
        'LAIdead_L_spatial',...
        'LAI_H_spatial',...
        'LAI_L_spatial',...
        'Lk_rock_spatial',...
        'Lk_spatial',...
        'Lk_wat_spatial',...
        'm_cell',...
        'NICe_spatial',...
        'NPP_H_spatial',...
        'NPP_L_spatial',...
        'n_cell',...
        'N_spatial',...
        'NDVI_spatial',...
        'OF_spatial',...
        'OH_spatial',...
        'OL_spatial',...
        'OS_spatial',...
        'O_spatial',...
        'PAR_spatial',...
        'Pre_spatial',...
        'Pr_liq_spatial',...
        'Pr_sno_spatial',...
        'Pr_spatial',...
        'QE_spatial',...
        'Qfm_spatial',...
        'Qlat_in_spatial',...
        'Qlat_out_spatial',...
        'Qv_spatial',...
        'Q_channel_spatial',...
        'q_runon_spatial',...
        'RA_H_spatial',...
        'RA_L_spatial',...
        'Rdark_H_spatial',...
        'Rdark_L_spatial',...
        'Rd_spatial',...
        'Rg_H_spatial',...
        'Rg_L_spatial',...
        'Rh_spatial',...
        'Rmc_H_spatial',...
        'Rmc_L_spatial',...
        'Rmr_H_spatial',...
        'Rmr_L_spatial',...
        'Rms_H_spatial',...
        'Rms_L_spatial',...
        'Rn_spatial',...
        'ros_spatial',...
        'Rsw_spatial',...
        'SAI_H_spatial',...
        'SAI_L_spatial',...
        'SAT_spatial',...
        'Smelt_spatial',...
        'SND_spatial',...
        'snow_albedo_spatial',...
        'SP_wc_spatial',...
        'SWE_avalanched_spatial',...
        'SWE_spatial',...
        'Ta_spatial',...
        'Tdamp_spatial',...
        'Tdew_spatial',...
        'Tdp_spatial',...
        'Tice_spatial',...
        'TsVEG_spatial',...
        'Ts_spatial',...
        'T_H_spatial',...
        'T_L_spatial',...
        'U_SWE_spatial',...
        'Vice_spatial',...
        'V_spatial',...
        'WAT_spatial',...
        'WIS_spatial',...
        'WR_IP_spatial',...
        'WR_SP_spatial',...
        'Ws_spatial',...
        'xllcorner',...
        'x_cell',...
        'yllcorner',...
        'y_cell',...
        'ZWT_spatial');
		%%% END <<OUTPUT_ZZZ_SPATIAL_YYY.dat>> %%%
    %%%%%%%%%%
    toutp = 1;
    %%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% SND OUTPUT DAILY MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if t==2
    toutp_snow=2;
else
    toutp_snow=toutp_snow+1;
end
if toutp_snow==2
	SND_spatial_daily     =	SND;
    ros_spatial_daily     =	ros;
    SSN_spatial_daily     = SSN;
    Smelt_spatial_daily   = Smelt;
    Ice_spatial_daily   = ICE;

else
    SND_spatial_daily         = ((toutp_snow-2)*SND_spatial_daily + SND)/(toutp_snow-1) ;
    ros_spatial_daily         = ((toutp_snow-2)*ros_spatial_daily + ros)/(toutp_snow-1) ;
    SSN_spatial_daily         = ((toutp_snow-2)*SSN_spatial_daily + SSN)/(toutp_snow-1) ;
    Smelt_spatial_daily         = ((toutp_snow-2)*Smelt_spatial_daily + Smelt)/(toutp_snow-1) ;
    Ice_spatial_daily         = ((toutp_snow-2)*Ice_spatial_daily + ICE)/(toutp_snow-1) ;

end
	
if  length(intersect(t,tstore_snow))==1 ||  t==N_time_step
	Title_save = strcat(outlocation,'OUTPUT_',TITLE_SAVE,'_SNOWMAP_',num2str(t));
    %%% START <<OUTPUT_ZZZ_SNOWMAP_YYY.dat>> %%%
	save(Title_save,...
        'SND_spatial_daily',...
        'ros_spatial_daily',...
        'SSN_spatial_daily',...
        'Smelt_spatial_daily',...
        'Ice_spatial_daily',...
        'DTM',...
        'x_cell',...
        'y_cell');
    %%% END <<OUTPUT_ZZZ_SNOWMAP_YYY.dat>> %%%

	toutp_snow = 1;
		
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ALB OUTPUT DAILY MAPS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    toutp_alb=2;
else
    toutp_alb=toutp_alb+1;
end
if toutp_alb==2
	snoalb_biweek_spatial     =	snow_albedo_out;
    suralb_biweek_spatial     = surface_albedo_out;
else
    snoalb_biweek_spatial         = ((toutp_alb-2)*snoalb_biweek_spatial + snow_albedo_out)/(toutp_alb-1) ;
    suralb_biweek_spatial         = ((toutp_alb-2)*suralb_biweek_spatial + surface_albedo_out)/(toutp_alb-1) ;

end
	
if  length(intersect(t,tstore_alb))==1 ||  t==N_time_step
	
	Title_save = strcat(outlocation,'OUTPUT_',TITLE_SAVE,'_SNOALB_',num2str(t));
	%%% START <<OUTPUT_ZZZ_ALBEDO_YYY.dat>> %%%
    save(Title_save,...
        'snoalb_biweek_spatial',...
        'suralb_biweek_spatial',...
        'DTM',...
        'x_cell',...
        'y_cell');
    %%% END <<OUTPUT_ZZZ_ALBEDO_YYY.dat>> %%%    
	toutp_alb = 1;
		
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TRACKED PIXELS TIME SERIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if output_manag(8) == 1

%%%%%%%% Variable names for pixel time series %%%%%%%

vars_pix = {'Date','Asur','alp_soil','Ca_S','Cice','Cicew','CK1','Csno','Csnow','cos_fst','Ct','DEB','DQ_S','dQ_S','Ds_S','DT_S',...
    'dw_SNO','e_sno','ea_S','EG','EICE','EIn_rock','EIn_urb','er','ESN','SSN','ESN_In','SSN_In','EWAT','f','FROCK','G',...
    'Gfin','H','ICE','ICE_D','Imelt','In_rock','In_SWE','In_urb','IP_wc','Lk','Lk_rock','Lk_wat','NIce','NIn_SWE','N_S','NDVI',...
    'OF','OS','PAR_space','PARB_S','PARD_S','Pre_S','Pr_liq','Pr_S','Pr_sno','QE','Qfm','QpointC','QpointH','Qv','Q_channel',...
    'q_runon','ra','Rd','Rh','Rn','ros','Rsw_space','r_soil','SAB1_S','SAB2_S','SAD1_S','SAD2_S','SE_rock','SE_urb','Slo_top_S','Smelt','SND','snow_albedo',...
    'SP_wc','SSN','surface_albedo','SWE','SWE_avalanched','Ta_S','Tdamp','Tdew_S','Tice','Ts','U_S','UpointC','UpointH','U_SWE','WAT','WIS','WR_IP','WR_SP','Ws_S','ZWT'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% [Output: columns in alphabetical order] %%%
for ipo=1:npoint
    ij = sub2ind(size(DTM),Yout(ipo),Xout(ipo));
    if t==2
        tit5{1,ipo}=strcat(outlocation,'OUTPUT_',SITE,'_PIXEL_',strrep(POI_names(ipo)," ",""),'.dat');
        fid5(1,ipo)=fopen(tit5{1,ipo},'a');
               
        % Add labels to column list
        for ii = 1:length(vars_pix)-1
           fprintf(fid5(1,ipo),'%s\t',vars_pix{ii});
        end 
        fprintf(fid5(1,ipo),'%s\t\n',vars_pix{length(vars_pix)});
    end

   if t==t1_reinit

   tit5{1,ipo}=strcat(outlocation,'OUTPUT_',SITE,'_PIXEL_',strrep(POI_names(ipo)," ",""),'.dat');

   % Check the line number after which delete everything (based on datestamp)
   fid5(1,ipo)=fopen(tit5{1,ipo},'r+');  % Open the file
   tline = fgetl(fid5(1,ipo)); %Get the first line
   lineCounter = 1;
   Date_str_m1 = datestr(datetime(Date_str) - hours(1)); %Get the datestamp of the last line to keep
    while ischar(tline)
       if length(tline) > 19 && strcmp(tline(1:20), Date_str_m1)
         break;
       end
      % Read next line
      tline = fgetl(fid5(1,ipo));
      lineCounter = lineCounter + 1;
    end

   % Copy all the lines until the line number found in the previous section
     fid5(1,ipo) = fopen(tit5{1,ipo},'r');
     your_text = cell(lineCounter,1);
     for ii = 1:lineCounter
       your_text(ii) = {fgetl(fid5(1,ipo))}; 
     end
     fclose(fid5(1,ipo));

    % Write all the lines until the line number found in the previous section
     fid5(1,ipo) = fopen(tit5{1,ipo},'w');
        for ii = 1:lineCounter
          fprintf(fid5(1,ipo),'%s\n',your_text{ii});
        end
      clear your_text     
    end

	%%% START <<OUTPUT_ZZZ_PIXEL_YYY.dat>> %%%
    fprintf(fid5(1,ipo),'%s\t',Date_str);
    fprintf(fid5(1,ipo),'%g\t',alp_soil(ij));
    fprintf(fid5(1,ipo),'%g\t',Asur(ij));
    fprintf(fid5(1,ipo),'%g\t',Ca_S(ij));	
    fprintf(fid5(1,ipo),'%g\t',Cice(ij));	
    fprintf(fid5(1,ipo),'%g\t',Cicew(ij));	
    fprintf(fid5(1,ipo),'%g\t',CK1(ij));	
    fprintf(fid5(1,ipo),'%g\t',Csno(ij));	
    fprintf(fid5(1,ipo),'%g\t',Csnow(ij));	
    fprintf(fid5(1,ipo),'%g\t',cos_fst(ij));	
    fprintf(fid5(1,ipo),'%g\t',Ct(ij));	
    fprintf(fid5(1,ipo),'%g\t',DEB_MAP(ij));	
    fprintf(fid5(1,ipo),'%g\t',DQ_S(ij));	
    fprintf(fid5(1,ipo),'%g\t',dQ_S(ij));	
    fprintf(fid5(1,ipo),'%g\t',Ds_S(ij));	
    fprintf(fid5(1,ipo),'%g\t',DT_S(ij));	
    fprintf(fid5(1,ipo),'%g\t',dw_SNO(ij));	
    fprintf(fid5(1,ipo),'%g\t',e_sno(ij));	
    fprintf(fid5(1,ipo),'%g\t',ea_S(ij));	
    fprintf(fid5(1,ipo),'%g\t',EG(ij));
    fprintf(fid5(1,ipo),'%g\t',EICE(ij));	
    fprintf(fid5(1,ipo),'%g\t',EIn_rock(ij));
    fprintf(fid5(1,ipo),'%g\t',EIn_urb(ij));
    fprintf(fid5(1,ipo),'%g\t',er(ij));
    fprintf(fid5(1,ipo),'%g\t',ESN(ij));
    fprintf(fid5(1,ipo),'%g\t',SSN(ij));
    fprintf(fid5(1,ipo),'%g\t',ESN_In(ij));
    fprintf(fid5(1,ipo),'%g\t',SSN_In(ij));
    fprintf(fid5(1,ipo),'%g\t',EWAT(ij));
    fprintf(fid5(1,ipo),'%g\t',f(ij));
    fprintf(fid5(1,ipo),'%g\t',FROCK(ij));
    fprintf(fid5(1,ipo),'%g\t',G(ij));
    fprintf(fid5(1,ipo),'%g\t',Gfin(ij));
    fprintf(fid5(1,ipo),'%g\t',H(ij));
    fprintf(fid5(1,ipo),'%g\t',ICE(ij));
    fprintf(fid5(1,ipo),'%g\t',ICE_D(ij));
    fprintf(fid5(1,ipo),'%g\t',Imelt(ij));
    fprintf(fid5(1,ipo),'%g\t',In_rock(ij)); 
    fprintf(fid5(1,ipo),'%g\t',In_SWE(ij));
    fprintf(fid5(1,ipo),'%g\t',In_urb(ij));
    fprintf(fid5(1,ipo),'%g\t',IP_wc(ij));
    fprintf(fid5(1,ipo),'%g\t',Lk(ij));
    fprintf(fid5(1,ipo),'%g\t',Lk_rock(ij));
    fprintf(fid5(1,ipo),'%g\t',Lk_wat(ij));
    fprintf(fid5(1,ipo),'%g\t',NIce(ij));	
    fprintf(fid5(1,ipo),'%g\t',NIn_SWE(ij));
    fprintf(fid5(1,ipo),'%g\t',N_S(ij));
    fprintf(fid5(1,ipo),'%g\t',NDVI(ij));
    fprintf(fid5(1,ipo),'%g\t',OF(ij));
    fprintf(fid5(1,ipo),'%g\t',OS(ij));
    fprintf(fid5(1,ipo),'%g\t',PAR_space(ij));
    fprintf(fid5(1,ipo),'%g\t',PARB_S(ij));
    fprintf(fid5(1,ipo),'%g\t',PARD_S(ij));
    fprintf(fid5(1,ipo),'%g\t',Pre_S(ij));
    fprintf(fid5(1,ipo),'%g\t',Pr_liq(ij));
    fprintf(fid5(1,ipo),'%g\t',Pr_S(ij));
    fprintf(fid5(1,ipo),'%g\t',Pr_sno(ij));
    fprintf(fid5(1,ipo),'%g\t',QE(ij));
    fprintf(fid5(1,ipo),'%g\t',Qfm(ij));
    fprintf(fid5(1,ipo),'%g\t',QpointC(ipo));
    fprintf(fid5(1,ipo),'%g\t',QpointH(ipo));
    fprintf(fid5(1,ipo),'%g\t',Qv(ij));
    fprintf(fid5(1,ipo),'%g\t',Q_channel(ij));
    fprintf(fid5(1,ipo),'%g\t',q_runon(ij));
    fprintf(fid5(1,ipo),'%g\t',ra(ij));
    fprintf(fid5(1,ipo),'%g\t',Rd(ij));
    fprintf(fid5(1,ipo),'%g\t',Rh(ij));
    fprintf(fid5(1,ipo),'%g\t',Rn(ij));
    fprintf(fid5(1,ipo),'%g\t',ros(ij));
    fprintf(fid5(1,ipo),'%g\t',Rsw_space(ij));
    fprintf(fid5(1,ipo),'%g\t',r_soil(ij));
    fprintf(fid5(1,ipo),'%g\t',SAB1_S(ij));
    fprintf(fid5(1,ipo),'%g\t',SAB2_S(ij));
    fprintf(fid5(1,ipo),'%g\t',SAD1_S(ij));
    fprintf(fid5(1,ipo),'%g\t',SAD2_S(ij));
    fprintf(fid5(1,ipo),'%g\t',SE_rock(ij));
    fprintf(fid5(1,ipo),'%g\t',SE_urb(ij));
    fprintf(fid5(1,ipo),'%g\t',Slo_top(ij));
    fprintf(fid5(1,ipo),'%g\t',Smelt(ij));
    fprintf(fid5(1,ipo),'%g\t',SND(ij));
    fprintf(fid5(1,ipo),'%g\t',snow_albedo(ij));
    fprintf(fid5(1,ipo),'%g\t',SP_wc(ij));
    fprintf(fid5(1,ipo),'%g\t',SSN(ij));
    fprintf(fid5(1,ipo),'%g\t',surface_albedo(ij));
    fprintf(fid5(1,ipo),'%g\t',SWE(ij));
    fprintf(fid5(1,ipo),'%g\t',SWE_avalanched(ij));
    fprintf(fid5(1,ipo),'%g\t',Ta_S(ij));
    fprintf(fid5(1,ipo),'%g\t',Tdamp(ij));
    fprintf(fid5(1,ipo),'%g\t',Tdew_S(ij));
    fprintf(fid5(1,ipo),'%g\t',Tice(ij));
    fprintf(fid5(1,ipo),'%g\t',Ts(ij));
    fprintf(fid5(1,ipo),'%g\t',U_S(ij));
    fprintf(fid5(1,ipo),'%g\t',UpointC(ipo));
    fprintf(fid5(1,ipo),'%g\t',UpointH(ipo));
    fprintf(fid5(1,ipo),'%g\t',U_SWE(ij));
    fprintf(fid5(1,ipo),'%g\t',WAT(ij)); 
    fprintf(fid5(1,ipo),'%g\t',WIS(ij));	
    fprintf(fid5(1,ipo),'%g\t',WR_IP(ij));
    fprintf(fid5(1,ipo),'%g\t',WR_SP(ij));
    fprintf(fid5(1,ipo),'%g\t',Ws_S(ij));
    fprintf(fid5(1,ipo),'%g\t\n',ZWT(ij));	
	%%% END <<OUTPUT_ZZZ_PIXEL_YYY.dat>> %%%
    %%%%%
    if t==N_time_step
        fclose(fid5(1,ipo));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%

%%%%% Vars for PIXEL SOIL

vars_soil = {'Date','O','Qi_in','Qi_out','Tdp','V'};

    if t==2
        tit6{ipo}=strcat(outlocation,'OUTPUT_',SITE,'_SOILPIXEL_',strrep(POI_names(ipo)," ",""),'.dat');
        fid6(ipo)=fopen(tit6{ipo},'a');
        % Add labels to column list
        fprintf(fid6(ipo),'%s\t',vars_soil{1});

           
            for kjk = 1:ms_max
                for ii = 2:length(vars_soil)
                    fprintf(fid6(ipo),'%s\t',[vars_soil{ii} '_lr' num2str(kjk)]);
                   
                end
            end 
          fprintf(fid6(ipo),'%s\t\n','0*CK1');

    end

   if t==t1_reinit

   tit6{ipo}=strcat(outlocation,'OUTPUT_',SITE,'_SOILPIXEL_',strrep(POI_names(ipo)," ",""),'.dat');

   % Check the line number after which delete everything (based on datestamp)
   fid6(ipo)=fopen(tit6{ipo},'r+');  % Open the file
   tline = fgetl(fid6(ipo)); %Get the first line
   lineCounter = 1;
   Date_str_m1 = datestr(datetime(Date_str) - hours(1)); %Get the datestamp of the last line to keep
    while ischar(tline)
       if length(tline) > 19 && strcmp(tline(1:20), Date_str_m1)
         break;
       end
      % Read next line
      tline = fgetl(fid6(ipo));
      lineCounter = lineCounter + 1;
    end

   % Copy all the lines until the line number found in the previous section
     fid6(ipo) = fopen(tit6{ipo},'r');
     your_text = cell(lineCounter,1);
     for ii = 1:lineCounter
       your_text(ii) = {fgetl(fid6(ipo))}; 
     end
     fclose(fid6(ipo));

    % Write all the lines until the line number found in the previous section
     fid6(ipo) = fopen(tit6{ipo},'w');
        for ii = 1:lineCounter
          fprintf(fid6(ipo),'%s\n',your_text{ii});
        end
      clear your_text  
    end

    fprintf(fid6(ipo),'%s\t',Date_str);
    for kjk = 1:ms_max
	    %%% START <<OUTPUT_ZZZ_PIXEL_SOIL_YYY.dat>> %%%
        fprintf(fid6(ipo),'%g\t',O(ij,kjk));
        fprintf(fid6(ipo),'%g\t',Qi_in(ij,kjk));
        fprintf(fid6(ipo),'%g\t',Qi_out(ij,kjk));
        fprintf(fid6(ipo),'%g\t',Tdp(ij,kjk));
        fprintf(fid6(ipo),'%g\t',V(ij,kjk));
    end
    fprintf(fid6(ipo),'%g\t\n',0*CK1(ij));
	%%% END <<OUTPUT_ZZZ_PIXEL_SOIL_YYY.dat>> %%%
    if t==N_time_step
        fclose(fid6(ipo));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% DELETED: OUTPUT FOR "OUTPUT_XX_PIXEL_YY_PFT_1.DAT"...
    
end

end
