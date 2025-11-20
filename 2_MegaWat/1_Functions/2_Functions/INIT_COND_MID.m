%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Saving condition of the model for potential restarting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = INIT_COND_MID(t_store, Directories, ...
     alp_soil,      b_soil,      Bam,             Bem,            BLit, ...
     Ccrown_t,      Cice,        Cicew,           CK1,            Csno,...
     Csnow,         dQ_S,        DQ_S,            dQVEG,          DT_S,...
     dw_SNO,        e_sno,       EG,              EICE,           EIn_rock,...
     EIn_urb,       EK,          ELitter,         er,             ESN_In,...
     SSN_In,        ESN,         SSN,             EWAT,           FROCK,...
     f,             Gfin,        G,               H,              HV,...
     ICE_D,         ICE,         Imelt,           In_H,           In_Litter,...
     In_L,          In_rock,     In_SWE,          In_urb,         IP_wc,...
     Lk_rock,       Lk_wat,      Lk,              Lpho,           NavlI,...
     NDVI,          NIce,        NIn_SWE,         OF,             Oice,...
     OS,            O,           POT,             Pr_liq,         Pr_sno,...
     Q_channel,     Q_exit,      q_runon,         QE,             QEV,...
     Qfm,           Qi_in,       Qi_out_Rout,     Qi_out,...
     Qsub_exit,     Qv,          r_litter,        r_soil,         ra,...
     Rd,            Rh,          Rn,              ros,            SE_rock,...
     SE_urb,        Slo_head,    Smelt,           SND,            snow_albedo, ...
     soil_albedo,   SP_wc,       surface_albedo,  SWE,            SWE_avalanched, ...
     t_sls,         tau_sno,     Tdamp,           Tdeb,           Tdp, ...
     Tdp_snow,      Tice,        Tstm0,           Ts,             Ts_under, ...
     TsVEG,         U_SWE,       Vice,            V,              WAT, ...
     WIS,           WR_IP,       WR_SP,           Ws_under,       ZWT, ...
     AgeDL_H,       AgeDL_L,     AgeL_H,          AgeL_L,         AgePl_H,...
     AgePl_L,       An_H,        An_L,            ANPP_H,         ANPP_L,...
     B_H,           B_L,         BA_H,            BA_L,           Bfac_dayH,...
     Bfac_dayL,     Bfac_weekH,  Bfac_weekL,      Ci_shdH,        Ci_shdL,...
     Ci_sunH,       Ci_sunL,     dflo_H,          dflo_L,         Dr_H,...
     Dr_L,          e_rel_H,     e_rel_L,         e_relN_H,       e_relN_L,...
     EIn_H,         EIn_L,       fapar_H,         fapar_L,        FNC_H,...
     FNC_L,         gsr_H,       gsr_L,           hc_H,           hc_L,...
     ISOIL_H,       ISOIL_L,     Jsx_H,           Jsx_L,          Jxl_H, ...
     Jxl_L,         Kleaf_H,     Kleaf_L,         Kreserve_H,     Kreserve_L, ...
     Kuptake_H,     Kuptake_L,   Kx_H,            Kx_L,           LAI_H, ...
     LAI_L,         LAIdead_H,   LAIdead_L,       ManIH,          ManIL, ...
     NBLeaf_H,      NBLeaf_L,    NBLI_H,          NBLI_L,         NPP_H, ...
     NPP_L,         NPPI_H,      NPPI_L,          Nreserve_H,     Nreserve_L, ...
     NuLit_H,       NuLit_L,     NupI_H,          NupI_L,         Nuptake_H, ...
     Nuptake_L,     OH,          OL,              PARI_H,         PARI_L, ...
     PHE_S_H,       PHE_S_L,     Preserve_H,      Preserve_L,     Psi_l_H, ...
     Psi_l_L,       Psi_s_H,     Psi_s_L,         Psi_x_H,        Psi_x_L, ...
     Puptake_H,     Puptake_L,   RA_H,            RA_L,           rap_H, ...
     rap_L,         RB_H,        rb_H,            RB_L,           rb_L, ...
     Rdark_H,       Rdark_L,     Rexmy_H,         Rexmy_L,        Rg_H, ...
     Rg_L,          rKc_H,       rKc_L,           Rmc_H,          Rmc_L, ...
     Rmr_H,         Rmr_L,       Rms_H,           Rms_L,          rNc_H, ...
     rNc_L,         rPc_H,       rPc_L,           Rrootl_H,       Rrootl_L, ...
     rs_shdH,       rs_shdL,     rs_sunH,         rs_sunL,        SAI_H, ...
     SAI_L,         Sfr_H,       Sfr_L,           SIF_H,          SIF_L, ...
     Slf_H,         Slf_L,       Sll_H,           Sll_L,          Sr_H, ...
     Sr_L,          SupK_H,      SupK_L,          SupN_H,         SupN_L, ...
     SupP_H,        SupP_L,      Swm_H,           Swm_L,          T_H, ...
     T_L,           TBio_H,      TBio_L,          Tden_H,         Tden_L, ...
     Tdp_H,         Tdp_L,       TdpI_H,          TdpI_L,         TexC_H, ...
     TexC_L,        TexK_H,      TexK_L,          TexN_H,         TexN_L, ...
     TexP_H,        TexP_L,      TNIT_H,          TNIT_L,         TPHO_H, ...
     TPHO_L,        TPOT_H,      TPOT_L,          Vl_H,           Vl_L, ...
     Vx_H,          Vx_L,        An_H_t,          An_L_t,         O_t, ...
     PAR_t,         Pr_sno_t,    Psi_l_H_t,       Psi_l_L_t,      Psi_x_H_t, ...
     Psi_x_L_t,     Rdark_H_t,   Rdark_L_t,       Ta_t,           Tdp_H_t, ...
     Tdp_L_t,       Tdp_t,       V_t,             Ared,           ICEym1, ...
     SNOWALB,       SWEym1 ...
        )

        %% Assigning variables 
        alp_soiltm1    = alp_soil;
        b_soiltm1      = b_soil;
        Bamtm1         = Bam;             
        Bemtm1         = Bem;             
        BLittm1        = BLit;            
        Ccrown_t_tm1   = Ccrown_t;	     
        Cicetm1        = Cice;            
        Cicewtm1       = Cicew;           
        CK1tm1         = CK1;            
        Csnotm1        = Csno;            
        Csnowtm1       = Csnow;           
        dQ_Stm1        = dQ_S;           
        DQ_Stm1        = DQ_S;            
        dQVEGtm1       = dQVEG;           
        DT_Stm1        = DT_S;           
        dw_SNOtm1      = dw_SNO;         
        e_snotm1       = e_sno;          
        EGtm1          = EG;        
        EICEtm1        = EICE;            
        EIn_rocktm1    = EIn_rock;	     
        EIn_urbtm1     = EIn_urb;     
        EKtm1          = EK;     
        ELittertm1     = ELitter;   
        ertm1          = er;     
        ESN_Intm1      = ESN_In;       
        SSN_Intm1      = SSN_In;      
        ESNtm1         = ESN;        
        SSNtm1         = SSN;               
        EWATtm1        = EWAT;            
        FROCKtm1       = FROCK;           
        ftm1           = f;        
        Gfintm1        = Gfin;            
        Gtm1           = G;          
        Htm1           = H;               
        HVtm1          = HV;              
        ICE_Dtm1       = ICE_D;           
        ICEtm1         = ICE;             
        Imelttm1       = Imelt;           
        In_Htm1        = In_H;            
        In_Littertm1   = In_Litter;     
        In_Ltm1        = In_L;           
        In_rocktm1     = In_rock;    
        In_SWEtm1      = In_SWE;         
        In_urbtm1      = In_urb;        
        IP_wctm1       = IP_wc;          
        Lk_rocktm1     = Lk_rock;    
        Lk_wattm1      = Lk_wat;        
        Lktm1          = Lk;              
        Lphotm1        = Lpho;          
        NavlItm1       = NavlI;         
        NDVItm1        = NDVI;           
        NIcetm1        = NIce;           
        NIn_SWEtm1     = NIn_SWE;	     
        OFtm1          = OF;             
        Oicetm1        = Oice;            
        OStm1          = OS;           
        Otm1           = O;               
        POTtm1         = POT;             
        Pr_liqtm1      = Pr_liq;          
        Pr_snotm1      = Pr_sno;          
        %Q_channel      = Q_channel;	     
        %Q_exit         = Q_exit;          
        %q_runon        = q_runon;         
        QEtm1          = QE;              
        QEVtm1         = QEV;            
        Qfmtm1         = Qfm;             
        %Qi_in          = Qi_in;           
        %Qi_out         = Qi_in_Ro; % Qi_in_Ro seems not used and it only overwrite Qi_out
        %Qi_out_Rout    = Qi_out_Rout;     
        Qi_outtm1      = Qi_out;          
        %Qsub_exit      = Qsub_exit;	     
        Qvtm1          = Qv;             
        r_littertm1    = r_litter;	     
        r_soiltm1      = r_soil;          
        ratm1          = ra;              
        Rdtm1          = Rd;              
        Rhtm1          = Rh;              
        Rntm1          = Rn;              
        rostm1         = ros;             
        SE_rocktm1     = SE_rock;    
        SE_urbtm1      = SE_urb;         
        %Slo_head      = Slo_head;	     
        Smelttm1       = Smelt;           
        SNDtm1         = SND;            
        snow_albedotm1 = snow_albedo;   
        soil_albedotm1 = soil_albedo;     
        SP_wctm1       = SP_wc;           
        surface_albedotm1 = surface_albedo;  
        SWEtm1         = SWE;
        SWE_avalanchedtm1 = SWE_avalanched; 
        t_slstm1       = t_sls;
        tau_snotm1     = tau_sno;	     
        Tdamptm1       = Tdamp;           
        Tdebtm1        = Tdeb;            
        Tdptm1         = Tdp;             
        Tdp_snowtm1    = Tdp_snow;        
        Ticetm1        = Tice;           
        %Tstm0          = Tstm0;           
        Tstm1          = Ts;             
        Ts_undertm1    = Ts_under;      
        TsVEGtm1       = TsVEG;     
        U_SWEtm1       = U_SWE;           
        Vicetm1        = Vice;          
        Vtm1           = V;          
        WATtm1         = WAT;             
        WIStm1         = WIS;             
        WR_IPtm1       = WR_IP;           
        WR_SPtm1       = WR_SP;          
        Ws_undertm1    = Ws_under;	    
        ZWTtm1         = ZWT;            
        
        % Specifications for high & low vegetation
        %------------------------------------------------------------------
        AgeDL_Htm1     = AgeDL_H;
        AgeDL_Ltm1     = AgeDL_L;	
        AgeL_Htm1      = AgeL_H;    
        AgeL_Ltm1      = AgeL_L;     
        AgePl_Htm1     = AgePl_H;	
        AgePl_Ltm1     = AgePl_L;	
        An_Htm1        = An_H;      
        An_Ltm1        = An_L;       
        ANPP_Htm1      = ANPP_H;    
        ANPP_Ltm1      = ANPP_L;     
        B_Htm1         = B_H;        
        B_Ltm1         = B_L;        
        BA_Htm1        = BA_H;       
        BA_Ltm1        = BA_L;       
        Bfac_dayHtm1   = Bfac_dayH;	
        Bfac_dayLtm1   = Bfac_dayL;	
        Bfac_weekHtm1  = Bfac_weekH;	
        Bfac_weekLtm1  = Bfac_weekL;	
        Citm1_shdH     = Ci_shdH;
        Citm1_shdL     = Ci_shdL;	
        Citm1_sunH     = Ci_sunH;	
        Citm1_sunL     = Ci_sunL;	
        dflo_Htm1      = dflo_H;    
        dflo_Ltm1      = dflo_L;     
        Dr_Htm1        = Dr_H;       
        Dr_Ltm1        = Dr_L;       
        e_rel_Htm1     = e_rel_H;
        e_rel_Ltm1     = e_rel_L;	
        e_relN_Htm1    = e_relN_H;	
        e_relN_Ltm1    = e_relN_L;	
        EIn_Htm1       = EIn_H;    
        EIn_Ltm1       = EIn_L;      
        fapar_Htm1     = fapar_H;	
        fapar_Ltm1     = fapar_L;	
        FNC_Htm1       = FNC_H;   
        FNC_Ltm1       = FNC_L;      
        gsr_Htm1       = gsr_H;      
        gsr_Ltm1       = gsr_L;      
        hc_Htm1        = hc_H;      
        hc_Ltm1        = hc_L;  
        ISOIL_Htm1     = ISOIL_H;	
        ISOIL_Ltm1     = ISOIL_L;	
        Jsx_Htm1       = Jsx_H;    
        Jsx_Ltm1       = Jsx_L;      
        Jxl_Htm1       = Jxl_H;      
        Jxl_Ltm1       = Jxl_L;      
        Kleaf_Htm1     = Kleaf_H;	
        Kleaf_Ltm1     = Kleaf_L;	
        Kreserve_Htm1  = Kreserve_H;	
        Kreserve_Ltm1  = Kreserve_L;	
        Kuptake_Htm1   = Kuptake_H;	
        Kuptake_Ltm1   = Kuptake_L;	
        Kx_Htm1        = Kx_H;      
        Kx_Ltm1        = Kx_L;       
        LAI_Htm1       = LAI_H;      
        LAI_Ltm1       = LAI_L;     
        LAIdead_Htm1   = LAIdead_H;
        LAIdead_Ltm1   = LAIdead_L;	
        ManIHtm1       = ManIH;     
        ManILtm1       = ManIL;      
        NBLeaf_Htm1    = NBLeaf_H;	
        NBLeaf_Ltm1    = NBLeaf_L;	
        NBLI_Htm1      = NBLI_H;    
        NBLI_Ltm1      = NBLI_L;     
        NPP_Htm1       = NPP_H;    
        NPP_Ltm1       = NPP_L;      
        NPPI_Htm1      = NPPI_H;     
        NPPI_Ltm1      = NPPI_L;     
        Nreserve_Htm1  = Nreserve_H;
        Nreserve_Ltm1  = Nreserve_L;	
        NuLit_Htm1     = NuLit_H;
        NuLit_Ltm1     = NuLit_L;	
        NupI_Htm1      = NupI_H;     
        NupI_Ltm1      = NupI_L;     
        Nuptake_Htm1   = Nuptake_H;
        Nuptake_Ltm1   = Nuptake_L;	
        OHtm1          = OH;         
        OLtm1          = OL;         
        PARI_Htm1      = PARI_H; 
        PARI_Ltm1      = PARI_L;     
        PHE_S_Htm1     = PHE_S_H;	
        PHE_S_Ltm1     = PHE_S_L;	
        Preserve_Htm1  = Preserve_H;
        Preserve_Ltm1  = Preserve_L;	
        Psi_l_Htm1     = Psi_l_H;
        Psi_l_Ltm1     = Psi_l_L;	
        Psi_s_Htm1     = Psi_s_H;	
        Psi_s_Ltm1     = Psi_s_L;	
        Psi_x_Htm1     = Psi_x_H;	
        Psi_x_Ltm1     = Psi_x_L;	
        Puptake_Htm1   = Puptake_H;
        Puptake_Ltm1   = Puptake_L;
        RA_Htm1        = RA_H;       
        RA_Ltm1        = RA_L;       
        rap_Htm1       = rap_H;      
        rap_Ltm1       = rap_L;     
        RB_Htm1        = RB_H;       
        rb_Htm1        = rb_H;       
        RB_Ltm1        = RB_L;       
        rb_Ltm1        = rb_L;       
        Rdark_Htm1     = Rdark_H;
        Rdark_Ltm1     = Rdark_L;	
        Rexmy_Htm1     = Rexmy_H;	
        Rexmy_Ltm1     = Rexmy_L;	
        Rg_Htm1        = Rg_H;     
        Rg_Ltm1        = Rg_L;       
        rKc_Htm1       = rKc_H;      
        rKc_Ltm1       = rKc_L;     
        Rmc_Htm1       = Rmc_H;     
        Rmc_Ltm1       = Rmc_L;     
        Rmr_Htm1       = Rmr_H;    
        Rmr_Ltm1       = Rmr_L;     
        Rms_Htm1       = Rms_H;     
        Rms_Ltm1       = Rms_L;     
        rNc_Htm1       = rNc_H;     
        rNc_Ltm1       = rNc_L;    
        rPc_Htm1       = rPc_H;
        rPc_Ltm1       = rPc_L;    
        Rrootl_Htm1    = Rrootl_H;	
        Rrootl_Ltm1    = Rrootl_L;	
        rs_shdHtm1     = rs_shdH;
        rs_shdLtm1     = rs_shdL;	
        rs_sunHtm1     = rs_sunH;	
        rs_sunLtm1     = rs_sunL;	
        SAI_Htm1       = SAI_H;   
        SAI_Ltm1       = SAI_L;     
        Sfr_Htm1       = Sfr_H;     
        Sfr_Ltm1       = Sfr_L;     
        SIF_Htm1       = SIF_H;     
        SIF_Ltm1       = SIF_L;     
        Slf_Htm1       = Slf_H;     
        Slf_Ltm1       = Slf_L;     
        Sll_Htm1       = Sll_H;     
        Sll_Ltm1       = Sll_L;     
        Sr_Htm1        = Sr_H;    
        Sr_Ltm1        = Sr_L;      
        SupK_Htm1      = SupK_H;     
        SupK_Ltm1      = SupK_L;     
        SupN_Htm1      = SupN_H;    
        SupN_Ltm1      = SupN_L;    
        SupP_Htm1      = SupP_H;    
        SupP_Ltm1      = SupP_L;    
        Swm_Htm1       = Swm_H;   
        Swm_Ltm1       = Swm_L;     
        T_Htm1         = T_H;   
        T_Ltm1         = T_L;       
        TBio_Htm1      = TBio_H;  
        TBio_Ltm1      = TBio_L;    
        Tden_Htm1      = Tden_H;    
        Tden_Ltm1      = Tden_L;    
        Tdp_Htm1       = Tdp_H;    
        Tdp_Ltm1       = Tdp_L;     
        TdpI_Htm1      = TdpI_H;     
        TdpI_Ltm1      = TdpI_L;     
        TexC_Htm1      = TexC_H;   
        TexC_Ltm1      = TexC_L;    
        TexK_Htm1      = TexK_H;    
        TexK_Ltm1      = TexK_L;   
        TexN_Htm1      = TexN_H;    
        TexN_Ltm1      = TexN_L;    
        TexP_Htm1      = TexP_H;    
        TexP_Ltm1      = TexP_L;    
        TNIT_Htm1      = TNIT_H;    
        TNIT_Ltm1      = TNIT_L;    
        TPHO_Htm1      = TPHO_H;    
        TPHO_Ltm1      = TPHO_L;    
        TPOT_Htm1      = TPOT_H;    
        TPOT_Ltm1      = TPOT_L;    
        Vl_Htm1        = Vl_H;   
        Vl_Ltm1        = Vl_L;      
        Vx_Htm1        = Vx_H;      
        Vx_Ltm1        = Vx_L;      

%% Storing
save([Directories.save 'Store/Initial_Conditions_' char(num2str(month(t_store))) '_' char(num2str(year(t_store))) '.mat'], ... 
'alp_soiltm1',        'b_soiltm1',         'Bamtm1',            'Bemtm1',       'BLittm1' ,...
'Ccrown_t_tm1',       'Cicetm1',           'Cicewtm1',          'CK1tm1',       'Csnotm1',...
'Csnowtm1',           'dQ_Stm1',           'DQ_Stm1',           'dQVEGtm1',     'DT_Stm1',...
'dw_SNOtm1',          'e_snotm1',          'EGtm1',             'EICEtm1',      'EIn_rocktm1',...
'EIn_urbtm1',         'EKtm1',             'ELittertm1',        'ertm1' ,       'ESN_Intm1',...
'SSN_Intm1',          'ESNtm1',            'SSNtm1' ,           'EWATtm1',      'FROCKtm1',...
'ftm1' ,              'Gfintm1' ,          'Gtm1' ,             'Htm1',         'HVtm1',...
'ICE_Dtm1',           'ICEtm1',            'Imelttm1',          'In_Htm1' ,     'In_Littertm1',...
'In_Ltm1',            'In_rocktm1',        'In_SWEtm1',         'In_urbtm1',    'IP_wctm1',...
'Lk_rocktm1',         'Lk_wattm1',         'Lktm1',             'Lphotm1',      'NavlItm1',...
'NDVItm1' ,           'NIcetm1' ,          'NIn_SWEtm1',        'OFtm1' ,       'Oicetm1',...
'OStm1',              'Otm1' ,             'POTtm1' ,           'Pr_liqtm1',    'Pr_snotm1',...
'Q_channel' ,         'Q_exit' ,           'q_runon' ,          'QEtm1',        'QEVtm1',...
'Qfmtm1',             'Qi_in',             'Qi_out',            'Qi_out_Rout',  'Qi_outtm1',...
'Qsub_exit',          'Qvtm1' ,            'r_littertm1' ,      'r_soiltm1' ,   'ratm1',...
'Rdtm1',              'Rhtm1',             'Rntm1',             'rostm1' ,      'SE_rocktm1',...
'SE_urbtm1' ,         'Slo_head' ,         'Smelttm1',          'SNDtm1',       'snow_albedotm1',...
'soil_albedotm1',     'SP_wctm1',          'surface_albedotm1', 'SWEtm1',       'SWE_avalanchedtm1',...
't_slstm1',           'tau_snotm1' ,       'Tdamptm1' ,         'Tdebtm1' ,     'Tdptm1',...
'Tdp_snowtm1',        'Ticetm1' ,          'Tstm0' ,            'Tstm1',        'Ts_undertm1',...
'TsVEGtm1',           'U_SWEtm1' ,         'Vicetm1' ,          'Vtm1' ,        'WATtm1',...
'WIStm1' ,            'WR_IPtm1',          'WR_SPtm1' ,         'Ws_undertm1',  'ZWTtm1',...   
'AgeDL_Htm1',         'AgeDL_Ltm1' ,       'AgeL_Htm1' ,        'AgeL_Ltm1' ,   'AgePl_Htm1',...
'AgePl_Ltm1',         'An_Htm1',           'An_Ltm1' ,          'ANPP_Htm1' ,   'ANPP_Ltm1',...
'B_Htm1' ,            'B_Ltm1' ,           'BA_Htm1' ,          'BA_Ltm1',      'Bfac_dayHtm1',...
'Bfac_dayLtm1',       'Bfac_weekHtm1' ,    'Bfac_weekLtm1' ,    'Citm1_shdH' ,  'Citm1_shdL',...
'Citm1_sunH' ,        'Citm1_sunL' ,       'dflo_Htm1' ,        'dflo_Ltm1',    'Dr_Htm1',...
'Dr_Ltm1' ,           'e_rel_Htm1',        'e_rel_Ltm1',        'e_relN_Htm1',  'e_relN_Ltm1',...
'EIn_Htm1',           'EIn_Ltm1' ,         'fapar_Htm1' ,       'fapar_Ltm1' ,  'FNC_Htm1',...
'FNC_Ltm1',           'gsr_Htm1',          'gsr_Ltm1' ,         'hc_Htm1' ,     'hc_Ltm1',...
'In_Htm1' ,           'In_Ltm1',           'ISOIL_Htm1' ,       'ISOIL_Ltm1' ,  'Jsx_Htm1',...
'Jsx_Ltm1',           'Jxl_Htm1',          'Jxl_Ltm1' ,         'Kleaf_Htm1',   'Kleaf_Ltm1',...
'Kreserve_Htm1',      'Kreserve_Ltm1',     'Kuptake_Htm1' ,     'Kuptake_Ltm1', 'Kx_Htm1',...
'Kx_Ltm1' ,           'LAI_Htm1' ,         'LAI_Ltm1' ,         'LAIdead_Htm1' ,'LAIdead_Ltm1',...
'ManIHtm1' ,          'ManILtm1' ,         'NBLeaf_Htm1',       'NBLeaf_Ltm1' , 'NBLI_Htm1',...
'NBLI_Ltm1',          'NPP_Htm1' ,         'NPP_Ltm1' ,         'NPPI_Htm1' ,   'NPPI_Ltm1',...
'Nreserve_Htm1',      'Nreserve_Ltm1',     'NuLit_Htm1' ,       'NuLit_Ltm1',   'NupI_Htm1',...
'NupI_Ltm1' ,         'Nuptake_Htm1' ,     'Nuptake_Ltm1' ,     'OHtm1',        'OLtm1',...
'PARI_Htm1',          'PARI_Ltm1' ,        'PHE_S_Htm1',        'PHE_S_Ltm1' ,  'Preserve_Htm1',...
'Preserve_Ltm1',      'Psi_l_Htm1' ,       'Psi_l_Ltm1',        'Psi_s_Htm1' ,  'Psi_s_Ltm1',...
'Psi_x_Htm1',         'Psi_x_Ltm1',        'Puptake_Htm1',      'Puptake_Ltm1', 'RA_Htm1',...
'RA_Ltm1',            'rap_Htm1' ,         'rap_Ltm1',          'RB_Htm1' ,     'rb_Htm1',...
'RB_Ltm1',            'rb_Ltm1' ,          'Rdark_Htm1' ,       'Rdark_Ltm1' ,  'Rexmy_Htm1',...
'Rexmy_Ltm1',         'Rg_Htm1',           'Rg_Ltm1' ,          'rKc_Htm1' ,    'rKc_Ltm1',...
'Rmc_Htm1' ,          'Rmc_Ltm1' ,         'Rmr_Htm1'  ,        'Rmr_Ltm1',     'Rms_Htm1',...
'Rms_Ltm1' ,          'rNc_Htm1',          'rNc_Ltm1' ,         'rPc_Htm1',     'rPc_Ltm1',...
'Rrootl_Htm1',        'Rrootl_Ltm1',       'rs_shdHtm1' ,       'rs_shdLtm1' ,  'rs_sunHtm1',...
'rs_sunLtm1',         'SAI_Htm1',          'SAI_Ltm1',          'Sfr_Htm1',     'Sfr_Ltm1',...
'SIF_Htm1',           'SIF_Ltm1',          'Slf_Htm1' ,         'Slf_Ltm1' ,    'Sll_Htm1',...
'Sll_Ltm1' ,          'Sr_Htm1',           'Sr_Ltm1'  ,         'SupK_Htm1' ,   'SupK_Ltm1',...
'SupN_Htm1',          'SupN_Ltm1' ,        'SupP_Htm1' ,        'SupP_Ltm1' ,   'Swm_Htm1',...
'Swm_Ltm1',           'T_Htm1' ,           'T_Ltm1'  ,          'TBio_Htm1'  ,  'TBio_Ltm1',...
'Tden_Htm1',          'Tden_Ltm1' ,        'Tdp_Htm1' ,         'Tdp_Ltm1',     'TdpI_Htm1',...
'TdpI_Ltm1',          'TexC_Htm1',         'TexC_Ltm1' ,        'TexK_Htm1' ,   'TexK_Ltm1',...
'TexN_Htm1' ,         'TexN_Ltm1' ,        'TexP_Htm1',         'TexP_Ltm1' ,   'TNIT_Htm1',...
'TNIT_Ltm1',          'TPHO_Htm1',         'TPHO_Ltm1' ,        'TPOT_Htm1' ,   'TPOT_Ltm1',...
'Vl_Htm1',            'Vl_Ltm1',           'Vx_Htm1' ,          'Vx_Ltm1',      'An_H_t', ...
'An_L_t',             'O_t',               'PAR_t',             'Pr_sno_t',     'Psi_l_H_t', ...
'Psi_l_L_t',          'Psi_x_H_t',         'Psi_x_L_t',         'Rdark_H_t',    'Rdark_L_t', ...
'Ta_t',               'Tdp_H_t',           'Tdp_L_t',           'Tdp_t',        'V_t', ...
'Ared',               'ICEym1',            'SNOWALB',           'SWEym1'...
            );

disp(['Variables correctly saved for ' char(num2str(day(t_store))) '_' char(num2str(month(t_store))) '_' char(num2str(year(t_store))) ])

end