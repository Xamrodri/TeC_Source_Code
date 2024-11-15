function [SPAVG] = load_catchment_average(dir_tcout,glacier,labelled_output)
%FIRST_FUNCTION Summary of this function goes here
%   Detailed explanation goes here

opts = delimitedTextImportOptions("NumVariables", 104); % To change if the OUTPUT manager used was different
opts.DataLines = [2, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["Date", "alp_soil_tg", "Ca_tg", "Cicew_tg", "Cice_tg", "CK1_tg", "Csnow_tg", "Csno_tg", "DQ_S_tg", "dQ_S_tg", "Dr_H_tg", "Dr_L_tg", "Ds_tg", "DT_S_tg", "dw_SNO_tg", "ea_tg", "EG_tg", "EICE_tg", "EIn_H_tg", "EIn_L_tg", "EIn_rock_tg", "EIn_tg", "EIn_urb_tg", "er_tg", "ESN_tg", "SSN_tg", "EWAT_tg", "Fract_sat_tg", "FROCK_tg", "f_tg", "Gfin_tg", "G_tg", "H_tg", "ICE_D_tg", "ICE_tg", "Imelt_tg", "Inveg_tg", "In_rock_tg", "In_SWE_tg", "In_tg", "In_urb_tg", "IP_wc_tg", "Lk_rock_tg", "Lk_tg", "Lk_wat_tg", "NIce_tg", "NIn_SWE_tg", "N_tg", "NDVI_tg", "OF_tg", "OS_tg", "O_tg", "PAR_tg", "Pre_tg", "Pr_liq_tg", "Pr_sno_tg", "Pr_tg", "QE_tg", "Qfm_tg", "Qlat_in_tg", "Qlat_out_tg", "Qsub_exit", "Qv_tg", "Q_channel_tg", "Q_exit", "q_runon_tg", "ra_tg", "Rd_tg", "Rh_tg", "Rn_tg", "ros_tg", "Rsw_tg", "r_soil_tg", "SE_rock_tg", "SE_urb_tg", "SLE", "SLnoise", "Smelt_tg", "SND_tg", "snow_albedo_tg", "SP_wc_tg", "SWE_avalanched_tg", "SWE_tg", "t", "Ta_tg", "Tdamp_tg", "Tdew_tg", "Tdp_tg", "Tice_tg", "TsVEG_tg", "Ts_tg", "T_H_tg", "T_L_tg", "T_tg", "U_SWE_tg", "Vice_tg", "V_tg", "WAT_tg", "WIS_tg", "WR_IP_tg", "WR_SP_tg", "Ws_tg", "ZWT_tg", "NA"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "NA", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Date", "InputFormat", "yyyy-MM-dd HH:mm:ss");

if labelled_output == 0 % Old version of OUTPUT MANAGER, after post-processing with Pascal R script
    SPAVG = readtable([dir_tcout '\POSTPROCESSED\TABLES\SPAVG.txt'], opts);
else   
    opts = setvaropts(opts, "Date", "InputFormat", "dd-MMM-yyyy HH:mm:ss");
    SPAVG = readtable([dir_tcout '\OUTPUT_' glacier '_AVG.dat'], opts);
    Date = SPAVG.Date;
    Date_d = unique(datetime(year(Date),month(Date),day(Date)));
end 
SPAVG.NA = [];

end

