textprogressbar('Loading T&C point-scale outputs:');

if distributed == 1

    cellfn = dir([dir_tcout '\OUTPUT_' glacier '_PIXEL_*']); % T&C tracked pixels
     
    % Get the poi names of the T&C output tables
    name_poi = (extractfield(cellfn,'name'))';

    for i = 1:length(name_poi)
        name_poi_i(i,1) = extractBetween(name_poi{i},['OUTPUT_' glacier '_PIXEL_'],'.dat');
    end
    name_poi_i = string(name_poi_i);
else 
    cellfn = dir([dir_tcout '\*results.txt']);
    cellfn_param = dir([dir_tcout '\*param.txt']);
    name_poi = (extractfield(cellfn,'name'))';
    for i = 1:length(name_poi)
        name_poi_i(i,1) = extractBetween(name_poi{i},'','_results');
    end
end 

% Remove clusters POI if there is no spatialization (i.e. only validation
% at station location
if spatialize~=1
    cellfn(contains(name_poi_i,'Cluster'),:) = [];
    name_poi_i(contains(name_poi_i,'Cluster'),:) = [];
end 

% Setttings to import tracked pixels T&C table output
if distributed == 1
    opts = delimitedTextImportOptions("NumVariables", 102);
    opts.DataLines = [2, Inf];
    opts.Delimiter = "\t";
    opts.VariableNames = ["Date", "Asur", "alp_soil", "Ca_S", "Cice", "Cicew", "CK1", "Csno", "Csnow", "cos_fst", "Ct", "DEB", "DQ_S", "dQ_S", "Ds_S", "DT_S", "dw_SNO", "e_sno", "ea_S", "EG", "EICE", "EIn_rock", "EIn_urb", "er", "ESN", "SSN", "ESN_In", "SSN_In", "EWAT", "f", "FROCK", "G", "Gfin", "H", "ICE", "ICE_D", "Imelt", "In_rock", "In_SWE", "In_urb", "IP_wc", "Lk", "Lk_rock", "Lk_wat", "NIce", "NIn_SWE", "N_S", "NDVI", "OF", "OS", "PAR_space", "PARB_S", "PARD_S", "Pre_S", "Pr_liq", "Pr_S", "Pr_sno", "QE", "Qfm", "QpointC", "QpointH", "Qv", "Q_channel", "q_runon", "ra", "Rd", "Rh", "Rn", "ros", "Rsw_space", "r_soil", "SAB1_S", "SAB2_S", "SAD1_S", "SAD2_S", "SE_rock", "SE_urb", "Slo_top_S", "Smelt", "SND", "snow_albedo", "SP_wc", "SSN1", "surface_albedo", "SWE", "SWE_avalanched", "Ta_S", "Tdamp", "Tdew_S", "Tice", "Ts", "U_S", "UpointC", "UpointH", "U_SWE", "WAT", "WIS", "WR_IP", "WR_SP", "Ws_S", "ZWT", "VarName102"];
    opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    opts = setvaropts(opts, "VarName102", "WhitespaceRule", "preserve");
    opts = setvaropts(opts, "VarName102", "EmptyFieldRule", "auto");
    opts = setvaropts(opts, "Date", "InputFormat", "dd-MMM-yyyy HH:mm:ss"); 
end

l = 0;

table_poi_i = cell(length(name_poi_i),3);
table_poi_dm = cell(length(name_poi_i),3);
table_poi_ds = cell(length(name_poi_i),3);

for i = 1:length(name_poi_i)
   textprogressbar((i/length(name_poi_i))*100)
   fn_poi = cellfn(i).name;

   if distributed == 1
       table_poi = readtable([dir_tcout '\' fn_poi],opts);
   else 
       table_poi = readtable([dir_tcout '\' fn_poi]);
   end 
       ind_dtm = find(strcmp(name_poi_i(i),POI_table.Name),1);   
       Date = table_poi.Date;

% To add back: SSN, snow_albedo, Imelt
    if distributed == 1; table_poi_i{i,1} = table2timetable(table_poi(:, {'Date','Ta_S','Ts', 'Pr_S','Pr_sno','SND','Rsw_space','PAR_space','Ws_S','U_S','Pre_S','N_S','SWE','ICE','snow_albedo','QpointC','Q_channel','Smelt','Imelt','ESN'})) ;
          table_poi_i{i,1}.Properties.VariableNames = {'Ta_poi','Ts_poi', 'Pr_poi','Pr_sno_poi','SND_poi','Rsw_poi','PAR_poi','WS_poi','RH_poi','Pre_poi','LWIN_poi','SWE','ICE','Albedo','QpointC','Q_channel','Smelt','Imelt','ESN'};
    else 
          table_poi_i{i,1} = table2timetable(table_poi(:, {'Date','Ta','Pre', 'Pr','Pr_sno','SND','Ws','U','N','SWE','ICE','Albedo','Smelt','ESN','EICE','SSN','Imelt','ET'}));
          table_poi_i{i,1}.Properties.VariableNames = {'Ta_poi','Pre_poi', 'Pr_poi','Pr_sno_poi','SND_poi','WS_poi','RH_poi','LWIN_poi','SWE','ICE','Albedo','Smelt_poi','ESN_poi','EICE','SSN_poi','Imelt','ET'};
          table_poi_i{i,1}.Rsw_poi = table_poi.Rsw;   
    end
    table_poi_i{i,2} = name_poi_i(i);
    table_poi_i{i,3} = POI_table.ELEV(ind_dtm);
    table_poi_i{i,4} = Asur(POI_table.ROW(ind_dtm),POI_table.COL(ind_dtm));

    table_poi_dm{i,1} = retime(table_poi_i{i,1},'daily',@nanmean);
    table_poi_dm{i,2} = table_poi_i{i,2};
    table_poi_dm{i,3} = table_poi_i{i,3};   

    table_poi_ds{i,1} = retime(table_poi_i{i,1},'daily',@nansum);
    table_poi_ds{i,2} = table_poi_i{i,2};
    table_poi_ds{i,3} = table_poi_i{i,3};
end

% table_poi_i = sortrows(table_poi_i,3);

table_poi_d = table_poi_dm;

for i = 1:length(name_poi_i)
    table_poi_d{i,1}.Pr_poi = table_poi_ds{i,1}.Pr_poi;
    table_poi_d{i,1}.Pr_sno_poi = table_poi_ds{i,1}.Pr_sno_poi;
end 

% Check which POI from the table were not found in the T&C multi-points output

name_poi_outputed = extractBefore(name_poi,'_results');
l = 1;
poi_not_found = [];

for ii = 1:size(POI_table,1)
   poi_found_id = find(strcmp(POI_table.Name{ii}, name_poi_outputed),1);

   if isempty(poi_found_id)
       poi_not_found{l,1} = ii;
       poi_not_found{l,2} = POI_table.Name{ii};
       l = l +1;
   end
end 

if ~isempty(poi_not_found)
    disp([num2str(size(poi_not_found,1)) ' POIs were not found: ' strjoin(poi_not_found(:,2),' '),' with loc: '  num2str([poi_not_found{:,1}]) ])
end 

textprogressbar('Finished')
