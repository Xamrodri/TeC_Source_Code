% NDSI_threshold = 0.40

opts = delimitedTextImportOptions("NumVariables", 5);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["dateTime", "se", "fCld", "fSnow", "fData"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, "dateTime", "InputFormat", "dd-MMM-yyyy");

% NDSI_threshold = 0.40
MOD_sca_40 = readtable([root '\Remote_sensing\MODIS\Processed_data\' glacier '\' glacier '_SE_MOD10A1_40.csv'], opts);
modis_fsc_40 = table2timetable(MOD_sca_40);

% NDSI_threshold = 0.25
MOD_sca_25 = readtable([root '\Remote_sensing\MODIS\Processed_data\' glacier '\' glacier '_SE_MOD10A1_25.csv'], opts);
modis_fsc_25 = table2timetable(MOD_sca_25);

% NDSI_threshold = 0.45
MOD_sca_45 = readtable([root '\Remote_sensing\MODIS\Processed_data\' glacier '\' glacier '_SE_MOD10A1_45.csv'], opts);
modis_fsc_45 = table2timetable(MOD_sca_45);

th_cloud = 0.1; % Discard scenes that have more than 10% cloud cover
th_data =0; % Discard scenes that have less  than th_data % data

cloud_filt = modis_fsc_40.fCld > th_cloud;
data_filt = modis_fsc_40.fData < th_data;
fsc_filt = modis_fsc_40.fSnow < 0.05; % For SLA, discard scenes with less than 5% snow cover

fsc_modis_filt_25 = modis_fsc_25.fSnow;
fsc_modis_filt_25(cloud_filt | data_filt) = NaN;
se_modis_filt_25 = modis_fsc_25.se;
se_modis_filt_25(cloud_filt | data_filt | fsc_filt) = NaN;

fsc_modis_filt_40 = modis_fsc_40.fSnow;
fsc_modis_filt_40(cloud_filt | data_filt) = NaN;% 
se_modis_filt_40 = modis_fsc_40.se;
se_modis_filt_40(cloud_filt | data_filt | fsc_filt) = NaN;

fsc_modis_filt_45 = modis_fsc_45.fSnow;
fsc_modis_filt_45(cloud_filt | data_filt) = NaN;
se_modis_filt_45 = modis_fsc_45.se;
se_modis_filt_45(cloud_filt | data_filt | fsc_filt) = NaN;

% Moving means
mv_window_modis = 10;

Modis_fsc_mv_25 = movmean(fsc_modis_filt_25,mv_window_modis,'omitnan');
Modis_fsc_mv_std_25 = movstd(fsc_modis_filt_25,mv_window_modis,'omitnan');
Modis_se_mv_25 = movmean(se_modis_filt_25,mv_window_modis,'omitnan');
Modis_se_mv_std_25 = movstd(se_modis_filt_25,mv_window_modis,'omitnan');

Modis_fsc_mv = movmean(fsc_modis_filt_40,mv_window_modis,'omitnan');
Modis_fsc_mv_std = movstd(fsc_modis_filt_40,mv_window_modis,'omitnan');
Modis_se_mv = movmean(se_modis_filt_40,mv_window_modis,'omitnan');
Modis_se_mv_std = movstd(se_modis_filt_40,mv_window_modis,'omitnan');

Modis_fsc_mv_45 = movmean(fsc_modis_filt_45,mv_window_modis,'omitnan');
Modis_fsc_mv_std_45 = movstd(fsc_modis_filt_45,mv_window_modis,'omitnan');
Modis_se_mv_45 = movmean(se_modis_filt_45,mv_window_modis,'omitnan');
Modis_se_mv_std_45 = movstd(se_modis_filt_45,mv_window_modis,'omitnan');
