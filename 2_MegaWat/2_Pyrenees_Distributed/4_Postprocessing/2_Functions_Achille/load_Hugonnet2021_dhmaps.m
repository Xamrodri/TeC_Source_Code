path_GMB = [root_B '\Study_sites\' glacier '\Spatial_data\GMB\Hugonnet_2021'];
path_gmb = [root '\Remote_sensing\Regional_global\GMB\Hugonnet_2021\Time-series'];

GMB_dem = GRIDobj([path_GMB '\' glacier '_AW3D_DEM_UTM_100m.tif']);
GMB_dem_all = GMB_dem;
demHugInfo = geotiffinfo([path_GMB '\' glacier '_AW3D_DEM_UTM_100m.tif']);
[demLonsHug,demLatsHug] = pixcenters(demHugInfo.RefMatrix,size(GMB_dem.Z),'makegrid');

GMB_hug_2000_2005 = GRIDobj([path_GMB '\Filled\' glacier '_2000_2005_dhdt_filled.tif']);
GMB_hug_2005_2010 = GRIDobj([path_GMB '\Filled\' glacier '_2005_2010_dhdt_filled.tif']);
GMB_hug_2010_2015 = GRIDobj([path_GMB '\Filled\' glacier '_2010_2015_dhdt_filled.tif']);
GMB_hug_2015_2020 = GRIDobj([path_GMB '\Filled\' glacier '_2015_2020_dhdt_filled.tif']);
GMB_hug_2000_2020 = GRIDobj([path_GMB '\Filled\' glacier '_2000_2020_dhdt_filled.tif']);

GMB_hug_all = cat(3,GMB_hug_2000_2005.Z,GMB_hug_2005_2010.Z, GMB_hug_2010_2015.Z, GMB_hug_2015_2020.Z, GMB_hug_2000_2020.Z);
clear GMB_hug_2000_2005 GMB_hug_2005_2010 GMB_hug_2010_2015 GMB_hug_2015_2020 GMB_hug_2000_2020

GLAID_hug = GRIDobj([path_GMB '\' glacier '_RGI_ID.tif']); 

GMB_hug.Z(GLAID_hug.Z ~= gla_id) = NaN; % Focus on a specific glacier, be defined further up !
GMB_dem.Z(GLAID_hug.Z ~= gla_id) = NaN;

dEL=50; % width of elevation bins
ELs_hug = nanmin(GMB_dem.Z,[],'all'):dEL:nanmax(GMB_dem.Z,[],'all');
    
for iel = 1:numel(ELs_hug)
     cur=(GMB_dem.Z<(ELs_hug(iel)+dEL/2))&(GMB_dem.Z>=(ELs_hug(iel)-dEL/2)); %current section of DEM
     pDH_hug(iel)=nanmean(GMB_hug.Z(cur)); % mean elevation change within that elevation band
end 

% Load glacier-wide mass balance and their uncertainty for each glacier

GLA_ids = unique(GLA_ID(GLA_ID>0));

% Determine the RGI zone number 
num_digits = floor(log10(GLA_ids(1))) + 1;
RGI_zone = floor(GLA_ids(1) / 10^(num_digits - 2));

opts = delimitedTextImportOptions("NumVariables", 16);
opts.DataLines = [2, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["rgiid", "period", "area", "dhdt", "err_dhdt", "dvoldt", "err_dvoldt", "dmdt", "err_dmdt", "dmdtda", "err_dmdtda", "perc_area_meas", "perc_area_res", "valid_obs", "valid_obs_py", "reg"];
opts.VariableTypes = ["categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts = setvaropts(opts, ["rgiid", "period"], "EmptyFieldRule", "auto");
Hug_gmb_table = readtable([root '\Remote_sensing\Regional_global\GMB\Hugonnet_2021\Time-series\dh_' num2str(RGI_zone) '_rgi60_pergla_rates.csv'], opts);
clear opts

Hug_rgis = string(Hug_gmb_table.rgiid);
Hug_period = string(Hug_gmb_table.period);

Hugonnet_periods_xlx = {'2000-01-01_2005-01-01','2005-01-01_2010-01-01',...
    '2010-01-01_2015-01-01','2015-01-01_2020-01-01','2000-01-01_2020-01-01'};

Hug_gmb =  [];

for ii = 1:numel(GLA_ids)

    GLA_id = num2str(GLA_ids(ii));
    ind_id =  contains(Hug_rgis,GLA_id(3:end));

    Hug_gmb_table_gla = Hug_gmb_table(ind_id,:);


    Hug_period_i = string(Hug_gmb_table.period(ind_id));
    ind_period = zeros(numel(Hug_period_i),1);

    for ss = 1:numel(Hugonnet_periods_xlx)
        ind_period = ind_period + strcmp(Hug_period_i, Hugonnet_periods_xlx{ss});
    end 
    Hug_gmb_table_gla_period = Hug_gmb_table_gla(ind_period==1,:);
    Hug_gmb.(['RGI_' num2str(GLA_id)]) = Hug_gmb_table_gla_period;
end 
