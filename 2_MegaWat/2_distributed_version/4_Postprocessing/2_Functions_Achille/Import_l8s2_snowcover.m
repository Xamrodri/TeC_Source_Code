%% Identify simulation dates for which we have Landsat/Sentinel snow cover maps

fn_maps = dir([path_L8S2 '/TIF/*RGB*.tif']);

% Extract the dates for which landsat/sentinel-2 snow cover maps are available
for ii = 1:size(fn_maps,1)
    fn_scene = fn_maps(ii).name;
    date_l8s2_char = char(extractBetween(fn_scene,[glacier '_'],'.tif'));
    sat_l8s2_char = extractBefore(fn_scene,'_RGB');
    sat_l8s2{ii} = sat_l8s2_char;
    date_l8s2(ii) = datetime(str2num(date_l8s2_char(1:4)),str2num(date_l8s2_char(6:7)),str2num(date_l8s2_char(9:10)));
end

% Extract the indices for which landsat/sentinel-2 snow maps we have T&C outputs
for ii = 1:numel(date_h)
    find_id_L8S2 = find(datetime(year(date_h(ii)),month(date_h(ii)),day(date_h(ii))) == date_l8s2',1);

    if ~isempty(find_id_L8S2)
        id_l8s2(ii,1) = find_id_L8S2;
        id_l8s2(ii,2) = ii;
    else
        id_l8s2(ii,1) = 0;
        id_l8s2(ii,2) = ii;
    end 
end
id_l8s2(id_l8s2(:,1) == 0,:) = [];

Date_l8s2 = date_h(id_l8s2(:,2));
Date_l8s2_d = date_l8s2(id_l8s2(:,1))';

% Load snow cover maps and RGB maps of the Landsat/Sentinel scenes
fn_snowmaps = dir([path_L8S2 '/TIF/*SnowCoverHybrid*.tif']);

% Load geospatial information from the first scene
[SCENE1, SCENE1_info] = readgeoraster([path_L8S2 '\TIF\' fn_snowmaps(1).name]);
SCENE1_INFO = georasterinfo([path_L8S2 '\TIF\' fn_snowmaps(1).name]);
[x_l8s2,y_l8s2] = worldGrid(SCENE1_info);
hls_GRIDobj = GRIDobj([path_L8S2 '\TIF\' fn_maps(1).name]);
DEM_l8s2 = project(DEM,hls_GRIDobj);
DEM_l8s2_orig = GRIDobj([extractBefore(path_L8S2,'_filtered') '\DEM_LS.tif']);

demInfo_l8s2 = geotiffinfo([extractBefore(path_L8S2,'_filtered') '\DEM_LS.tif']);
[demLons_l8s2,demLats_l8s2] = pixcenters(demInfo_l8s2.RefMatrix,size(DEM_l8s2_orig.Z),'makegrid');  

[l8s2_r, l8s2_g, l8s2_b, l8s2_sno] = deal(NaN(SCENE1_INFO.RasterSize(1),SCENE1_INFO.RasterSize(2), size(id_l8s2,1)));

for ii = 1:size(id_l8s2,1)
   l8s2_rgb_load=  imread([path_L8S2 '\TIF\' fn_maps(id_l8s2(ii,1)).name]);
   l8s2_r(:,:,ii) = l8s2_rgb_load(:,:,1);
   l8s2_g(:,:,ii) = l8s2_rgb_load(:,:,2);
   l8s2_b(:,:,ii) = l8s2_rgb_load(:,:,3);

   l8s2_sno_load=  imread([path_L8S2 '\TIF\' fn_snowmaps(id_l8s2(ii,1)).name]);
   l8s2_sno(:,:,ii) = l8s2_sno_load;   
end 

MASK_hls = nansum(l8s2_g,3)>0;