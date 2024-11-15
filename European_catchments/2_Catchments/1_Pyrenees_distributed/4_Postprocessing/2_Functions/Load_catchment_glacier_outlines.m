if strcmp(glacier,'Kyzylsu')
    gla_outline_path = [root '\Study_sites\' glacier '\Spatial_data\Shapefiles\glacier_Kyzylsu_catchment_RGI60_UTM42N.shp'];
elseif strcmp(glacier,'Rolwaling')
    gla_outline_path = [root_B '\Study_sites\' glacier '\Spatial_data\Shapefiles\Rolwaling_glacier_mask_200m_32645.shp'];
    %gla_outline_path = [root '\Study_sites\' glacier '\Spatial_data\Shapefiles\Gandam_RolwalingRegion_UTM45N.shp'];
elseif strcmp(glacier,'Parlung4')
    gla_outline_path = [root_B '\Study_sites\' glacier '\Spatial_data\Shapefiles\Parlung4_RGI60_glacier_outlines_UTM47N.shp'];
elseif strcmp(glacier,'Mugagangqiong')
    gla_outline_path = [root '\Study_sites\' glacier '\Spatial_data\GlacierOutlines\Mugagangqiong_glacier_outlines_RGI60_utm45.shp'];
end 

if exist('gla_outline_path','var')
    glaciers = shaperead(gla_outline_path); 

% Add a column with the glacier ID 
if isfield(glaciers,'RGIID') 
    for i = 1:size(glaciers,1)
      RGI_Id_str = glaciers(i).RGIID;
      RGI_Id(i) = str2num(RGI_Id_str([7:8 10:end]));
    end
    C = num2cell(RGI_Id'); [glaciers.('RGIID_num')] = C{:};
elseif  isfield(glaciers,'RGIId') 
    for i = 1:size(glaciers,1)
      RGI_Id_str = glaciers(i).RGIId;
      RGI_Id(i) = str2num(RGI_Id_str([7:8 10:end]));
    end
    C = num2cell(RGI_Id'); [glaciers.('RGIID_num')] = C{:};
elseif isfield(glaciers,'FID_1')
     [glaciers.('RGIID_num')] = glaciers.FID_1;
end 
end

% Load catchment shapefile (native resolution, not super useful here)

if strcmp(glacier,'Kyzylsu')
    catchShp_path = [root '\Study_sites\Kyzylsu\Spatial_data\Masks\Kyzylsu_catchment_mask_100m_32642.shp'];
elseif strcmp(glacier,'Rolwaling')
    catchShp_path = [root_B '\Study_sites\' glacier '\Spatial_data\Masks\Rolwaling_catchment_mask_200m_32645.shp'];
elseif strcmp(glacier,'Parlung4')
    catchShp_path = [root_B '\Study_sites\' glacier '\Spatial_data\Masks\Basin_Parlung4_100m_32647.shp'];
elseif strcmp(glacier,'Mugagangqiong')
    catchShp_path = [root '\Study_sites\' glacier '\Spatial_data\Masks\Basin_Mugagangqiong_100m_32645.shp'];
elseif strcmp(glacier,'Langtang')
    catchShp_path = [root_B '\Study_sites\' glacier '\Spatial_data\Shapefiles\Basin_Langtang_100m_UTM45.shp'];
end
catchShp_utm = shaperead(catchShp_path);
catchShp_info = shapeinfo(catchShp_path);
