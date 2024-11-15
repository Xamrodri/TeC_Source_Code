function IGM_netcdf(path_out,x,y,GLA,GLH,DEM_bedrock, DEM, SMB)

%%%%%%%% INPUT %%%%%%%%%%%%%%%%%
% path_out = absolute path of where should be stored the .nc file
% x : vector of X coordinates in UTM
% y : vector of Y coordinates in UTM
% GLA: glacier mask array
% GLH: ice thickness array
% DEM_bedrock: bedrock DEM
% DEM: Surface DEM
% SMB : surface map balance map

%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%
% Create a the "geology.nc" file needed as input for IGM

col = size(GLH,2);
row = size(GLH,1);

% Create .nc file or open it if already existing
if exist(path_out, 'file') == 0
    ncid = netcdf.create(path_out,'NETCDF4');
else 
    delete(path_out);
    ncid = netcdf.create(path_out,'NETCDF4');
end 

% Define dimensions
dimid2 =  netcdf.defDim(ncid,'y',col);
dimid1 =  netcdf.defDim(ncid,'x',row);

% Define variables
y_id = netcdf.defVar(ncid,"y","NC_DOUBLE",dimid2);
x_id = netcdf.defVar(ncid,"x","NC_DOUBLE",dimid1);
icemask_id = netcdf.defVar(ncid,"icemask","NC_DOUBLE",[dimid1 dimid2]);
topg_id = netcdf.defVar(ncid,"topg","NC_DOUBLE",[dimid1 dimid2]);
usurf_id = netcdf.defVar(ncid,"usurf","NC_DOUBLE",[dimid1 dimid2]);
thk_id = netcdf.defVar(ncid,"thk","NC_DOUBLE",[dimid1 dimid2]);
smb_id = netcdf.defVar(ncid,"smb","NC_DOUBLE",[dimid1 dimid2]);


% Add values to variables
netcdf.putVar(ncid,y_id,y)
netcdf.putVar(ncid,x_id,x)
netcdf.putVar(ncid,icemask_id,GLA)
netcdf.putVar(ncid,topg_id,DEM_bedrock)
netcdf.putVar(ncid,usurf_id,DEM)
netcdf.putVar(ncid,thk_id,GLH)
netcdf.putVar(ncid,smb_id,SMB)


% Close .nc file
netcdf.close(ncid)