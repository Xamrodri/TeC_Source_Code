fnsMap= dir([dir_tcout '\*SNOWMAP*.mat']);
spatial_id = str2double(string(extractBetween([fnsMap.name],"SNOWMAP_",".mat")));
[~,idx] = sort(spatial_id);
fnsMap = fnsMap(idx);
% fnsMap(end,:) = [];

nFiles = length(fnsMap);

% Allocate memory here TO DO
date_h = NaT(nFiles,1);
scas = NaN(nFiles,1);
[SMS_mean] = deal(NaN(nFiles,1));
[snow_pres, snow_depth, SWE_map, SMS_map, ROS_map, ICEd_map] = deal(NaN(size(DTM,1),size(DTM,2), nFiles));

textprogressbar('Loading T&C daily spatial outputs:');

for iFile = 1:nFiles
    textprogressbar((iFile/nFiles)*100)
    load([dir_tcout '\' fnsMap(iFile).name],'ros_spatial_daily','Smelt_spatial_daily','SND_spatial_daily','SSN_spatial_daily','Ice_spatial_daily')
    hr = string(extractBetween(fnsMap(iFile).name,"SNOWMAP_",".mat"));
    date_h(iFile) = startDate + hours(str2double(hr)-1);
    
    SMS_mean(iFile) = nanmean(Smelt_spatial_daily(mask));
  
    % Compute daily snow cover maps
    SND = reshape(SND_spatial_daily,[nRows,nCols]); SND(~MASK) = NaN;
    SWE = reshape(SND_spatial_daily.*ros_spatial_daily,[nRows,nCols]); SWE(~MASK) = NaN;
    SMS = reshape(Smelt_spatial_daily,[nRows,nCols]); SMS(~MASK) = NaN;
    ROS = reshape(ros_spatial_daily,[nRows,nCols]); ROS(~MASK) = NaN;
    ICE_d = reshape(Ice_spatial_daily,[nRows,nCols]); ROS(~MASK) = NaN;
    
    snow_pres(:,:,iFile) = SND.*1000 > 50;
    snow_depth(:,:,iFile) = SND;
    SWE_map(:,:,iFile) = SWE;
    SMS_map(:,:,iFile) = SMS;
    ROS_map(:,:,iFile) = ROS;
    ICEd_map(:,:,iFile) = ICE_d;
    scas(iFile) = (sum(snow_pres(:,:,iFile),'all'))/nCatchPix; % Daily snow cover fraction
end
textprogressbar('Finished')
