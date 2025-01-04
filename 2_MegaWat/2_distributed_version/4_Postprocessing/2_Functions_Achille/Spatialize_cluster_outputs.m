

parpool('threads')
tic
for ii = 1:length(VAR_map)    
    [MAPSi{ii}, DATESi{ii}] = Spatialize_cluster(Cluster_map,table_cluster,dt_spatialize{ii},...
        var_to_spatialize{ii});
end
toc

for ii = 1:length(VAR_map)
    MAPS.(VAR_map{ii}) = MAPSi{ii};
    DATES.(VAR_map{ii}) = DATESi{ii};
end
clear MAPSi DATESi

if ~isempty(find(strcmp(fieldnames(MAPS),{'ESN_map'}),1))
    if strcmp(dt_spatialize(strcmp(var_to_spatialize,'ESN_poi')),'monthly')
        for ii = 1:size(MAPS.ESN_map,3)
           MAPS.ESN_map(:,:,ii) = MAPS.ESN_map(:,:,ii).*eomday(year(DATES.ESN_map(ii)),month(DATES.ESN_map(ii))).*24;
        end 
    elseif strcmp(dt_spatialize(strcmp(var_to_spatialize,'ESN_poi')),'daily')
           MAPS.ESN_map = MAPS.ESN_map.*24;
    end
    ESN_map = MAPS.ESN_map;
end

if ~isempty(find(strcmp(fieldnames(MAPS),{'EICE_map'}),1))
    if strcmp(dt_spatialize(strcmp(var_to_spatialize,'EICE')),'monthly')
        for ii = 1:size(MAPS.EICE_map,3)
           MAPS.EICE_map(:,:,ii) = MAPS.EICE_map(:,:,ii).*eomday(year(DATES.EICE_map(ii)),month(DATES.EICE_map(ii))).*24;
        end 
    elseif strcmp(dt_spatialize(strcmp(var_to_spatialize,'EICE')),'daily')
           MAPS.EICE_map = MAPS.EICE_map.*24;
    end
    EICE_map = MAPS.EICE_map;
end

if ~isempty(find(strcmp(fieldnames(MAPS),{'SSN_map'}),1))
    if strcmp(dt_spatialize(strcmp(var_to_spatialize,'SSN_poi')),'monthly')
        for ii = 1:size(MAPS.SSN_map,3)
           MAPS.SSN_map(:,:,ii) = MAPS.SSN_map(:,:,ii).*eomday(year(DATES.SSN_map(ii)),month(DATES.SSN_map(ii))).*24;
        end 
    elseif strcmp(dt_spatialize(strcmp(var_to_spatialize,'SSN_poi')),'daily')
           MAPS.SSN_map = MAPS.SSN_map.*24;
    end
    SSN_map = MAPS.SSN_map;
end

if ~isempty(find(strcmp(fieldnames(MAPS),{'SMG_map'}),1))
    if strcmp(dt_spatialize(strcmp(var_to_spatialize,'Imelt')),'monthly')
        for ii = 1:size(MAPS.SSN_map,3)
           MAPS.SMG_map(:,:,ii) = MAPS.SMG_map(:,:,ii).*eomday(year(DATES.SMG_map(ii)),month(DATES.SMG_map(ii))).*24;
        end 
    elseif strcmp(dt_spatialize(strcmp(var_to_spatialize,'Imelt')),'daily')
           MAPS.SMG_map = MAPS.SSN_map.*24;
    end
    SMG_map = MAPS.SMG_map;
end

if ~isempty(find(strcmp(fieldnames(MAPS),{'SMS_map'}),1))
    if strcmp(dt_spatialize(strcmp(var_to_spatialize,'Smelt_poi')),'monthly')
        for ii = 1:size(MAPS.SMS_map,3)
           MAPS.SMS_map(:,:,ii) = MAPS.SMS_map(:,:,ii).*eomday(year(DATES.SMS_map(ii)),month(DATES.SMS_map(ii))).*24;
        end 
    elseif strcmp(dt_spatialize(strcmp(var_to_spatialize,'Smelt_poi')),'daily')
           MAPS.SMS_map = MAPS.SMS_map.*24;
    end
    SMS_map = MAPS.SMS_map;
end

if ~isempty(find(strcmp(fieldnames(MAPS),{'snow_depth'}),1))
    snow_depth = MAPS.snow_depth;
    snow_pres = snow_depth.*1000 > 50;
    scas = (squeeze(sum(snow_pres,[1 2])))/nCatchPix; % Daily snow cover fraction
end 

if ~isempty(find(strcmp(fieldnames(MAPS),{'SWE_map'}),1))
    SWE_map = MAPS.SWE_map;
end 

if ~isempty(find(strcmp(fieldnames(MAPS),{'TA_map'}),1))
    TA_map = MAPS.TA_map;
end 

if ~isempty(find(strcmp(fieldnames(MAPS),{'WS_map'}),1))
    WS_map = MAPS.WS_map;
end 

if ~isempty(find(strcmp(fieldnames(MAPS),{'RH_map'}),1))
    RH_map = MAPS.RH_map;
end 

if ~isempty(find(strcmp(fieldnames(MAPS),{'PRECIP_map'}),1))
    PRECIP_map = MAPS.PRECIP_map;
end 

if ~isempty(find(strcmp(fieldnames(MAPS),{'PSNOW_map'}),1))
    PSNOW_map = MAPS.PSNOW_map;
end 


if ~isempty(find(strcmp(fieldnames(MAPS),{'ET_map'}),1))
    if strcmp(dt_spatialize(strcmp(var_to_spatialize,{'ET'})),'monthly')
        for ii = 1:size(MAPS.ET_map,3)
           MAPS.ET_map(:,:,ii) = MAPS.ET_map(:,:,ii).*eomday(year(DATES.ET_map(ii)),month(DATES.ET_map(ii))).*24;
        end 
    elseif strcmp(dt_spatialize(strcmp(var_to_spatialize,'ET')),'daily')
           MAPS.ET_map = MAPS.ET_map.*24;
    end
    ET_map = MAPS.ET_map;
end

clear MAPS