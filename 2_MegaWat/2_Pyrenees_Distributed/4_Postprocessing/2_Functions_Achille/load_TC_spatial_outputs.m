fnsMap= dir([dir_tcout '\*SPATIAL*.mat']);

spatial_id = str2double(string(extractBetween([fnsMap.name],"SPATIAL_",".mat")));
[~,idx] = sort(spatial_id);
fnsMap = fnsMap(idx);
nFiles = length(fnsMap);

% Allocate memory 
date_m = NaT(nFiles,1);
[Psnow_mean, SMG_mean, SMS_mean, Pliq_mean, SMSm_mean] = deal(NaN(nFiles,1));

[SWEm_map, SWin_map, SMSm_map, SMG_map, AVA_map, QFM_neg_map, QFM_pos_map,...
    ESN_map, EICE_map, SSN_map, TA_map, ICE_map, WS_map, PSNOW_map, PRAIN_map,...
    PRECIP_map, ET_map, RH_map] = deal(NaN(size(DTM,1),size(DTM,2), nFiles));

textprogressbar('Loading T&C monthly spatial outputs:');

hr = [];
hr(1) = 0;

for iFile = 1:nFiles

    textprogressbar((iFile/nFiles)*100)

    load([dir_tcout '\' fnsMap(iFile).name],'Ta_spatial','Pr_spatial','Pr_liq_spatial','Pr_sno_spatial','Rsw_spatial',...
                                          'Rn_spatial','ESN_spatial','SSN_spatial','SND_spatial','EICE_spatial','SWE_spatial','SWE_avalanched_spatial',...
                                          'Qfm_spatial','Ws_spatial','Imelt_spatial','Smelt_spatial','ICE_spatial','T_H_spatial', 'T_L_spatial',...
                                          'EIn_H_spatial', 'EIn_L_spatial', 'EG_spatial', 'EWAT_spatial','EIn_rock_spatial','ea_spatial','Ds_spatial')

    hr(iFile+1) = string(extractBetween(fnsMap(iFile).name,"SPATIAL_",".mat"));
    hr_duration(iFile) = hr(iFile+1) - hr(iFile);
    num_hour = hr_duration(iFile); 

    cumHrs(iFile) = str2double(hr(iFile+1));    date_m(iFile) = startDate + hours(str2double(hr));
    date_m(iFile) = startDate + hours(hr(iFile+1));
    
    Psnow_mean(iFile) = mean(Pr_sno_spatial(mask),'omitnan');
    Pliq_mean(iFile) = mean(Pr_liq_spatial(mask),'omitnan');
    SMG_mean(iFile) = mean(Imelt_spatial(mask).*num_hour,'omitnan');
    SMSm_mean(iFile) = mean(Smelt_spatial(mask).*num_hour,'omitnan');
    
    %compute ET
    ET_spatial = T_H_spatial + T_L_spatial + EIn_H_spatial + EIn_L_spatial + EG_spatial + EWAT_spatial + EIn_rock_spatial + ...
        ESN_spatial + EICE_spatial - SSN_spatial;

    %Compute RH
    RH_spatial = ea_spatial./(ea_spatial+Ds_spatial);

    ESN = reshape(ESN_spatial.*num_hour,[nRows,nCols]); ESN(~MASK) = NaN;
    SSN = reshape(SSN_spatial.*num_hour,[nRows,nCols]); ESN(~MASK) = NaN;
    SMG = reshape(Imelt_spatial.*num_hour,[nRows,nCols]); SMG(~MASK) = NaN;
    SMS_m = reshape(Smelt_spatial.*num_hour,[nRows,nCols]); SMG_m(~MASK) = NaN;
    SND_m = reshape(SND_spatial,[nRows,nCols]); SND_m(~MASK) = NaN;
    TA = reshape(Ta_spatial,[nRows,nCols]); TA(~MASK) = NaN;
    WS = reshape(Ws_spatial,[nRows,nCols]); WS(~MASK) = NaN;
    PSNOW = reshape(Pr_sno_spatial,[nRows,nCols]); PSNOW(~MASK) = NaN;
    PRAIN = reshape(Pr_liq_spatial,[nRows,nCols]); PRAIN(~MASK) = NaN;
    PRECIP = reshape(Pr_spatial,[nRows,nCols]); PRECIP(~MASK) = NaN;
    EICE = reshape(EICE_spatial.*num_hour,[nRows,nCols]); EICE(~MASK) = NaN;
    ICE = reshape(ICE_spatial,[nRows,nCols]);ICE(~MASK) = NaN;
    SWE_m = reshape(SWE_spatial,[nRows,nCols]); SWE_m(~MASK) = NaN;
    AVA = reshape(SWE_avalanched_spatial.*num_hour,[nRows,nCols]); AVA(~MASK) = NaN;
    QFM = reshape(Qfm_spatial,[nRows,nCols]); QFM(~MASK) = NaN;
    QFM_neg =  QFM; QFM_neg(QFM_neg>0)=0;
    QFM_pos =  QFM; QFM_pos(QFM_pos<0)=0;
    SWin = reshape(Rsw_spatial,[nRows,nCols]); SWin(~MASK) = NaN;  
    ET = reshape(ET_spatial.*num_hour,[nRows,nCols]); ET(~MASK) = NaN;
    RH = reshape(RH_spatial,[nRows,nCols]); RH(~MASK) = NaN;

    SWEm_map(:,:,iFile) = SWE_m;
    SWin_map(:,:,iFile) = SWin;
    SMG_map(:,:,iFile) = SMG;
    SMSm_map(:,:,iFile) = SMS_m;
    SNDm_map(:,:,iFile) = SND_m;
    AVA_map(:,:,iFile) = AVA;
    QFM_neg_map(:,:,iFile) = QFM_neg;
    QFM_pos_map(:,:,iFile) = QFM_pos;
    ESN_map(:,:,iFile) = ESN;
    EICE_map(:,:,iFile) = EICE;
    ET_map(:,:,iFile) = ET;
    SSN_map(:,:,iFile) = SSN;
    TA_map(:,:,iFile) = TA;
    ICE_map(:,:,iFile) = ICE;
    WS_map(:,:,iFile) = WS;
    PSNOW_map(:,:,iFile)= PSNOW;
    PRAIN_map(:,:,iFile)= PRAIN;
    PRECIP_map(:,:,iFile)= PRECIP;
    RH_map(:,:,iFile)= RH;

end
textprogressbar('Finished')