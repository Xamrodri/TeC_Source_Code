function [se] = Snowline_Krajci(DTM,snow_pres)
%SNOWLINE_KRAJCI Summary of this function goes here
%   Detailed explanation goes here


scIm = double(snow_pres == 1); % Map of snow covered pixels
scIm(isnan(DTM)) = NaN;
nscIm = double(snow_pres == 0); % Map of non-snow covered pixels
nscIm(isnan(DTM)) = NaN;

zMin = min(DTM(:));
zMax = max(DTM(:));

possSles = zMin:50:zMax; 
nPossSles = length(possSles);
sB = nan(1,nPossSles); % Snow below
nsA = sB; % Non-snow above
for iPossSle = 1:nPossSles
    sB(iPossSle) = sum(scIm(DTM < possSles(iPossSle)),'all');
    nsA(iPossSle) = sum(nscIm(DTM >= possSles(iPossSle)),'all');
end
[~,ind] = min(sB+nsA);
se = possSles(ind); % snow line altitude

end

