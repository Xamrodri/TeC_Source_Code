function rowCol = getrowcolfromdem(demR,pixelSize,xy)

% Get row and column relative to DEM
rowCol(1) = floor((demR.YWorldLimits(1,2)+0.5*pixelSize-xy(2))/...
    pixelSize)+1;
rowCol(2) = floor((xy(1)-(demR.XWorldLimits(1,1)-0.5*pixelSize))/...
    pixelSize)+1;

end