function [imOut,imOutR] = cropraster(imIn,imInR,xyLims)

% Get pixel centres
[imInX,imInY] = pixcenters(imInR,imInR.RasterSize);

% Need to flip y values so they correspond to image
imInY = flip(imInY);

% Get rid of data to west of xyLims
toDelete = imInX < xyLims(1);
imOut = imIn(:,~toDelete);
imOutX = imInX(~toDelete);

% Get rid of data to east of xyLims
toDelete = imOutX > xyLims(2);
imOut = imOut(:,~toDelete);
imOutX = imOutX(~toDelete);

% Get rid of data to south of xyLims
toDelete = imInY < xyLims(3);
imOut = imOut(~toDelete,:);
imOutY = imInY(~toDelete);

% Get rid of data to north of xyLims
toDelete = imOutY > xyLims(4);
imOut = imOut(~toDelete,:);
imOutY = imOutY(~toDelete);

% Get x and y mins and maxs
xMin = min(min(imOutX));
xMax = max(max(imOutX));
yMin = min(min(imOutY));
yMax = max(max(imOutY));

% Make new referencing object
outRes = round(imInR.CellExtentInWorldX);
imOutSize = size(imOut);
xLims = [xMin-0.5*outRes xMax+0.5*outRes];
yLims = [yMin-0.5*outRes yMax+0.5*outRes];
imOutR = maprefcells(xLims,yLims,imOutSize);

% Note, x and y limits going out will be slightly different from those 
% going in if the pixel size is not a factor of those going in or of pixel
% centres

end
