% Create debris cover mask for satellite image

% Notes
% - Script by Mike McCarthy (2018)
% - Uses debris shapefile (debShp) and satellite image info (imgInfo)

function [mask,maskIndex] = createmask(shp,imgR)

% Get shapefile into inpolygons format
rx = {shp.X};
rx = horzcat(rx{:});
ry = {shp.Y};
ry = horzcat(ry{:});

% Get image size
imgSize = imgR.RasterSize;

% Get satellite image pixel centers
[xPix,yPix] = pixcenters(imgR,imgSize);
yPix = flip(yPix);
[xPixQ,yPixQ] = meshgrid(xPix,yPix);

% Create debris cover mask where 1 is debris, 0 is not debris
[mask,maskIndex] = inpolygons(xPixQ,yPixQ,rx,ry);

end