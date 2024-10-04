% Can speed this up by using interp2 instead of scatteredInterpolant.
% Instead of converting lat-lon pixel centres to UTM, then interpolating 
% from these in UTM, calculate new pixel centres in lat-lon, then 
% interpolate from these in lat-lon

function [imOut,imOutR] = projectraster(imIn,imInInfo,outRes)

% Get image size
imInSize = size(imIn);

% Get referencing matrix
imInRM = imInInfo.RefMatrix;

% Get pixel locations in lat lon
[lon,lat] = pixcenters(imInRM,imInSize);
[lon,lat] = meshgrid(lon,lat);

% Get UTM zone of image centre
utmZone = utmzone(mean(lat,'omitnan'),mean(lon,'omitnan'));

% Set up UTM coordinate system
utmStruct = defaultm('utm'); 
utmStruct.zone = utmZone;  
utmStruct.geoid = wgs84Ellipsoid;
utmStruct = defaultm(utmStruct);

% Get UTM pixel locations
[x,y] = mfwdtran(utmStruct,lat(:),lon(:));
% [x,y] = ll2utm(lat(:),lon(:));
x = reshape(x,imInSize);
y = reshape(y,imInSize);

% Make UTM grid for output raster
xMax = max(max(x));
yMax = max(max(y));
xMin = min(min(x));
yMin = min(min(y));
imOutX = xMin:outRes:xMax;
imOutY = flip(yMin:outRes:yMax);
[xQ,yQ] = meshgrid(imOutX,imOutY);

% Resample the input raster to the UTM grid
imOut = griddata(x,y,imIn,xQ,yQ,'cubic');
% F = scatteredInterpolant(x(:),y(:),imIn(:),'natural','none');
% imOut = F(xQ,yQ);

% Make new referencing object
imOutSize = size(imOut);
xLims = [xMin-0.5*outRes xMax+0.5*outRes];
yLims = [yMin-0.5*outRes yMax+0.5*outRes];
imOutR = maprefcells(xLims,yLims,imOutSize);

% % Make new referencing matrix
% imOutRM = makerefmat(xMin,yMax,outRes,outRes);

end