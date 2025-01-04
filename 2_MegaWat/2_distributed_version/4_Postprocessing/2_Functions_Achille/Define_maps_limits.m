function [xmin_map, xmax_map, ymin_map, ymax_map] = Define_maps_limits(x,y,MASK)
%DEFINE_MAPS_LIMITS Summary of this function goes here
%   Detailed explanation goes here
X = repmat(x,length(y),1); X(flipud(MASK)==0) = NaN;
xmin_map = min(X,[],'all') - 1000;
xmax_map = max(X,[],'all') + 1000;
X_range = xmax_map - xmin_map;

Y = repmat(flipud(y),1,length(x)); Y(flipud(MASK)==0) = NaN;
ymin_map = min(Y,[],'all') - 1000;
ymax_map = max(Y,[],'all') + 1000;
Y_range = ymax_map - ymin_map;

if Y_range > X_range
    xmin_map = xmin_map - 0.5*(Y_range - X_range);
    xmax_map = xmax_map + 0.5*(Y_range - X_range);
elseif Y_range < X_range
    ymin_map = ymin_map - 0.5*(X_range - Y_range);
    ymax_map = ymax_map + 0.5*(X_range - Y_range);
end
end

