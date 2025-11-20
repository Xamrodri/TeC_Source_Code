%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTHOR: Maximiliano Rodriguez
% This code is for the interpolation of a map
% It consider that missing values are as specified in the missing parameter
% A matrix must be given.
%
% Inputs:
%   map: a matrix.
%   missing: The value where to interpolate. It can be 0 or NaN. 
%   plotImage: If the plot is needed for the before and after. Parameters
%   is "yes" for plotting and anything else for no plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Debugger
%A = [1 2 0 4; 5 0 7 8; 9 10 11 0];
%map = PSAN_map.Z 
%map = A;
%--------------------------------------------------------------------------

%% MAIN FUNCTION
function result=Interpolation_map(map, missing, plotImage)

%missing = 0; 
result = map; 
%imagesc(result)

% Find indices of non-missing elements
[row_non_zero, col_non_zero, val_non_zero] = find(map);

% Find missing values
[row_zero, col_zero] = find(map == missing);

% Create an interpolant from the non-zero data points
F = scatteredInterpolant(col_non_zero, row_non_zero, double(val_non_zero), 'linear');

% Use the interpolant to predict values for the zero-valued points
interpolated_values = F(col_zero, row_zero);

% Replace the zeros in the original matrix with the interpolated values
result(result == missing) = interpolated_values;

if plotImage == 'yes'
figure();
subplot(1, 2, 1);
imagesc(map);
title('Original map');
colorbar;

subplot(1, 2, 2);
imagesc(result);
title('Interpolated map');
colorbar;
end

end
