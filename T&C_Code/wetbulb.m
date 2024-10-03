function [wetbulbTemp] = wetbulb(P,Td,T)
%WETBULB calculates the wetbulb temperature given pressure, dewpoint, and
%temperature.
%   WETBULB(P,Td,T) calculates the wetbulb temperature given pressure,
%   dewpoint, and temperature where in input variables are:
%      P: pressure in hectopascals - scalar or matrix
%      Td: dewpoint in degrees C - scalar or matrix
%      T: temperature in degrees C - scalar or matrix
%
%   If given matrices, the sizes must match.
%
%   Examples
%      P = 1013.4;  % pressure in hPa
%      T = 13.2;  % temperature in degrees C
%      Td = 8.2;  % temperature in degrees C
%      wetbulbTemp = wetbulb(P,Td,T);
%
%      wetbulbTemp =
%          10.4637
%
%      P = 1013.4;  % pressure in hPa
%      T = linspace(10,15,20);  % 20 temperatures from 10-15 in degrees C
%      Td = 8.2;  % temperature in degrees C
%      wetbulbTemp = wetbulb(P,Td,T);
%
%References and Notes
%Function to calculate wetbulb temperature given pressure, dewpoint, and
%temperature. Uses the psychrometric formula from the American
%Meteorological Society glossary. Vapor pressure calculations use the
%improved August-Roche-Magnus approximation, that is, equation 21 from
%Alduchov, O.A. and R.E. Eskridge, 1996: 
%Improved Magnus Form Approximation of Saturation Vapor Pressure.
%J. Appl. Meteor., 35, 601-609,
%https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2
%
%Origin Date: 6/21/2018
%Last major revision: 6/22/2021 by Matthew Miller
%
%Change notes:
%6/22/2021
%- Now used fzero for dramatic speed increase
%- Now uses arrayfun to handle matrix inputs
%- Can handle combination of vector and scalar inputs
%
%Written by: Daniel Hueholt
%North Carolina State University
%Undergraduate Research Assistant at Environment Analytics
%
%See also FZERO, ARRAYFUN
%
%handle a combination of scalar and matrix inputs and return output in shape
%matching input(s)
if ~isscalar(P) || ~isscalar(Td) || ~isscalar(T)
    for ii = {P,Td,T}
        if ~isscalar(ii{1})
            vec_size = size(ii{1});
            break
        end
    end
    if isscalar(P)
        P = repmat(P,vec_size);
    end
    if isscalar(Td)
        Td = repmat(Td,vec_size);
    end
    if isscalar(T)
        T = repmat(T,vec_size);
    end
end
if sum(Td>T) > 0
    warning('Td > T. Theoretically possible but unlikely. Twb value returned may not be realistic. Check your data!');
end
wb_fun = @(Var_P,Var_Td,Var_T) wb_solver(Var_P,Var_Td,Var_T);
wetbulbTemp = arrayfun(wb_fun,P,Td,T);
end
%%%%%% LOCAL FUNCTION(S) %%%%%%
function [wetbulbTemp_single] = wb_solver(P_single,Td_single,T_single)
%if Td_single > T_single
%    wetbulbTemp_single = NaN;
%else
%epsilon = 0.622;
%Lv = 2.5*10^6; %J/kg
%Cp = 1005; %J/kg
%psychro = (1005.*P_single)./(0.622.*(2.5*10^6)); %Psychrometric constant
eAct = 6.1094.*exp((17.625.*Td_single)./(243.04+Td_single)); % Actual vapor pressure calculated from Td using improved ARM
if ~isfinite(P_single) || ~isfinite(Td_single) || ~isfinite(T_single)
    wetbulbTemp_single = NaN;
else
    if T_single > 0
        Fun = @(Tw) 6.1094.*exp((17.625.*Tw)./(243.04+Tw))-6.60.*10^(-4).*(1+0.00115.*Tw).*P_single.*(T_single-Tw) - eAct;
    else
        Fun = @(Tw) 6.1094.*exp((17.625.*Tw)./(243.04+Tw))-5.82.*10^(-4).*(1+0.00115.*Tw).*P_single.*(T_single-Tw) - eAct;
    end
    
    wetbulbTemp_single = fzero(Fun,0.667*T_single + 0.333*Td_single); %Solves the wetbulb equation numerically
end
end