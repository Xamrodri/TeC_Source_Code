function Tw = Wetbulb_Temp( Ta, rh, pres )

%==================================================
%
%                             Wet Bulb Temperature
%
%       This scheme is to calculate the wet bulb temperature by using
%       Tetens's empirical formula of saturated vapor pressure. 
%
%       Reference: 
%           Murray, F.W., 1967. On the computation of saturation 
%           vapor pressure. J Appl Meteorol 6, 203-204.
%
%==================================================
%
%       Author: Baohong DING
%       22/11/2012
%
%==================================================

%% --------------------Arguments Announcement-----------------------

%{

    rh                         % ( input data ) relative humidity [range: 0~1]
    pres                     % ( input data ) air pressure [hPa]
    Ta                        % ( input data ) air temperature [degree C]

    Tw                       % ( output data ) wetbulb temperature [degree C]

    esat                      % ( local variable ) saturated water vapor pressure [hPa]
    desat                    % ( local variable ) 
    e                          % ( local variable ) vapor pressure [hPa]

%}

%-----------------------------------End Argument List-----------------------------------

%% ----------------------------Calculation---------------------------------

    Lv = 2.5104e6;                      % [ J/kg ]
    cp = 1004;                            % [J/(kg K)]
    ibsi = 0.622;

    B = cp * pres / Lv / ibsi;

    esat = 6.1078 * exp( 17.27 * Ta ./ (Ta + 237.3) );      % saturated vapor pressure [hPa], Tetens
    desat = esat * 17.27 * 237.3 ./ (Ta + 237.3) ./ (Ta + 237.3);

    e = esat .* rh;                                     % vapor pressure [hPa]

    Tw = Ta - ( esat - e ) ./ ( B + desat );    % ['C]

end