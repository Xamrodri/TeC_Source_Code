%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% CHECK OF VARIABLE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{

AUTHOR: MAXIMILIANO RODRIGUEZ
DATE: JULY/2025

Names of Variables must be inserted in this order:
    Temp = Temperature
    Dew_Temp = Dew Point temperature
    Prep = Precipitation
    SW = Incident short wave radiation
    LW = Downdward Long wave radiation
    WS = Wind speed
    Pres = Air pressure
    es = saturation vapor pressure
    ea = actual vapor pressure
    RH = Relative humidity
    SAD1 = SAD1
    SAD2 = SAD2
    SAB1 = SAB1
    SAB2 = SAB2
    PARB = PARB
    PARD = PARD
    N = Clouds

%}

%% DEBUGGER

%               1     2    3     4      5      6     7    8    9    10         
%forc_names = ["t2m" "d2m" "tp" "ssrd" "strd" "ws10" "sp" "es" "ea" "RH" ...
%"SAD1" "SAD2" "SAB1" "SAB2" "PARB" "PARD" "N"]
% 11     12     13     14    15      16  17  

%% FUNCTION
function check_var(forcing, var_names, Point)

if ~(length(var_names) == 17 ) 
disp("You have to add [Temp Dew_Temp Prep SW LW WS Pres es ea RH SAD1 SAD2 SAB1 SAB2 PARB PARD N]")
disp("Just this number of variables")
error('check_var:Variables', 'Error: You need to add 17 variables');
end

%% CHECK FOR NANS
errors = 0;
mistakes = cell(1, 1); 
warnings_code = cell(1, 1); 
z = 1; 
z2 = 1; 

if sum(isnan(forcing.(var_names(1)))) > 0   
errors = errors + 1;    
mistakes{z,1} = "Error: There are NaNs in your forcings for Temperature";
z =z+1;
end

if sum(isnan(forcing.(var_names(2)))) > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for Dew Point Temperature";
z =z+1;
end

if sum(isnan(forcing.(var_names(3)))) > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for Precipitation";
z =z+1;
end

if sum(isnan(forcing.(var_names(4)))) > 0   
mistakes{z,1} = "Error: There are NaNs in your forcings for Short Wave radiation";
z =z+1;
end

if sum(isnan(forcing.(var_names(5)))) > 0  
mistakes{z,1} = "Error: There are NaNs in your forcings for Long wave radiation";
z =z+1;   
end

if sum(isnan(forcing.(var_names(6)))) > 0 
mistakes{z,1} = "Error: There are NaNs in your forcings for Wind speed";
z =z+1;  
end

if sum(isnan(forcing.(var_names(7)))) > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for Air pressure";
z =z+1;      
end

if sum(isnan(forcing.(var_names(8)))) > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for saturation vapor pressure";
z =z+1; 
end

if sum(isnan(forcing.(var_names(9)))) > 0 
mistakes{z,1} = "Error: There are NaNs in your forcings for actual vapor pressure";
z =z+1;   
end

if sum(isnan(forcing.(var_names(10)))) > 0  
mistakes{z,1} = "Error: There are NaNs in your forcings for relative humidity";
z =z+1; 
end

if sum(isnan(forcing.(var_names(11)))) > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for SAD1";
z =z+1; 
end

if sum(isnan(forcing.(var_names(12)))) > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for SAD2";
z =z+1; 
end

if sum(isnan(forcing.(var_names(13)))) > 0   
mistakes{z,1} = "Error: There are NaNs in your forcings for SAB1";
z =z+1; 
end

if sum(isnan(forcing.(var_names(14)))) > 0  
mistakes{z,1} = "Error: There are NaNs in your forcings for SAB2";
z =z+1; 
end

if sum(isnan(forcing.(var_names(15)))) > 0    
mistakes{z,1} = "Error: There are NaNs in your forcings for PARB";
z =z+1; 
end

if sum(isnan(forcing.(var_names(16)))) > 0 
mistakes{z,1} = "Error: There are NaNs in your forcings for PARD";
z =z+1; 
end

if sum(isnan(forcing.(var_names(17)))) > 0 
mistakes{z,1} = "Error: There are NaNs in your forcings for N";
z =z+1; 
end

%% CHECK TEMPERATURES

Tmin = min(forcing.(var_names(1)));
Tmax = max(forcing.(var_names(1)));

TDPmin = min(forcing.(var_names(2)));
TDPmax = max(forcing.(var_names(2)));

if Tmin < -100
mistakes{z,1} = "Error: Temperature seems to not be in celsius degrees. There values of temperature lower than 100 celsius";
z =z+1; 
end

if Tmax > 100
mistakes{z,1} = "Error: Temperature seems to not be in celsius degrees. There values of temperature higher than 100 celsius";
z =z+1; 
end

if TDPmin < -100 
mistakes{z,1} = "Error: Dew Point Temperature seems to not be in celsius degrees. There values of Dew Point temperature lower than 100 celsius";
z =z+1; 
end

if TDPmax > 100
mistakes{z,1} = "Error: Dew Point Temperature seems to not be in celsius degrees. There values of Dew Point temperature higher than 100 celsius";
z =z+1; 
end

if any(forcing.(var_names(2)) > forcing.(var_names(1)))
mistakes{z,1} = "Error: Some values of Dew Point Temperature are higher than Air Temperature. Check.";
z =z+1; 
end


%% CHECK PRECIPITATION

if any(forcing.(var_names(3)) > 100)
mistakes{z,1} = "Error: Precipitation seems to be very high. More than 100 mm in a single hour or Inf";
z =z+1;
end

if any(forcing.(var_names(3)) < 0)
mistakes{z,1} = "Error: Precipitation is negative";
z =z+1;
end

%% CHECK solar radiation

if any(forcing.(var_names(4)) > 3000)
mistakes{z,1} = "Error: Downward shortwave solar radiation is over 3000 W m-2.";
z =z+1;
end

if any(forcing.(var_names(4)) < 0)
mistakes{z,1} = "Error: Downward shortwave solar radiation is below zero.";
z =z+1;
end

if any(forcing.(var_names(5)) > 1000)
mistakes{z,1} = "Error: Long wave radiation is higher than 1000 W m-2, check. Downward longwave radiation is usually between 200 to 450 W m-2";
z =z+1;
end

%% CHECK SOLAR COMPONENTS

if any(forcing.(var_names(11)) > 1000)
mistakes{z,1} = "Error: SAD1 from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(forcing.(var_names(12)) > 1000)
mistakes{z,1} = "Error: SAD2 from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(forcing.(var_names(13)) > 1000)
mistakes{z,1} = "Error: SAB1 from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(forcing.(var_names(14)) > 1000)
mistakes{z,1} = "Error: SAB2 from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(forcing.(var_names(15)) > 1000)
mistakes{z,1} = "Error: PARB from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end

if any(forcing.(var_names(16)) > 1000)
mistakes{z,1} = "Error: PARD from the radiation partition has values higher than 1000 W m-2 or values are Inf";
z =z+1;
end


%Less than zero
if any(forcing.(var_names(11)) < 0)
warnings_code{z,1} = "Warning: SAD1 from the radiation partition has negative values";
z =z+1;
end

if any(forcing.(var_names(12)) < 0)
warnings_code{z,1} = "Warning: SAD2 from the radiation partition has negative values";
z =z+1;
end

if any(forcing.(var_names(13)) < 0)
warnings_code{z,1} = "Warning: SAB1 from the radiation partition has negative values";
z =z+1;
end

if any(forcing.(var_names(14)) < 0)
warnings_code{z,1} = "Warning: SAB2 from the radiation partition has negative values";
z =z+1;
end

if any(forcing.(var_names(15)) < 0)
warnings_code{z,1} = "Warning: PARB from the radiation partition has negative values";
z =z+1;
end

if any(forcing.(var_names(16)) < 0)
warnings_code{z,1} = "Warning: PARD from the radiation partition has negative values";
z =z+1;
end


%% Vapor pressure
if any(forcing.(var_names(8)) > 10000)
mistakes{z,1} = "Error: Vapor pressure must be in Pa. There are values over 10000 Pa in es." + ...
    " A Temperature over 50 C is needed to get es>10000 Pa. Typical values are within 600 Pa and 5000 Pa";
z =z+1;
end

if any(forcing.(var_names(9)) > 10000)
mistakes{z,1} = "Error: Vapor pressure must be in Pa. There are values over 10000 Pa in ea. " + ...
    "A Temperature over 50 C is needed to get ea>10000 Pa. Typical values are within 600 Pa and 5000 Pa";
z =z+1;
end

%% Wind speed
if any(forcing.(var_names(6)) > 100)
mistakes{z,1} = "Error: There are wind speeds over 100 m/s. Wind speed must be in m/s";
z =z+1;
end

if any(forcing.(var_names(6)) < 0)
mistakes{z,1} = "Error: There are negative values. Wind speed must be in m/s";
z =z+1;
end

%% Relative humidity
if any(forcing.(var_names(10)) > 100)
mistakes{z,1} = "Error: Relative humidity must be between 0 and 100. Relative humidity is higher than 100%";
z =z+1;
end

%if (forcing.(var_names(10)) < 0) 
%mistakes{z,1} = "Error: Relative humidity must be between 0 and 100. Relative humidity is lower than 0%";
%z =z+1;
%end

%if any(forcing.(var_names(10)) > 110)
%warnings_code{z2,1} = "Warning: Relative humidity is higher than 110% in some days, but less than 200%";
%z2 =z2+1;
%end

%% FINAL WARNINGS
if length(warnings_code{1}) > 0
cellfun(@disp, warnings_code);
warning('check_var:Var', ['Warnings: Some data is out of range: ' char(Point) '. Please check details']);
end

%% FINAL ERROR
if length(mistakes{1}) > 0
cellfun(@disp, mistakes);
error('check_var:Var', ['Error: Errors found in the forcings for: ' char(Point) '. Please check details']);
else
disp(['No errors found in your forcings for ' char(Point)])    
end


end % end function

