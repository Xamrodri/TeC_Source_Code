% ===============================================================================
%  CORRECT FOR WIND-UNDERCATCH (MASUDA ET AL., 2019)
% ===============================================================================
% Masuda et al., 2019 ESS (http://dx.doi.org/10.1029/2019EA000659)
% Tab. 2:

function [PPcorr,solid_fraction, CR] = UnderCatchCorr_Masuda(PP,U,T,RH,PRESS,ELE,Ztb,Zu,Z0,MODE,THRESHOLD)

% PP = Precipitation (mm/hr)
% U = Wind at given height (m s-1)
% T = Temperature (Deg C)
% RH = Relative Humidity (%)
% PRESS = Air Pressure (hPa)
% ELE = Elevation of point (m a.s.l.)
% Ztb = Height of tipping bucket instrument (m)
% Zu = Height of anemometer (m)
% Z0  = Relative roughness element (0-1);

% Convert air temperature to Kelvin
Tk = T + 273.15;

switch MODE
        % for CR (catch ratio) and m (correction parameter) 
    % (Yokoyama et al., 2003)
    
    case('TB') % For Unshielded RT-1 (Tipping Bucket)  

    cc_s = 0.213;  %correction coefficient snow 
    cc_r = 0.0454; %correction coefficient rain
    
    case('WG') % Weighing Gauge (with heating or anti-freeze solution) - non shielded
        
         cc_s = 0.30;  %correction coefficient snow (default:   0.346)
        cc_r = 0.0856; %correction coefficient rain (default: 0.0856)
        
    case('WGS') % Weighing Gauge (with heating or anti-freeze solution) - shielded    
        
         cc_s = 0.128;  %correction coefficient snow 
        cc_r = 0.0192; %correction coefficient rain   
        
end


% Wind speed at height of rain gauge
Ucorr = U.*(log(Ztb)-log(Z0))./(log(Zu)-log(Z0));  %in R: log() = ln()

%Type of precipitation: solid or liquid 
switch THRESHOLD
    case('STATIC')
        
        solid_fraction = Tk < 275.15;
        
    case('DING')

        [type,solid_fraction,Snow_acc,Tw] = Precip_Type(PP,Tk,RH,PRESS,ELE);

end
% Catch ratio (partition between liquid and snow)

for tt = 1:length(PP)
    if PP(tt) == 0 
        PPcorr(tt) = 0;
        CR(tt) = NaN;
        solid_fraction(tt) = NaN;
    elseif isnan(PP(tt))
        PPcorr(tt) = NaN;
        CR(tt) = NaN;
        solid_fraction(tt) = NaN;
    else 
        if solid_fraction(tt) < 0.5
                 m = cc_r;
        else
                 m = cc_s;
        end

        CR(tt) = 1/(1+m*Ucorr(tt));
        
        PPcorr(tt) = PP(tt)/CR(tt); %adjusted precipitation [mm/hr]
    end
end

PPcorr = transpose(PPcorr);
CR = 1./transpose(CR);