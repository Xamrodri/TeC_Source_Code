%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Partition Precipitation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Pr_sno,Pr_liq, Tmin, Tmax, T0, Twb,solid_fraction]=Precipitation_partition_function(Pr,Ta,ea,Pre,parameterize_phase,Ta_Ding_d,Pre_Ding_d,ea_Ding_d,Z)
%%%INPUTS
%%% Precipitation [mm/h]
%%% Air Temperature  [�C]
%%% OUTPUTS
%Pr_sno  solid precipitation [mm/h]
%Pr_liq  liquid precipitation [mm/h]
%%%%%%%%%%%%

    OPT_Pr_Part= parameterize_phase.OPT_Pr_Part; % 1 = 2-threshold, 2 = Ding, 3 = single-threshold, 4 = Pomeroy 2013
    Tmin = parameterize_phase.Tmin;
    Tmax = parameterize_phase.Tmax;
    Tconst = parameterize_phase.Tconst;

    %Use daily mean for Ding parametrization
%     Ta = Ta_Ding_d;
%     ea = ea_Ding_d;
%     Pre = Pre_Ding_d;

    % Coefficient for Pomeroy et al. 2013 parametrization
    b = 2.50286;
    c = 0.125006;
    %%%%%
switch OPT_Pr_Part
    case 1
    Pr_sno = Pr.*NaN;
    Pr_liq = Pr.*NaN;
    Pr_type = Pr.*NaN;
        %%REFERENCES %%  Wigmosta et al., 1994
        %Tmin= -1.1; Tmax = 3.3; %% General
        %Tmin = -1.1 ; Tmax = 2.5 ; %% RCW optimized
        case_1 = (Ta > Tmin) & (Ta < Tmax);
            Pr_sno(case_1) = Pr(case_1).*(Tmax - Ta(case_1))./(Tmax- Tmin);
            Pr_liq(case_1) = Pr(case_1)- Pr_sno(case_1);
            Pr_type(case_1) = 1;
        
        case_2 = Ta <= Tmin;
            Pr_sno(case_2) = Pr(case_2);
            Pr_liq(case_2)= 0;
            Pr_type(case_2) = 2;
        
        case_3 =  Ta >= Tmax;
            Pr_sno(case_3) = 0;
            Pr_liq(case_3)= Pr(case_3);
            Pr_type(case_3) = 3;

            [T0, Twb,solid_fraction] = deal(NaN(1,1));
        
    case 2
        %%%%%%
        % Ding, et , 2014.   J. Hydrol  http://dx.doi.org/10.1016/j.jhydrol.2014.03.038.
        %%%%%%%%%%%%%
        esat=611.*exp(17.27.*Ta./(237.3+Ta)); %% [Pa] Vapor pressure saturation
        U=ea./esat; %% Relative Humidity
        Laten= 1000.*(2501.3 - 2.361.*(Ta)); %%%% Latent heat vaporization/condensaition [J/kg]
        cp=1005 + ((Ta +23.15).^2)./3364; %% specific heat air  [J/kg K]
        gam=cp.*100.*Pre./(0.622.*Laten); %% [Pa/C] psycrometric constant
        del=(4098.*esat)./((237.3+Ta).^2); %% Pa/C
        %Twb = Ta - ( esat - ea )./( gam + del);    % [C] %% Wet bulb temperature
        Twb = Ta - esat.*(1-U)./(0.000643.*Pre + 6.1078.*del);
        %%%%%
        %%%%%%
        %Wetbulb temperature following Sadeghi et.al 2013 DOI: 10.1175/JTECH-D-12-00191.1
        %         if Ta > 0
        %             a = 611;  b = 17.368; c = 238.88; gam = 6.42*10^-4;
        %         else
        %             a = 611; b = 17.966; c = 247.15; gam = 5.68*10^-4;
        %         end
        %         xi = (-3*10^-7)*Ta.^3 - (10^-5)*Ta.^2 +(2*10^-5).*Ta + 4.44*10^-2;%%empirical coefficient
        %         phi = xi + gam*Pre/10;%%empirical coefficient
        %         lam = 0.0014*exp(0.027*Ta);%%empirical coefficient
        %         psi = a - gam*Pre/10*Ta-ea/1000;%%empirical coefficient
        %         Twb = (-phi+sqrt(phi^2-4*lam*psi))/(2*lam); %wetbulb temperature
        %         esat = a*exp(b*(Ta)/(Ta+c));% saturation vapor pressure Pa
        %         U=ea./esat; %% Relative Humidity
        %%%%%%%%%%%%%%%%
        
        g= 9.81; %% [m/s^2] gravity acceleration
        P_Ref= 1013.25; %% [Pa] reference pressure
        Rd =287.05; %% [J/kg K] dry air gas constant
        %%%%% Reference Elevation
%         if isempty(Z)
%           Zref= -((Ta+15)./2+273.15).*(Rd/g).*log(Pre./P_Ref); %% [m]
%         else
            Zref = Z;
%         end 
        Zref =Zref./1000; % elevation, [m] => [km]
%  Zref = 7;
        
        %%% Thresholds for discrimination precipitation
        dT = 0.215 - 0.099.*U + 1.018.*U.^2;
        dS = 2.374 - 1.634.*U;
        T0 = -5.87 - 0.1042.*Zref + 0.0885.*Zref.^2 + 16.06.*U - 9.614.*U.^2;
        %%%%%%%%%%
        Tmin = T0;
        Tmax = T0;
        %%%%%%%%%%%
        % Condition dT./dS > log(2)
        cond = dT./dS > log(2);
            Tmin(cond) = T0(cond) - dS(cond).*log(exp(dT(cond)./dS(cond)) - 2.*exp(-dT(cond)./dS(cond)));
            Tmax(cond) = 2.*T0(cond) - Tmin(cond);
            
        %%% Calculation of solid fraction of precipitation
        solid_fraction = 1./(1 + exp((Twb - T0)./dS));    % sleet
        %%% Discrimination of rain, sleet, and snow
%         %% Scalar version
%         if (Twb > Tmin) & (Twb < Tmax)
%             % sleet
%             Pr_sno = Pr.*solid_fraction;
%             Pr_liq = Pr.*(1-solid_fraction);
%         end
%         if Twb <= Tmin
%             % snow
%             Pr_sno = Pr;
%             Pr_liq= 0;
%         end
%         if Twb >= Tmax
%             % rain
%             Pr_sno = 0;
%             Pr_liq= Pr;
%         end
        %% Vector version
        
        case_1 = (Twb > Tmin) & (Twb < Tmax);
        case_2 = Twb <= Tmin;
        case_3 = Twb >= Tmax;
        case_4 = isnan(Twb);
% 
%          figure; plot(Ta); hold on; grid on; plot(Twb)
%          figure; plot(Tmin); hold on; grid on; plot(Tmax)
        
            % sleet
            Pr_sno(case_1) = Pr(case_1).*solid_fraction(case_1);
            Pr_liq(case_1) = Pr(case_1).*(1-solid_fraction(case_1));
            Pr_type(case_1) = 1;

            % snow
            Pr_sno(case_2) = Pr(case_2);
            Pr_liq(case_2) = 0;
            Pr_type(case_2) = 2;
            
            % rain
            Pr_sno(case_3) = 0;
            Pr_liq(case_3)= Pr(case_3);
            Pr_type(case_3) = 3;
            
            % NaN
            Pr_sno(case_4) = NaN;
            Pr_liq(case_4)= NaN;
            Pr_type(case_4) = NaN;

            nansum(Pr_sno)./nansum(Pr_sno + Pr_liq);
            
            
    case 3
             Pr_sno = Pr.*NaN;
             Pr_liq = Pr.*NaN;
             Pr_type = Pr.*NaN;

        case_1 = Ta < Tmin;
            Pr_sno(case_1) = Pr(case_1);
            Pr_liq(case_1)= 0;
            Pr_type(case_1) = 2;
        
        case_2 =  Ta >= Tmin;
            Pr_sno(case_2) = 0;
            Pr_liq(case_2)= Pr(case_2);
            Pr_type(case_2) = 3;
    case 4
        
        if length(Ta) == 1 % Scalar case 
          D = 2.06.*10^(-5)*((Ta+273.15)./273.15).^1.75; %% diffusivity of water vapour in air [m2s-1]
          p = (0.01801528.*ea)./(8.31441.*(Ta+273.15)); % water vapour density in the free atmosphere[kg m-3]
          lambda = 0.000063.*(Ta+273.15) + 0.00673; % thermal conductivity of air [J m-1s-1K-1]
          if Ta < 0
            Ls = 1000*(2834.1-0.29.*Ta-0.004.*Ta.^2); % latent heat of sublimation [J kg-1]
          else
            Ls = 1000*(2501-(2.361.*Ta)); % latent heat of vaporisation [J kg-1]
          end 
          fun = @(x)x - Ta-273.15 - (D./lambda).*Ls.*(p -(0.01801528.*0.611.*exp(17.27*(x-273.15)./(237.3+x-273.15)).*1000)./(8.31441.*(x)));
          Ti = fzero(fun,[230 350]); % Iterative solution methods such as the Newton–Raphston approach. Ti = hydrometeor temperature [°C]
          fr = 1./(1 + b.*c.^(Ti-273.15)); % Rainfall fraction [-]        
        else
          fr = Ta.*NaN; % Initialize fr vector
          D = 2.06.*10^(-5)*((Ta+273.15)./273.15).^1.75; %% diffusivity of water vapour in air [m2s-1]
          p = (0.01801528.*ea)./(8.31441.*(Ta+273.15)); % water vapour density in the free atmosphere[kg m-3]
          lambda = 0.000063.*(Ta+273.15) + 0.00673; % thermal conductivity of air [J m-1s-1K-1]
          Ls(Ta <0) = 1000*(2834.1-0.29.*Ta(Ta <0)-0.004.*Ta(Ta <0).^2); % latent heat of sublimation [J kg-1]
          Ls(Ta >= 0) = 1000*(2501-(2.361.*Ta(Ta >= 0))); % latent heat of vaporisation [J kg-1]
          for i = 1:length(Ta)
             if ~isnan(Ta(i))
                 fun = @(x)x - Ta(i)-273.15 - (D(i)./lambda(i)).*Ls(i).*(p(i) -(0.01801528.*0.611.*exp(17.27*(x-273.15)./(237.3+x-273.15)).*1000)./(8.31441.*(x)));
                 Ti = fzero(fun,[230 350]);
                 fr(i) = 1./(1 + b.*c.^(Ti-273.15));
             else
                 fr(i) = NaN;
             end
         end 
        end
            Pr_sno = Pr.*(1-fr); % Snowfall [mm]
            Pr_liq= Pr.*fr; % Rainfall [mm]
end
return
