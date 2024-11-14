%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Conductivity_Suction_O  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Ko,Po]=Conductivity_Suction(SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,lVG,Ks_mac,Omac,alpVGM,nVGM,lVGM,Phy1,s_SVG,bVG,O)
%%Loading 
%disp('Im in Conductivity suction')

%{
try %if load fails, it just continue as normal
    load('M:/19_ISTA/1_TC/3_Model_Source/TeC_Source_Code/Case_study/Kyzylsu_distributed/Store_par/Param_saved.mat')
    %Here we are checking if the values in the store struct were run
    %before. If so, then the values Ko and Po are returned.
    for z = 1:length(stored)
        if stored(z).SPAR==SPAR && stored(z).Ks == Ks && stored(z).Osat==Osat && stored(z).Ohy==Ohy && ...
            stored(z).L==L && stored(z).Pe == Pe && stored(z).O33==O33 && stored(z).alpVG==alpVG && ...
            stored(z).nVG==nVG && stored(z).lVG == lVG && stored(z).Ks_mac==Ks_mac && stored(z).Omac==Omac && ...
            stored(z).alpVGM==alpVGM && stored(z).nVGM == nVGM && stored(z).lVGM==lVGM && stored(z).Phy1c==Phy1 && ...
            stored(z).s_SVG==s_SVG && stored(z).bVG == bVG && stored(z).O==O
    
            Ko=stored(z).Ko;
            Po=stored(z).Po;
            disp('shortened')
            return
        end
    end
end

%}
%%Old code

%%REFERENCES %%   Saxton and Rawls 2006
mVG= 1-1./nVG;
mVGM= 1-1./nVGM;
%%%INPUTS
%%% Osat [] Saturation moisture 0 kPa
%%% L % Slope of logaritimc tension-moisture curve
%%% Pe % Tension at air antry (bubbling pressure) [kPa]
%%% Ks  % saturation conductivty [mm/h]
%%% O33 %% 33 kPa Moisture
%%% O Soil Moisture -- []
%%%dt time step [s]
%%% OUTPUTS
%%%%%% Ko  % hydraulic conductivty at O [mm/h]
%%% Po Tension at O [mm]
%Fe  Desorption rate [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch SPAR
    case 1
        %%% Van-Genuchten, 1980 Corrected 
        Se = (O-Ohy)./(Osat-Ohy);
        Phy = Phy1*101.9368; %% [mm] 
        P1 = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        P2 = -Phy*exp(Se.*-abs(bVG)); % [mm] 
        P = (Se<s_SVG).*P2 +(Se>=s_SVG).*P1; % [mm] 
        P(Se==0)=P2(Se==0); 
        Po = -P; 
        Ko= Ks.*((Se).^(lVG)).*(1-(1-(Se).^(1./mVG)).^mVG).^2; %%% [mm/h]
    case 2
        %%%%% Saxton and Rawls 1986 
        gw= 9810; %% specific weight water [N/m^3]
        B=1/L;
        A= exp(log(33)+B*log(O33)); % Coefficient of moisture tension
        %%%%%%%%%%%%%%%%%%%%
        Ko = Ks*(O/Osat)^(3+(2/L)); %%% [mm/h]
        %Ko = Ks*power((O/Osat),(3+(2/L))); %%% [mm/h]
        if O < O33
            Psi = A*(O^-B); %% [kPa]
            %Psi = A*power(O,-B); %% [kPa]
        else
            Psi =33-((O-O33)*(33-Pe)/(Osat-O33));%% [kPa]
        end
        %%%%%%%%%5
        Po =1000*1000*Psi/(gw); %%[mm]% Tension at O
        %%%%%%%%%%%%%%%
    case 3
         %%% Van-Genuchten, + Soil Structure effects 
        if Omac(1) == 0
            Se = (O-Ohy)./(Osat-Ohy);
            Po = -(1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
            h=-Po;
        else
            Se = (O-Ohy)./(Osat-Omac-Ohy);
            Se(Se>1)=1;
            h1 = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
            h1(Se==1 & Omac>0)=NaN; h1max = nanmax(h1);
            %%%%%%%%%%%%
            Se2 = (O-(Osat-Omac))./(Osat-(Osat-Omac)); Se2(Se2<0)=0;
            h2 = (1./alpVGM).*((Se2).^(-1./mVGM)-1).^(1./nVGM); %%% [mm]
            h2(Se2==0)=NaN;  h2(h2<h1max)=h1max;
            %%%%%%%%%
            h=[h1 ; h2]; h=nanmean(h);
            
            numn = length(Osat);
            for i = 1:numn
                if h(i) < - 10 && h(i) > - 500
                    hp=-10.^(0:0.1:3);  %% [mm]
                    Op = (Omac(i))./((1+abs(alpVGM(i)*hp).^nVGM(i)).^mVGM(i)) +  Ohy(i) + (Osat(i)-Omac(i)-Ohy(i))./((1+abs(alpVG(i)*hp).^nVG(i)).^mVG(i));
                    h(i)=interp1(Op,hp,O(i),'linear');
                end
            end
            Po = -h;
        end
        
        
        Kom=Ks.*(((1 - ((abs(alpVG.*h)).^(nVG-1)).*(1+(abs(alpVG.*h)).^nVG).^-mVG).^2)./((1+(abs(alpVG.*h)).^nVG).^(lVG.*mVG)));
        Kop = Ks_mac.*(((1 - ((abs(alpVGM.*h)).^(nVGM-1)).*(1+(abs(alpVGM.*h)).^nVGM).^-mVGM).^2)./((1+(abs(alpVGM.*h)).^nVGM).^(lVGM.*mVGM)));
        %%%%
        Ko =  Kom + Kop ;
    case 4
        %%% Van-Genuchten, 1980
        Se = (O-Ohy)./(Osat-Ohy);
        Po = -(1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        Ko= Ks.*((Se).^(lVG)).*(1-(1-(Se).^(1./mVG)).^mVG).^2; %%% [mm/h]
end

%%Saving
%{
if exist('stored', 'var')
    nrow = length(stored); %I calculate the next row available
    stored(nrow).SPAR = SPAR;
    stored(nrow).Ks = Ks;
    stored(nrow).Osat = Osat;
    stored(nrow).Ohy = Ohy;
    stored(nrow).L = L;
    stored(nrow).Pe = Pe;
    stored(nrow).O33 = O33;
    stored(nrow).alpVG = alpVG;
    stored(nrow).nVG = nVG;
    stored(nrow).lVG = lVG;
    stored(nrow).Ks_mac = Ks_mac;
    stored(nrow).Omac = Omac;
    stored(nrow).alpVGM = alpVGM;
    stored(nrow).nVGM = nVGM;
    stored(nrow).lVGM = lVGM;
    stored(nrow).Phy1 = Phy1;
    stored(nrow).s_SVG = s_SVG;
    stored(nrow).bVG = bVG;
    stored(nrow).O = O;
    stored(nrow).Ko = Ko;
    stored(nrow).Po = Po;
else %if the struc stored does not exist before, then I create it. 
    
    stored = struct('SPAR', {SPAR},'Ks',{Ks},'Osat',{Osat},'Ohy',{Ohy},'L', {L}, ...
                'Pe', {Pe},'O33', {O33},'alpVG', {alpVG},'nVG', {nVG},'lVG', {lVG}, ...
                'Ks_mac',{Ks_mac},'Omac', {Omac},'alpVGM', {alpVGM},'nVGM', {nVGM}, ...
                'lVGM', {lVGM},'Phy1', {Phy1},'s_SVG',{s_SVG},'bVG',{bVG},'O',{O}, ...
                'Ko',{Ko}, 'Po',{Po})
end
save('M:/19_ISTA/1_TC/3_Model_Source/TeC_Source_Code/Case_study/Kyzylsu_distributed/Store_par/Param_saved.mat','stored')
%}
return