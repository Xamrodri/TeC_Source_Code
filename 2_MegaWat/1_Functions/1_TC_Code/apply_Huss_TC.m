function [THX,GLA_ID] = apply_Huss_TC(TH,THym1,DEM,GLA_ID,Gla_nEl_nDH,choice,resol,period)
%%% TH -  Ice thickness current timestep [m]
%%% THym1 - Ice thickness one year before [m]
%%% DTM - digital terrain model [m asl]
%%% GLA_ID - glacier ID based on RGI inventory
%%% mask_update - updated outlines from previous years
%%% Gla_nEl_nDH - normalized elevation/dh pattern
%%% Huss_curve - choice of 'small', 'med', 'large' glacier parameters



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters to use in case of glacier advance
P1 = 10; % P1 is N terminus pixels to average the ice thickness over - best guess is 40
P2 = 0.1; % P2 is partition of mass gain that is used for glacier advance (rest is used for glacier thickening) 
P3 = 5; % P3 is the number of pixels to advance into - depends on pixel size and glacier width
P4 = 10; %P4 is the minimum terminus thickness allowing advance - guess is 10m
% Elevation interval
%dEL=10; %10m interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
glacier_id_all = unique(GLA_ID); glacier_id_all(glacier_id_all == 0 | isnan(glacier_id_all)) = [];
THX = THym1;  %distributed initial ice thickness

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% compute mass balance and apply Huss parameterization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start the loop for individual glaciers
for i = 1:length(glacier_id_all)

 
mask = GLA_ID == glacier_id_all(i) & ~isnan(DEM); % initial mask of each individual glacier
param_idx=find([Gla_nEl_nDH.Glacier_id]==glacier_id_all(i));

if sum(mask,"all")==0 || any(isnan(Gla_nEl_nDH(param_idx).(period)))
    continue
end



%nELs = Gla_nEl_nDH(i).nELs; % Normalized elevation
mMB =  nanmean(TH(mask)-THym1(mask)); % T&C OUTPUT HERE,; % annual glacier-wide MB (in m w.e. per year)
THX_i = THX; THX_i(~mask) = 0;  % ice thickness for the glacier i

gla_pixs = length(mask>1);
gla_area = resol^2*gla_pixs/10^6;
if choice == 1 % Normalized ice thickness change, standard parameters
        % standard parameters based on glacier size    
        if gla_area <= 5
            param = [-.3,.6,.09,2];
        elseif gla_area>5 && gla_area<20 
            param = [-.05,.19,.01,4];
        elseif gla_area >=20
            param = [-.02,.12,0,6];
        end    
        dEL=25;
        ELs = nanmin(DEM(mask)):dEL:nanmax(DEM(mask));
        nELs = 1-(ELs-nanmin(DEM(mask)))/(nanmax(DEM(mask))-nanmin(DEM(mask)));
        relDH = (nELs+param(1)).^param(4)+param(2).*(nELs+param(1))+param(3);

    elseif  choice == 2  % Normalized ice thickness change from dh/dt derived parameters
        dEL=(max(DEM(mask))-min(DEM(mask)))./(length(Gla_nEl_nDH(param_idx).(period))-1);
        ELs = nanmin(DEM(mask)):dEL:nanmax(DEM(mask));
        nELs = 1-(ELs-nanmin(DEM(mask)))/(nanmax(DEM(mask))-nanmin(DEM(mask)));
        relDH = Gla_nEl_nDH(param_idx).(period); 
    else 
        fprintf('Huss parameters not defined'); 
end 


[newDHg,newDEM] = apply_Huss(relDH,dEL,nELs,ELs,mMB,DEM,THX_i,P1,P2,P3,P4);

% The following lines are needed to deal with ice disappearance

    %THX_i(THX_i <= 0)=0;     
    cTHXt=THX_i+newDHg;  
    t1=cTHXt<0;  % Identify areas where SMB is more negative that ice thickness left
    newDHg(t1)=-THX_i(t1); %Set the maximum SMB to the ice thickness left

%     cSMB = gMB_all; cSMB(~mask) = NaN;
%     t2=cSMB>0; 
%     newDHg(t2&(cTHXt<=0))=0; % Orginally made to ensure that positive SMB areas remain glacier, with thickness of xx m..
%     %but here set to 0 (the glacier can retreat even in zones of positive SMB).
%     newDHg=newDHg.*mask;

    %rescale cDH to conserve mass (loss) after above corrections
    % This is not strictly mass conversative, but it also makes sense since
    % the disappareance of ice reduces the full mass loss theorically caused by SMB
    % (sometimes the ice would disappear before the end of the timestep)
    
   % newDHg(~(t1|t2))=newDHg(~(t1|t2)).*mMB./nanmean(newDHg(mask&~(t1|t2)));
    newDHg(~(t1))=newDHg(~(t1)).*mMB./nanmean(newDHg(mask&~(t1)));
    
    newDHg(isnan(newDHg)) = 0;  
    newDHg(isinf(newDHg)) = 0;
%     nansum(newDHg(mask))./.9
    NewTHX = THX_i+ newDHg;
    Newmask = NewTHX > 0;
    GLA_ID(Newmask) = glacier_id_all(i); % Updated GLA_ID map
    GLA_ID(GLA_ID == glacier_id_all(i) & ~Newmask & ~isnan(DEM)) = 0; % Updated GLA_ID map
%     GLA_ID(isnan(DEM))=NaN; %trim to catchment/domain mask
    THX(Newmask | mask) = NewTHX(Newmask | mask);  % Updated ice thickness
    DEM(Newmask | mask) = DEM(Newmask | mask)+newDHg(Newmask | mask);   % Updated surface DEM
end
end
%% actual redistribution
function [newDHg,newDEM] = apply_Huss(relDH,dEL,nELs,ELs,gMB,DEM,THX,P1,P2,P3,P4)

MASK=THX>0;

if nansum(THX(:))>0 
    if nansum(MASK,'all') > 1 && length(unique(DEM(MASK))) > 2 % To avoid problems with 1-pixel glaciers or glacier with constant elevation

%     rho_g = 0.85; %glacier-wide specific gravity of ice at the glacier surface (see Huss et al 2013);

    % normalize elevation and determine hypsometry
    aDEM=DEM;
    DEM(~MASK) =NaN; 
    nDEM=(DEM-nanmin(DEM(MASK)))./(nanmax(DEM(MASK))-nanmin(DEM(MASK)));
    nDEM=(1-nDEM); %Huss2010 uses inverted normalized elevation!
    nDEM(MASK==0)=NaN;

    gArea=sum(MASK(:)); %just npixels, to normalize area!

    %dEL=(max(DEM(MASK))-min(DEM(MASK)))./(length(relDH)-1); %10m interval
    %ELs = min(DEM(MASK)):dEL:max(DEM(MASK));
    %nELs=1-(ELs-min(DEM(MASK)))/(max(DEM(MASK))-min(DEM(MASK)));
    nArea = NaN.*nELs;

    for iel=1:numel(nELs)-1
%        cur=(nDEM>nELs(iel+1))&(nDEM<=nELs(iel)); %current section of DEM
        cur=(DEM<(ELs(iel)+dEL/2))&(DEM>=(ELs(iel)-dEL/2)); %current section of DEM
        nArea(iel)=nansum(cur(:))./gArea;
    end

    fs=gMB./nansum(relDH.*nArea);
    newDH=fs.*relDH; %altitudinal

    % distribute spatially
    newDHv=interp1(nELs,newDH,nDEM(MASK));
    newDHg=double(MASK);newDHg(MASK)=newDHv;

    newDHg(isnan(newDHg))=0;
    newDHg(MASK==0)=0;

%      figure; imagesc(nDEM);colorbar
%      figure; imagesc(newDHg);colorbar

    % deal with positive mass balance; similar to Huss and Hock 2015 but distinct routine for when a glacier advances: terminus must be 10m thick

     %determine thickness of lowest 10 pixels - this is the thickness to be added
     [~,il10]=sort(DEM(MASK));THXv=THX(MASK);ie = min(length(THXv),P1);
     l10thx=nanmean(THXv(il10(1:ie)));

        if l10thx>P4 && gMB > 0  %if the terminus is thicker than 10m and MB is positive, then advance. otherwise just thicken
            advVol=min(P2.*gMB.*nansum(THX(:)>0),l10thx.*P3); %volume attirbutable to advance is 1/10 the total mass gain or 5m thick for 10 adjacent pixels; note in pix-m
            gMBn=gMB-advVol./(nansum(THX(:)>0)); %remaining mass gain to distribute; note advVol is in pix-m
            newDHg=newDHg.*gMBn./gMB; %scale the DH pattern accordingly

            %identify terminus pixels as lowest-elevation 20 pixels bordering glacier
            PotentialPix = logical(imdilate(MASK,strel('diamond',1))-MASK);
            [~,iRank]=sort(aDEM(PotentialPix)); %the index to the ranked values
            [~,Rank]=sort(iRank); %the rank of each value
            rankedPotentialPix=0.*aDEM;rankedPotentialPix(PotentialPix)=Rank;
            Advance=(PotentialPix &(rankedPotentialPix<=P3));

            newDHg(Advance)=advVol./P3; %distribute volume attributable to advance 

        end

    newDEM=DEM+newDHg;
    else
      newDHg = 0.*THX;  % Deals with 1-pixel glaciers
      newDHg(MASK) = gMB;
      newDEM=DEM+newDHg;
    end
else
    newDEM=DEM;
    newDHg=0.*THX; 
        
end
end
