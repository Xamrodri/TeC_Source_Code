function [newDHg,newDEM] = apply_Huss_Hug(relDH,gMB,DEM,THX)

MASK=THX>0;

if nansum(THX(:))>0 
    if nansum(MASK,'all') > 1 & length(unique(DEM(MASK))) > 2 % To avoid problems with 1-pixel glaciers or glacier with constant elevation

    rho_g = 0.85; %glacier-wide specific gravity of ice at the glacier surface (see Huss et al 2013);

    %% normalize elevation and determine hypsometry
    DEM(~MASK) =NaN; 
    nDEM=(DEM-nanmin(DEM(MASK)))./(nanmax(DEM(MASK))-nanmin(DEM(MASK)));
    nDEM=(1-nDEM); %Huss2010 uses inverted normalized elevation!
    nDEM(MASK==0)=NaN;

    gArea=sum(MASK(:)); %just npixels, to normalize area!

    dEL=(max(DEM(MASK))-min(DEM(MASK)))./(length(relDH)-1); %10m interval
    ELs = min(DEM(MASK)):dEL:max(DEM(MASK));
    nELs=1-(ELs-min(DEM(MASK)))/(max(DEM(MASK))-min(DEM(MASK)));
    nArea = NaN.*nELs;

    for iel=1:numel(nELs)-1
%        cur=(nDEM>nELs(iel+1))&(nDEM<=nELs(iel)); %current section of DEM
        cur=(DEM<(ELs(iel)+dEL/2))&(DEM>=(ELs(iel)-dEL/2)); %current section of DEM
        nArea(iel)=nansum(cur(:))./gArea;
    end

    fs=gMB./nansum(relDH.*nArea);
    newDH=fs.*relDH; %altitudinal

    %% distribute spatially
    newDHv=interp1(nELs,newDH,nDEM(MASK));
    newDHg=double(MASK);newDHg(MASK)=newDHv;

    newDHg(isnan(newDHg))=0;
    newDHg(MASK==0)=0;

    % figure; 
    % subaxis(1,2,1)
    % plot(nELs,newDH)
    % subaxis(1,2,2)
    % imagesc(newDHg);colorbar
    %% deal with positive mass balance; similar to Huss and Hock 2015 but distinct routine for when a glacier advances: terminus must be 10m thick

        %determine thickness of lowest 10 pixels - this is the thickness to be added
        [~,il10]=sort(DEM(MASK));THXv=THX(MASK);ie = min(length(THXv),10);
        l10thx=nanmean(THXv(il10(1:ie)));

        if l10thx>10 & gMB > 0 %if the terminus is thicker than 10m and MB is positive, then advance. otherwise just thicken
            advVol=min(0.1.*gMB.*nansum(THX(:)>0),l10thx.*20); %volume attirbutable to advance is 1/10 the total mass gain or 5m thick for 10 adjacent pixels; note in pix-m
            gMBn=gMB-advVol./(nansum(THX(:)>0)); %remaining mass gain to distribute; note advVol is in pix-m
            newDHg=newDHg.*gMBn./gMB; %scale the DH pattern accordingly

            %identify terminus pixels as lowest-elevation 20 pixels bordering glacier
            PotentialPix = logical(imdilate(MASK,strel('diamond',1))-MASK);
            [~,iRank]=sort(DEM(PotentialPix)); %the index to the ranked values
            [~,Rank]=sort(iRank); %the rank of each value
            rankedPotentialPix=0.*DEM;rankedPotentialPix(PotentialPix)=Rank;
            Advance=(PotentialPix &(rankedPotentialPix<=20));

            newDHg(Advance)=advVol./20; %distribute volume attributable to advance 

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