%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute_Land_Cover_Fractions         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
function[]=Check_Land_Cover_Fractions(Crock,Curb,Cwat,Cbare,Ccrown) 
%%%INPUTS
% Csno - Crock - Curb - Cwat - Cbare - Ccrown (1...n)  
%%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if  round(Cwat + Curb + Crock + Cbare + sum(Ccrown)) ~=1
    disp('LAND USE COVER INPUTS INCONSISTENT')
    clear Ccrown Curb Crock Cbare Cwat 
    return 
end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
