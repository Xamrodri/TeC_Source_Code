%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[VI]=Interp_Smooth_Cycle_Pr(V,Date)
%%%
V=reshape(V,1,length(V));
Date=reshape(Date,1,length(V));
%%%%%%%%%% Intepolating Missing Value with replacement with the cycle of
%%%%%%%%%% that day
cur_dir=cd;
[Yr,Mo,Da,Hr,Mi]=datevec(Date');
Datam = [Yr, Mo, Da, Hr];
jDay=Date*0;
cd('D:\T&C\TC');
for i=1:length(Date);
    [jDay(i)]= JULIAN_DAY(Datam(i,:));
end
cd(cur_dir);
Hr=Hr';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ac=NaN*ones(366,24);
for i=1:366; %%%
    for k=0:23
        I= find(jDay==i & Hr==k);
        if not(isempty(I))
            Ac(i,k+1) = nanmean(V(I));
        end
    end
end
%%%%%
AcD = sum(Ac,2);
Ac2=Ac*0;
for i=1:366
    Ac2(i,19:22)=AcD(i)/4;
end
%%%%
for i=1:length(Date);
    if isnan(V(i))
        V(i)=Ac2(jDay(i),Hr(i)+1);
    end
end
VI=V;
return