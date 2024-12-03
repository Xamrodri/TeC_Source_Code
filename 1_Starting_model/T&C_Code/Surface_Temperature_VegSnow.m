%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CLOSURE ENERGY BUDGET  DETERMINATION SURFACE TEMPERATURE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[DQ]=Surface_Temperature_VegSnow(Ts,Ts2,dt,Ta,ea,Latm,SvF,Pre,...
    Csno,Ccrown,...
    hc_H,hc_L,SnoDep,LAI_H,LAI_L,SAI_H,SAI_L,...
    RabsbSun_vegH,RabsbShd_vegH,FsunH,FshdH,...
    FsunL,FshdL,e_sno,e_gr,...
    dw_H,dw_SNO,In_H,...
    rs_sunH,rs_shdH,d_leaf_H,d_leaf_L,...
    zatm,disp_h,zom,zoh,zom_under,disp_h_H,zom_H,disp_h_L,zom_L,Ws,Lpho,Vavail_plant_H)
%%%%%%%%%%%%%%%%%%%%%
%%% Ts -  Vegetation 
%%% Ts2 - surface (snow)
%%%%%%%
%%%%
[ra]=Aerodynamic_Resistence(Ta,Ts,Pre,zatm,disp_h,zom,zoh,Ws,ea);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rb_H=zeros(1,length(Ccrown)); rb_L=zeros(1,length(Ccrown));
%rap_H=zeros(1,length(Ccrown)); rap_L=zeros(1,length(Ccrown));
%for i=1:length(Ccrown)
    %%%%%%% Undercanopy and Leaf  resitence
%    if  (hc_L(i) == 0) && (hc_H(i) == 0)
%        rap_H(i) = 0; rap_L(i) = 0 ;
%        rb_H(i)=0; rb_L(i)=0;
%    else
%        [rap_H(i),rap_L(i),rb_H(i),rb_L(i)]=Undercanopy_Leaf_Resistence(Ws,Ta,Ts,hc_H(i),hc_L(i),...
%            (LAI_H(i)+SAI_H(i)),(LAI_L(i)+SAI_L(i)),d_leaf_H(i),d_leaf_L(i),...
%            zatm,disp_h,zom,zom_under,SnoDep,disp_h_H(i),zom_H(i),disp_h_L(i),zom_L(i));
%    end
%end
[rap_H,rap_L,rb_H,rb_L]=Undercanopy_Leaf_Resistence2(Ws,Ta,Ts,Ccrown,hc_H,hc_L,...
    (LAI_H+SAI_H),(LAI_L+SAI_L),d_leaf_H,d_leaf_L,...
    zatm,disp_h,zom,zom_under,SnoDep,disp_h_H,zom_H,disp_h_L,zom_L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Net Radiation
[Rnet]=Net_Radiation_Manager_VegSnow(Ts,Ts2,Latm,SvF,...
Csno,Ccrown,hc_L,SnoDep,LAI_H,LAI_L,SAI_H,SAI_L,...
    RabsbSun_vegH,RabsbShd_vegH,FsunH,FshdH,...
    FsunL,FshdL,e_sno,e_gr) ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Heat Fluxes
[H,QE]=Heat_fluxes_VegSnow(dt,...
    Ta,Ts,ea,Pre,Csno,...
    dw_H,dw_SNO,Ccrown,FsunH,FshdH,...
    LAI_H,SAI_H,In_H,...
    ra,rs_sunH,rs_shdH,rb_H,Vavail_plant_H);
%%%%%%% Energetic Balance
DQ = Rnet-H-QE-Lpho;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return