%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% CODE FOR BIAS CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  CORRECTION APPLICATION  %%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
AUTHOR: MAXIMILIANO RODRIGUEZ
DATE:07/August/2025

NOTES:
This code changes the format of the big matrices for better storage. In
a compressed format
It changes double to single or integer depending of the needs.

Variables as defined by ERA5

%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mat_mod = Data_compressor(mat, decision)

mat_mod = mat; 

if decision == "compress"

mat_mod.t2m = int16(mat.t2m*100);
mat_mod.d2m = int16(mat.d2m*100);
mat_mod.tp = int16(mat.tp*100);
mat_mod.ssrd = int16(mat.ssrd*10);
mat_mod.strd = int16(mat.strd*10);
mat_mod.ws10 = int16(mat.ws10*10);
mat_mod.sp = int32(mat.sp);
mat_mod.es = int16(mat.es*10);
mat_mod.ea = int16(mat.ea*10);
mat_mod.RH = int16(mat.RH*10);

        if ismember("SAD1", mat.Properties.VariableNames)
        mat_mod.SAD1 = int16(mat.SAD1*10);
        end
        if ismember("SAD2", mat.Properties.VariableNames)
        mat_mod.SAD2 = int16(mat.SAD2*10);
        end
        if ismember("SAB1", mat.Properties.VariableNames)
        mat_mod.SAB1 = int16(mat.SAB1*10);
        end
        if ismember("SAB2", mat.Properties.VariableNames)
        mat_mod.SAB2 = int16(mat.SAB2*10);
        end
        if ismember("PARB", mat.Properties.VariableNames)
        mat_mod.PARB = int16(mat.PARB*10);
        end
        if ismember("PARD", mat.Properties.VariableNames)
        mat_mod.PARD = int16(mat.PARD*10);
        end
        if ismember("N", mat.Properties.VariableNames)
        mat_mod.N = int16(mat.N*100);
        end



elseif decision == "decompress"

mat_mod.t2m = single(mat.t2m)/100;
mat_mod.d2m = single(mat.d2m)/100;
mat_mod.tp = single(mat.tp)/100;
mat_mod.ssrd = single(mat.ssrd)/10;
mat_mod.strd = single(mat.strd)/10;
mat_mod.ws10 = single(mat.ws10)/10;
mat_mod.sp = int32(mat.sp);
mat_mod.es = single(mat.es)/10;
mat_mod.ea = single(mat.ea)/10;
mat_mod.RH = single(mat.RH)/10;

        if ismember("SAD1", mat.Properties.VariableNames)
        mat_mod.SAD1 = single(mat.SAD1)/10;
        end
        if ismember("SAD2", mat.Properties.VariableNames)
        mat_mod.SAD2 = single(mat.SAD2)/10;
        end
        if ismember("SAB1", mat.Properties.VariableNames)
        mat_mod.SAB1 = single(mat.SAB1)/10;
        end
        if ismember("SAB2", mat.Properties.VariableNames)
        mat_mod.SAB2 = single(mat.SAB2)/10;
        end
        if ismember("PARB", mat.Properties.VariableNames)
        mat_mod.PARB = single(mat.PARB)/10;
        end
        if ismember("PARD", mat.Properties.VariableNames)
        mat_mod.PARD = single(mat.PARD)/10;
        end
        if ismember("N", mat.Properties.VariableNames)
        mat_mod.N = single(mat.N)/100;
        end

else
    disp('Choose between compress or decompress')
    error('check_var:Var', ['Error: Errors found in the selection. Choose "compress" or "decompress"']);
end



end