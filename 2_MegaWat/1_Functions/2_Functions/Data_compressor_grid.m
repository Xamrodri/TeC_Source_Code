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

function mat_mod = Data_compressor_grid(mat, var, decision)

if decision == "compress"
    
    switch var
        case {"t2m", "d2m", "tp", "es", "ea", "RH", "N"}
            mat_mod = int16(floor(mat*100));

        case {"ssrd", "strd", "ws10", "SAD1", "SAD2", "SAB1", "SAB2", "PARB", "PARD"}
            mat_mod = int16(mat*10);
        case "sp"
            mat_mod = int32(mat);
        otherwise
            disp("Choose the right variable")
            error('check_var:Var', 'Error: Errors found in the selection.');
    end

elseif decision == "decompress"


    switch var
        case {"t2m", "d2m", "tp", "es", "ea", "RH", "N"}
            mat_mod = single(mat)/100;

        case {"ssrd", "strd", "ws10", "SAD1", "SAD2", "SAB1", "SAB2", "PARB", "PARD"}
            mat_mod = single(mat)/10;

        case "sp"
            mat_mod = single(mat);

        otherwise
            disp("Choose the right variable")
            error('check_var:Var', 'Error: Errors found in the selection.');
    end

else
    disp('Choose between compress or decompress')
    error('check_var:Var', 'Error: Errors found in the selection. Choose "compress" or "decompress"');
end



end