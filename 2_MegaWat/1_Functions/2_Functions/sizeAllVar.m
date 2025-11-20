function sizeAllVar(vars, text)

%% Memory

% Iterate through the variables and sum their sizes
%--------------------------------------------------------------------------
totalSize = 0;
for i = 1:length(vars)
    totalSize = totalSize + vars(i).bytes;
end

% Display variables size
%--------------------------------------------------------------------------
disp(['Memory used by variables at ' char(text) ': ', num2str(totalSize/1e6), ' MB']);

end