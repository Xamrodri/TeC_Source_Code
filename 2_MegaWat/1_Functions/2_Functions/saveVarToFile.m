function saveVarToFile(filename, variableValue)
% SAVEVARTOFILE Saves a single variable to a MAT-file, inferring the variable name.
%
%   SAVEVARTOFILE(FILENAME, VARIABLEVALUE) saves the contents of 
%   VARIABLEVALUE to the specified FILENAME, using the variable's name
%   as the name inside the .mat file.

    % 1. Get the name of the variable passed as the second argument (position 2)
    variableNameStr = inputname(2); 

    if isempty(variableNameStr)
        % This is a safety check: inputname returns empty if the input was an expression.
        error('saveVarToFile:InvalidInput', ...
              'The second argument must be a simple variable name (not an expression or constant).');
    end

    % 2. Use the functional form of 'save' to save the variable's VALUE 
    %    and assign it the extracted NAME string.
    
    % The trick is to temporarily create a variable with the desired name in 
    % the current function's workspace, and then use the functional save form.
    eval([variableNameStr ' = variableValue;']);

    % Save the file using the extracted name string. 
    % Note: The variableNameStr must exist in the workspace when this command is run.
    save(filename, variableNameStr);
end