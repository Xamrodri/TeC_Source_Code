function [startIdx, endIdx] = LongestSeriesOnes(vector)
    maxLen = 0;      % Maximum length of continuous 1s
    currentLen = 0;  % Current length of continuous 1s
    startIdx = 0;    % Start index of the longest series
    endIdx = 0;      % End index of the longest series
    currentStart = 0;% Current start index of the series

    for i = 1:length(vector)
        if vector(i) == 1
            if currentLen == 0
                currentStart = i; % Start of a new series of 1s
            end
            currentLen = currentLen + 1;
        else
            if currentLen > maxLen
                maxLen = currentLen;
                startIdx = currentStart;
                endIdx = i - 1;
            end
            currentLen = 0; % Reset the current length for the next series
        end
    end

    % Check the last series
    if currentLen > maxLen
        startIdx = currentStart;
        endIdx = length(vector);
    end
end