function isLeap = is_leap_year(year)
% IS_LEAP_YEAR Determines if a given year is a leap year.
%   isLeap = IS_LEAP_YEAR(YEAR) returns true if YEAR is a leap year, and
%   false otherwise.

    if ~isscalar(year) || ~isnumeric(year) || year ~= fix(year) || year <= 0
        error('Input must be a positive integer scalar year.');
    end

    if mod(year, 400) == 0
        isLeap = true;
    elseif mod(year, 100) == 0
        isLeap = false;
    elseif mod(year, 4) == 0
        isLeap = true;
    else
        isLeap = false;
    end
end