%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATE GENERATOR
%{
Author: Maximiliano Rodriguez
Date: October 27, 2025

Note: 
This code generates a vector from date 1 to date 2 by hour.

Variable:
    yy: year as numeric
    mth: month as numeric

Returns:
    Date_Vector: A date vector for that year and that month
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dedubgger
%yy_entry = 2000;
%mth_entry = 5;

%% FUNCTION
function [Date_Vector] =  Date_generator(yy_entry, mth_entry)

yy = char(num2str(yy_entry));
mth = char(num2str(mth_entry));

% Vector date
%--------------------------------------------------------------------------
Start_String = [yy, '-', mth, '-01 00:00:00'];
days_in_month = eomday(str2double(yy), str2double(mth));

End_String = [yy, '-', mth, '-', num2str(days_in_month), ' 23:00:00'];
input_format = 'yyyy-MM-dd HH:mm:ss';

StartTime = datetime(Start_String, 'InputFormat', input_format);
EndTime = datetime(End_String, 'InputFormat', input_format);

Date_Vector = StartTime : hours(1) : EndTime;
Date_Vector = Date_Vector';

end