function [rician_param, raw] = get_rician_parameter_from_excel_file(filename, xlsrange)
% get rician parameter from excel file
%
% [output]
% - rician_param: rician_param cell array whose element is double array
% - raw: rician_param cell array whose element is string, 
%   used in fm_show_fading_parameter_call function(uimenu callback)
% [usage]
% get_rician_parameter_from_excel_file('rician_param.xlsx')

% sheet = 1;
% xlsrange = 'b3:d7'; % see excel file(blue bold cell)

[num,txt,raw] = xlsread(filename, xlsrange);
% num,txt,raw

% ######## row meaning of raw: sensor number
% ######## column meaning of raw: 
% column 1, k(rician k factor), number(scalar) or text(vector) or empty
% column 2, tau delta ratio(path delay increment), number(scalar) or text(vector) or empty(when k is empty)
% column 3, pdb(path gain in db), number(scalar) or text(vector) or empty(when k is empty)

[sensor_length] = size(raw, 1);

% ######## row meaning of rician_param: sensor number
% ######## column meaning of rician_param: 
% column 1, k(rician k factor), scalar or vector or empty
% column 2, tau delta ratio(path delay increment), scalar or vector or empty(when k is empty)
% column 3, pdb(path gain in db), scalar or vector or empty(when k is empty)

rician_param = cell(sensor_length, 3);

for n = 1 : sensor_length
    % k factor
    k = raw{n, 1};
    if isnan(k)
        continue;
    end
    
    if isnumeric(k)
        rician_param{n, 1} = k;
    else
        rician_param{n, 1} = str2num(k);        
    end
    k_factor = rician_param{n, 1};
    
    % path delay
    tau = raw{n, 2};
    if isnumeric(tau)
        rician_param{n, 2} = tau;
    else
        rician_param{n, 2} = str2num(tau);
    end
    tau_sec = rician_param{n, 2};
    
    % path gain
    pdb = raw{n, 3};
    if isnumeric(pdb)
        rician_param{n, 3} = pdb;
    else
        rician_param{n, 3} = str2num(pdb);
    end
    pdb_db = rician_param{n, 3};
    
    k_length = length(k_factor);
    if k_length > 1
        if length(tau_sec) ~= k_length || length(pdb_db) ~= k_length
            error('###### [sensor %d] length of path delay and path gain must be same as k factor length', n);
        end
    else
        if length(tau_sec) ~= length(pdb_db)
            error('###### [sensor %d] length of path delay must be same as length of path gain', n);
        end
    end
end

rician_param;

end

%%
function [] = learn_xlsread

filename = 'rician_param.xlsx';

sheet = 1;
xlsrange = 'b3:e6';

[num,txt,raw] = xlsread(filename, sheet, xlsrange);
num,txt,raw

whos

end
