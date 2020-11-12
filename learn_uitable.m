function [rician_param, raw] = learn_uitable(filename, xlsrange)
% learn 2-d graphic table ui component
%
% [usage]
% learn_uitable('rician_all_k,tau_3vec.xlsx', 'b3:d7');
% learn_uitable('rician_part_k=10.xlsx', 'b3:d7');

% t = uitable;
%
% set(t,'Data',magic(10));

% f = figure('Position',[200 200 400 150]);
% dat = rand(3);
% cnames = {'X-Data','Y-Data','Z-Data'};
% rnames = {'First','Second','Third'};
% t = uitable('Parent',f,'Data',dat,'ColumnName',cnames,...
%     'RowName',rnames,'Position',[20 20 360 100]);
% 

% filename = 'rician_all_k,tau_3vec.xlsx';
% xlsrange = 'b3:d7';

[num,txt,raw] = xlsread(filename, xlsrange);
raw
sensor_length = size(raw, 1);
rnames = cell(1, sensor_length);
for n = 1 : sensor_length
    rnames{n} = sprintf('sensor%d', n);
end

f = figure('Position',[200 200 450 200], ...
    'menubar','none',...
    'name','fading parameters',...
    'numbertitle','off',...
    'resize','off');
dat = raw;
cnames = {'K factor','path delay(delta ratio)','path loss(dB)'};
t = uitable('Parent',f,'Data',dat,'ColumnName',cnames,...
    'RowName',rnames,'Position',[20 20 370 170]);

[rician_param, rician_param_raw] = get_rician_parameter_from_excel_file(filename, xlsrange)
sensor_length = size(rician_param, 1);

rnames = cell(1, sensor_length);
for n = 1 : sensor_length
    rnames{n} = sprintf('sensor%d', n);
end

f = figure('Position',[400 400 450 200], ...
    'menubar','none',...
    'name','fading parameters',...
    'numbertitle','off',...
    'resize','off');
% dat = rician_param;
dat = rician_param_raw;
cnames = {'K factor','path delay(delta ratio)','path loss(dB)'};
t = uitable('Parent',f,'Data',dat,'ColumnName',cnames, ...
    'RowName',rnames,'Position',[20 20 370 170]);

end

