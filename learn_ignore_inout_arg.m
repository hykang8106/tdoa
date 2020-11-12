function [] = learn_ignore_inout_arg()

a = 1;
b = 2;
c = 3;

[d, ~, f] = sub_ignore(a, b);

end

%%
function [d, e, f] = sub_ignore(a, b, c)

d = a + b

e = a - b;

f = a * b;

% d = a + b + c;
% 
% e = a - b - c;
% 
% f = a * b * c;

end
