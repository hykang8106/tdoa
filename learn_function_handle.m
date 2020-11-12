function [] = learn_function_handle
% Function handles can give you access to functions you might not be able to execute. 
% For instance, with a handle you can call a function even if it is no longer on your MATLAB¢ç path. 
% You can also call a local function from outside of the file that defines that function.

p = makeParabola(1.3, .2, 30);

x = 25;
y = p(x)

% matlab command, fplot: plot function
fplot(p, [-25, 25]);

% #### sub_my_func is outside of this file(see sub_my_func.m file), also working!!
d = 10; e = 30;
[f, g] = sub_my_func(d, e);
f, g

% fh: function handle
fh = @sub_my_func;
[h, i] = fh(d, e);
h, i

% If the function being called takes no input arguments, 
% then you must call the function with empty parentheses placed after the handle name. 
% If you use only the handle name, MATLAB just identifies the name of the function.
fh

% anonymous function with multiple output
F = @(x)find(x);

m = [3 2 0; -5 0 7; 0 0 1]

% find nonzero element in m
[row, col, val] = F(m);
row, col, val

a = 100; b = 20;
% ### sub_get_handle is outside of this file(see sub_get_handle.m file)
% you can use function handles to call a function that may otherwise be hidden or out of scope
[f1, f2, c] = sub_get_handle(a, b);

e = 10; f = 300;
d = f1(e, f);
d

h = -2; i = -20;
g = f2(h, i);
g

c

end

%%
function p = makeParabola(a,b,c)

p = @parabola;

    function y = parabola(x)
        y = a*x.^2 + b*x + c;
    end

end

%% #### comment out: move to outside of this file

% function [f, g] = sub_my_func(d, e)
% 
% f = d + e;
% g = d - e;
% 
% end

