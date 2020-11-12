function [fh1, fh2, c] = sub_get_handle(a, b)
% use in learn_function_handle.m
% you can use function handles to call a function that may otherwise be hidden or out of scope

c = a * b;

fh1 = @subfun1;
fh2 = @subfun2;

    function d = subfun1(e, f)
        
        d = e + f;
        
    end

    function g = subfun2(h, i)
        
        g = h - i;
        
    end

end
