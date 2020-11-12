function [] = learn_nested_function

disp('This is the parent function');
nestedfx;

x = 5;
nestfun1;

nestfun2;

    function nestedfx
        disp('This is the nested function');
    end

    function nestfun1
        x = x + 1;
        y = 5;
    end

    function nestfun2
        % #### z is local because parent function dont use it
        z = 3;
        z = z + 12;
        z
    end

x
y = y + 3;
y

nestfun;

    function w = nestfun
        w = x * 10;
        w
    end

q = nestfun;
q = q + 19;
q

end
