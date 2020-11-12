%hPositioningPosition eNodeB position
%   P = hPositioningPosition(I,N) is a random station position P = [px py]
%   where the coordinates px and py are chosen such that for stations I =
%   0...N-1 the stations will be well spread out.

%   Copyright 2011-2013 The MathWorks, Inc.

function p = hPositioningPosition(i,n)
    
    % Position eNodeB randomly within a I*2*pi/N radian sector
    phi = i*2*pi/n+rand(1, 1)*2*pi/(2*n)-2*pi/(2*n);
    
    % Position eNodeB randomly between 4000 + (I*5000/N) and 5000 +
    % (I*5000/N) from the UE
    r = randi([0, 1000], 1, 1)+4000+(i*5000/n);
    
    % Convert coordinates
    [x, y] = pol2cart(phi, r);
    p = [x, y];
    
end
