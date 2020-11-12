%hPositioningTDOA Time difference of arrival
%   Y = hPositioningTDOA(X,SR) is a TDOA matrix given vector of sample
%   arrival times X and sampling rate SR.

%   Copyright 2011-2013 The MathWorks, Inc.

function y = hPositioningTDOA(x,sr)

    % Compute the number of samples delay between arrivals
    y = zeros(length(x));
    for j = 1:length(x)
        for i = (j+1):length(x)
            y(i,j) = x(i)-x(j);
        end
    end
    
    y
    % Convert to time
    y = y./sr;

end