%hPositioningTDOACurve TDOA curve
%   [X,Y] = hPositioningTDOACurve(P1,P2,DD) computes a constant TDOA
%   hyperbolic curve as a series of [X,Y] locations given station positions
%   P1 and P2 and the distance difference DD.

%   Copyright 2011-2013 The MathWorks, Inc.

function [x, y] = hPositioningTDOACurve(p1, p2, dd)

% Calculate vector between stations
delta = p1 - p2;

% Express in polar form
[phi, r] = cart2pol(delta(1), delta(2));

% Split the vector between stations such that the distance differs by
% dd
rd = (r + dd) / 2;

% Radius rd gives distance between hyperbola focus and transverse
% axis crossing. Radius r give the distance between foci. Use this
% to compute hyperbola parameters a,b,c and e.
a = (r / 2 )- rd;
e = r / (2 * a); % eccentricity
c = a * e;
b = sqrt(c^2 - a^2);

% Use parametric equation for the hyperbola in terms of a, b and
% parameter mu. [x,y] is initially calculated for a hyperbola whose
% transverse axis is the x-axis and whose foci are +/- c. The foci
% is returned to the origin by subtracting c, the hyperbola is rotated
% such that the transverse axis lies along vector delta joining the
% stations, and finally the hyperbola is translated such that the
% foci lies on the position of station 1.

% ##### 499 is constant, consider uca_radius_meter and radius_ratio
% ##### see pah.m which i wrote
x = zeros(1, 499);
y = zeros(1, 499);
for k = 1 : 499
    mu = (k - 249) / 50;
    x(k) = a * cosh(mu);
    y(k) = b * sinh(mu);
    x(k) = x(k) - c;
    y(k) = y(k);
    [phi2, r2] = cart2pol(real(x(k)), real(y(k)));
    [x(k), y(k)] = pol2cart(phi2 + phi, r2);
end

x = x + p1(1);
y = y + p1(2);

end


