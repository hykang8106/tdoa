function PLdB = hPositioningPathLoss(d, fc)
% hPositioningPathLoss Urban macro line of sight path loss
%
% PLDB = hPositioningPathLoss(D,FC) is the urban macro line of sight path loss in dB per TR36.814 B.1.2.1,
% where D is the distance in meters and FC is the carrier frequency in Hertz.
%
% [usage]
% PLdB = hPositioningPathLoss(3e3, 2.1e9)

%   Copyright 2011-2013 The MathWorks, Inc.
%

speedOfLight = 299792458.0; % Speed of light in m/s

hBS = 25;           % Height of base station antenna
hUT = 1.5;          % Height of user terminal antenna
hdashBS = hBS - 1.0;  % Effective base station antenna height
hdashUT = hUT - 1.0;  % Effective user terminal antenna height

dBP = 4 * hdashBS * hdashUT * fc / speedOfLight; % Break point distance

% Urban macro line of sight path loss in dB per TR36.814 B.1.2.1
if (d < dBP)
    PLdB = 20.0 * log10(d) + 28.0 + 20.0 * log10(fc / 1e9);
else
    PLdB = 40.0 * log10(d) + 7.8 - 18.0 * log10(hdashBS) - ...
        18.0 * log10(hdashUT) + 2.0 * log10(fc / 1e9);
end

end