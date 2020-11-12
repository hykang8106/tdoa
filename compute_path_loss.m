function [path_loss_db] = compute_path_loss(freq_mhz, distance_km, mobile_antenna_height, base_antenna_height)
% compute path loss only for metropolitan area
%
% [input]
% - freq_mhz: 150 ~ 2000
% - distance_km: 1 ~ 10
% - mobile_antenna_height: 1 ~ 10
% - base_antenna_height: 30 ~ 200
%
% [usage]
% compute_path_loss(800, 5, 10, 50)
%

% check freq
if freq_mhz < 150 || freq_mhz > 2000
    fprintf('error: freq = 150 ~ 2000 mhz\n');
    return;
end

% check distance
if distance_km < 1 || distance_km > 10
    fprintf('error: distance = 1 ~ 10 km\n');
    return;
end

% for hata urban model 
% when 0, small or medium size city
large_city = 1;

% for cost hata model
% C = 3, metropolitan area. C = 0, for suburban or rural area.
C = 3;

if freq_mhz >= 150 && freq_mhz <= 1500
    fprintf('using hata urban model\n');
    [path_loss_db] = hata_urban(freq_mhz, distance_km, mobile_antenna_height, base_antenna_height, large_city);
else
    fprintf('using cost hata model\n');
    [path_loss_db] = cost_hata(freq_mhz, distance_km, mobile_antenna_height, base_antenna_height, C);
end

fprintf('metropolitan area path loss = %f db\n', path_loss_db);
fprintf('free space path loss = %f db\n', fspl(distance_km * 10^3, physconst('LightSpeed')/(freq_mhz * 10^6)));

end

%%
function [path_loss_db] = cost_hata(freq_mhz, distance_km, mobile_antenna_height, base_antenna_height, C)
% [cost hata model]
% https://en.wikipedia.org/wiki/COST_Hata_model
% extends the urban Hata model to cover a more elaborated range of frequencies. 
% It is the most often cited of the COST 231 models
% http://www.lx.it.pt/cost231/final_report.htm
%
% Coverage
% Frequency: 1500?2000 MHz
% Mobile station antenna height: 1?10 m
% Base station Antenna height: 30?200 m
% Link distance: 1?20 km
%
% Limitations
% This model requires that the base station antenna is higher than all adjacent rooftops.

antenna_height_correction_factor = (1.1 * log10(freq_mhz) - 0.7) * mobile_antenna_height ...
    - (1.56 * log10(freq_mhz) - 0.8);

path_loss_db = 46.3 + 33.9 * log10(freq_mhz) - 13.82 * log10(base_antenna_height) ...
    - antenna_height_correction_factor + (44.9 - 6.55 * log10(base_antenna_height)) * log10(distance_km) + C;

end

%%
function [path_loss_db] = hata_urban(freq_mhz, distance_km, mobile_antenna_height, base_antenna_height, large_city)
% [hata model for urban area = okumura-hata model]
% https://en.wikipedia.org/wiki/Hata_model_for_urban_areas
%
% Coverage
% Frequency: 150?1500 MHz
% Mobile Station Antenna Height: 1?10 m
% Base station Antenna Height: 30?200 m
% Link distance: 1?10 km.

if large_city
    if freq_mhz >= 200 && freq_mhz <= 1500 
        % freq = 200 ~ 1500 mhz
        antenna_height_correction_factor = 3.2 * log10(11.75 * mobile_antenna_height)^2 - 4.97;        
    else
        % freq = 150 ~ 200 mhz
        antenna_height_correction_factor = 8.29 * log10(1.54 * mobile_antenna_height)^2 - 1.1;
    end
else
    antenna_height_correction_factor = 0.8 + (1.1 * log10(freq_mhz) - 0.7) * mobile_antenna_height ...
        - 1.56 * log10(freq_mhz);
end

path_loss_db = 69.55 + 26.16 * log10(freq_mhz) - 13.82 * log10(base_antenna_height) ...
    - antenna_height_correction_factor + (44.9 - 6.55 * log10(base_antenna_height)) * log10(distance_km);

end
