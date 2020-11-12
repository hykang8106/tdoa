function [bw_mhz, fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb, subframe_length)
% nicer version than old_get_bw_from_prs_spec_db
%
% [usage]
% get_bw_from_prs_spec_db(15, 2, 1)

% subcarrier spacing = 15e3 hz
% prs bw: 1.4 MHz (6RB) | 3 MHz (15RB) | 5 MHz (25RB) | 10 MHz (50RB) | 15 MHz (75RB) | 20 MHz (100RB)
% reference = http://rfmw.em.keysight.com/wireless/helpfiles/n7624b/Content/Main/Positioning%20Reference%20Signals%20(PRS)%20-%20(Advanced%20LTE%20FDD%20Downlink).htm

% bw_mhz = [];

% ###### column(spec) meaning:
% ###### ndlrb, nfft(subcarrier number), bw_mhz
prs_spec_db = [
    6, 128, 1.4;
    15, 256, 3;
    25, 512, 5;
    50, 1024, 10;
    75, 2048, 15;
    100, 2048, 20;
    ];

subcarrier_spacing_hz = 15e3;

idx = find(prs_spec_db(:, 1) == ndlrb);
if isempty(idx)
    % tip: when fprintf fid is 2(std error), message color is red!
    fprintf(2, '###### check number of downlink resource block: %d not in list\n', ndlrb);
    return;
end

nfft = prs_spec_db(idx, 2);
fs = nfft * subcarrier_spacing_hz;

sample_length_per_subframe = fs / 1e3;

bw_mhz = prs_spec_db(idx, 3) * nprsrb / ndlrb;

sample_length = sample_length_per_subframe * subframe_length;
    
% ndlrb nfft  fs(mhz) sample(1 subframe)
% 6	    128   1.92    1920
% 15	256   3.84    3840
% 25	512   7.68    7680
% 50	1024  15.36   15360
% 75	2048  30.72   30720
% 100	2048  30.72   30720

% ### blind touch elephant ###
%
% [ndlrb = 6]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
%
% [ndlrb = 15]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
% 8        1.5      4
% 10       2        4
% 12       2.4      4
% 14       2.6      4
%
% [ndlrb = 25]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
% 8        1.5      4
% 10       2        4
% 12       2.4      4
% 14       2.6      4
% 16       3        4
% 18       3.2      4
% 20       3.4      4
% 22       4        4
% 24       4.2      4

% bw_mhz = [];
% pulse_length = [];
% 
% % ###### column(spec) meaning:
% % ###### ndlrb, nprsrb, bw_mhz, pulse_length
% % ###### bw_mhz is rough value by eye inspection, NOT EXACT
% prs_spec_db = [
%     6, 2, 0.4, 2;
%     6, 4, 0.8, 2;
%     6, 6, 1, 2;
%     15, 2, 0.4, 2;
%     15, 4, 0.8, 2;
%     15, 6, 1, 2;
%     15, 8, 1.5, 2;
%     15, 10, 2, 2;
%     15, 12, 2.4, 2;
%     15, 14, 2.6, 2;
%     25, 2, 0.4, 2;
%     25, 4, 0.8, 2;
%     25, 6, 1, 2;
%     25, 8, 1.5, 4;
%     25, 10, 2, 4;
%     25, 12, 2.4, 4;
%     25, 14, 2.6, 4;
%     25, 16, 3, 4;
%     25, 18, 3.2, 4;
%     25, 20, 3.4, 4;
%     25, 22, 4, 4;
%     25, 24, 4.2, 4;
%     ];
% 
% db_length = size(prs_spec_db, 1);
% for n = 1 : db_length
%     if (prs_spec_db(n, 1) == ndlrb) && (prs_spec_db(n, 2) == nprsrb)
%         bw_mhz = prs_spec_db(n, 3);
%         pulse_length = prs_spec_db(n, 4);
%         break;
%     else
%         continue;
%     end
% end

end

%%
function [bw_mhz, pulse_length] = old_get_bw_from_prs_spec_db(ndlrb, nprsrb)

% subcarrier spacing = 15e3 hz
% Choice: 1.4 MHz (6RB) | 3 MHz (15RB) | 5 MHz (25RB) | 10 MHz (50RB) | 15 MHz (75RB) | 20 MHz (100RB)
% http://rfmw.em.keysight.com/wireless/helpfiles/n7624b/Content/Main/...
% Positioning%20Reference%20Signals%20(PRS)%20-%20(Advanced%20LTE%20FDD%20Downlink).htm

% ndlrb nfft  fs(mhz) sample
% 6	    128   1.92    1920
% 15	256   3.84    3840
% 25	512   7.68    7680
% 50	1024
% 75	2048
% 100	2048

% ### blind touch elephant ###
%
% [ndlrb = 6]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
%
% [ndlrb = 15]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
% 8        1.5      4
% 10       2        4
% 12       2.4      4
% 14       2.6      4
%
% [ndlrb = 25]
% nprsrb   bw_mhz   pulse_length
% 2        0.4      2
% 4        0.8      2
% 6        1        2
% 8        1.5      4
% 10       2        4
% 12       2.4      4
% 14       2.6      4
% 16       3        4
% 18       3.2      4
% 20       3.4      4
% 22       4        4
% 24       4.2      4

bw_mhz = [];
pulse_length = [];

% ###### column(spec) meaning:
% ###### ndlrb, nprsrb, bw_mhz, pulse_length
% ###### bw_mhz is rough value by eye inspection, NOT EXACT
prs_spec_db = [
    6, 2, 0.4, 2;
    6, 4, 0.8, 2;
    6, 6, 1, 2;
    15, 2, 0.4, 2;
    15, 4, 0.8, 2;
    15, 6, 1, 2;
    15, 8, 1.5, 2;
    15, 10, 2, 2;
    15, 12, 2.4, 2;
    15, 14, 2.6, 2;
    25, 2, 0.4, 2;
    25, 4, 0.8, 2;
    25, 6, 1, 2;
    25, 8, 1.5, 4;
    25, 10, 2, 4;
    25, 12, 2.4, 4;
    25, 14, 2.6, 4;
    25, 16, 3, 4;
    25, 18, 3.2, 4;
    25, 20, 3.4, 4;
    25, 22, 4, 4;
    25, 24, 4.2, 4;
    ];

db_length = size(prs_spec_db, 1);
for n = 1 : db_length
    if (prs_spec_db(n, 1) == ndlrb) && (prs_spec_db(n, 2) == nprsrb)
        bw_mhz = prs_spec_db(n, 3);
        pulse_length = prs_spec_db(n, 4);
        break;
    else
        continue;
    end
end

end

