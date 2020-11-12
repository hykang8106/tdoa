function [H] = make_torrieri_H_matrix(len, ref)
% make torrieri H matrix
%
% [input]
% - len: sensor length
% - ref: reference sensor number
% [usage]
% make_torrieri_H_matrix(4, 1)
%

H = zeros(len - 1, len);

i = 1 : len - 1;
j = setdiff(1 : len, ref);

for n = 1 : len - 1
    H(i(n), j(n)) = -1;
end

H(:, ref) = 1;
H;

end