function [H] = make_position_error_histo_array(error_cell_array, x)

[nprsrb_len, snr_len, nsubframe_len] = size(error_cell_array);

histo_bin_length = length(x);
H = zeros(histo_bin_length, nprsrb_len, snr_len, nsubframe_len);

for i = 1 : nprsrb_len
    for j = 1 : snr_len
        for k = 1 : nsubframe_len
            H(:, i, j, k) = histc(error_cell_array{i, j, k}, x);
        end
    end
end
size(H);

end

