function [] = sub_plot_tdoa_result_bw_corr_snr(filename, snr_idx, nsubframe_idx)

load(filename);

sensor_position, target_position, trial_length;

ndlrb, snr_vec, nprsrb_vec, nsubframe_vec

% ###################################################################################
% ### new mat(after 1611011130) have 'rician_param'
% ### old mat(before 1611011130) have 'fading_param_filename', not 'rician_param
% ###################################################################################
if exist('rician_param', 'var')
    rician_param
end

if exist('target_position_limit', 'var')
    target_position_limit
end

if exist('radius_ratio', 'var')
    radius_ratio
end

% ### reminder ###
% ### in batch_simulate_tdoa_bw_corr_snr.m
% ######### tip: for dimension rearrange, use "permute" function
% position_error_array = zeros(nprsrb_len, snr_len, nsubframe_len);

[nprsrb_len, snr_len, nsubframe_len] = size(position_error_array);
if snr_idx > snr_len
    snr_idx = snr_len;
    fprintf(2, '##### snr_idx: set to %d\n', snr_len);
end
    
if nsubframe_idx > nsubframe_len
    nsubframe_idx = nsubframe_len;
    fprintf(2, '##### nsubframe_idx: set to %d\n', nsubframe_len);
end

bw_mhz_vec = zeros(1, nprsrb_len);
bw_mhz_vec_cell = cell(1, nprsrb_len);
for n = 1 : nprsrb_len
    [bw_mhz_vec(n), fs, nfft, sample_length] = get_bw_from_prs_spec_db(ndlrb, nprsrb_vec(n), nsubframe_vec(1));
    bw_mhz_vec_cell{n} = sprintf('%.1f', bw_mhz_vec(n));
end
bw_mhz_vec;
    
corr_length_vec = zeros(1, nsubframe_len);
for n = 1 : nsubframe_len
    [bw_mhz, fs, nfft, corr_length_vec(n)] = get_bw_from_prs_spec_db(ndlrb, nprsrb_vec(1), nsubframe_vec(n));
end
corr_length_vec;

legend_text = cell(1, snr_len);
for n = 1 : snr_len
    legend_text{n} = sprintf('snr %d db', snr_vec(n));
end

if ~exist('rician_param', 'var') || isempty(rician_param)
    for n = 1 : nsubframe_len
        position_error = position_error_array(:, :, n);
        
        figure;
        plot(bw_mhz_vec, position_error, '-s');
        % prevent xtick label from overlapping
        if nprsrb_len <= 15
            set(gca, 'xtick', bw_mhz_vec);
        end
        xlabel('target signal bw in mhz');
        ylabel('meter');
        grid on;
        legend(legend_text);
        title_text = ...
            sprintf('[location error, trial = %d, ndlrb = %d] correlation = %d sample', ...
            trial_length, ndlrb, corr_length_vec(n));
        title(title_text);
    end
end

if exist('error_cell_array', 'var') && exist('rician_param', 'var') && ~isempty(rician_param)
    
    size(error_cell_array)
    
    histo_bin_length = 50;
    
    [max_error] = get_max_in_cell_array(error_cell_array);
    
    x = linspace(0, max_error, histo_bin_length);
%     x = linspace(0, target_position_limit, histo_bin_length);
    
    [H] = make_position_error_histo_array(error_cell_array, x);
    % H dimension = histo_bin_length x nprsrb_len x snr_len x nsubframe_len
    size(H);
    
    B = H(:, :, snr_idx, nsubframe_idx);
    size(B);
    B;
    
    figure('Position',[360 389 744 533], 'name', filename);
    %     width = .5;
    width = 1;
    h = bar3(x, B, width);
    xlabel('bw(mhz)');
    ylabel('meter');
    zlabel('count');
    set(gca, 'XTickLabel', bw_mhz_vec_cell);
    get(gca,'PlotBoxAspectRatio');
    R = get(gca,'DataAspectRatio');
    R(2) = R(2) * 1.2;
    set(gca, 'DataAspectRatio', R);
    title(sprintf('[location error histogram] snr = %d db, correlation = %d samples, trial = %d', ...
        snr_vec(snr_idx), corr_length_vec(nsubframe_idx), trial_length));
    for k = 1:length(h)
        zdata = get(h(k),'ZData');
        set(h(k),'CData',zdata,'FaceColor','interp','FaceLighting','phong');
    end
    axis tight;
%     view(-60, 30);
    view(-70, 70);
    brighten(.7);
    
    %     h = colorbar;
    %     bar_pos = get(h, 'position');
    %     bar_pos(4) = bar_pos(4) / 2;
    %     bar_pos(3) = bar_pos(3) / 3;
    %     bar_pos(1) = bar_pos(1) * 1.1;
    %     bar_pos(2) = .3;
    %     bar_pos;
    %     set(h, 'position', bar_pos);
    %     set(get(h, 'YLabel'), 'String', 'count');
    %
    %     % https://www.mathworks.com/matlabcentral/answers/15491-strange-colorbar-behaviour-must-be-a-bug
    %     % prevent colorbar strange behavior
    %     set(gcf,'Renderer', 'zbuffer');
end

end

%%
function [max_error] = get_max_in_cell_array(error_cell_array)

[nprsrb_len, snr_len, nsubframe_len] = size(error_cell_array);

max_error_array = zeros(nprsrb_len, snr_len, nsubframe_len);

for i = 1 : nprsrb_len
    for j = 1 : snr_len
        for k = 1 : nsubframe_len
            max_error_array(i, j, k) = max(error_cell_array{i, j, k});
        end
    end
end

max_error = max(max(max(max_error_array)));

end

