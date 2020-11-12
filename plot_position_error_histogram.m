function [] = plot_position_error_histogram(position_error_excl, bin_length, title_text, ...
    fading_param_filename)

error_max = max(position_error_excl);

x = linspace(0, error_max, bin_length);

N = histc(position_error_excl, x);

if ~isempty(fading_param_filename)
    figure('name', fading_param_filename);
else
    figure;
end
bar(x, N);
grid on;
xlabel('position error in meter');
ylabel('count');

title(title_text);

end

