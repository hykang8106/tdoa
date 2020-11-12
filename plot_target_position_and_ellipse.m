function [] = plot_target_position_and_ellipse(target_position, ...
    x_est_torrieri, cep, gdop, lambda1, lambda2, theta_rad, legendstr)

trial_length = size(x_est_torrieri, 2);

plot(target_position(1), target_position(2), 'ko', 'MarkerSize', 7, 'MarkerFaceColor', 'r', ...
    'LineWidth', 2);
legendstr = [legendstr, 'target'];
legend(legendstr, 'AutoUpdate', 'off');

for n = 1 : trial_length
    hold on;
    plot(x_est_torrieri(1, n), x_est_torrieri(2, n), 'b+', 'LineWidth', 2);
end

% concentration ellipse corresponding to probability Pe: 
% major axes length = 2 * sqrt(k * lambda1) 
% minor axes length = 2 * sqrt(k * lambda2)
% where k = -2 * log(1 - Pe)
% see eq (60)

Pe = .5;
k = -2 * log(1 - Pe);
ra = 2 * sqrt(k * lambda1);
rb = 2 * sqrt(k * lambda2);
my_ellipse(ra, rb, theta_rad, x_est_torrieri(1, :), x_est_torrieri(2, :), 'r');

end

