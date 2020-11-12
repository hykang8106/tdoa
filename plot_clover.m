function [] = plot_clover(leaf_number)

t = (0 : .001: 2 * pi)';
x = cos(leaf_number * t) .* cos(t);
y = cos(leaf_number * t) .* sin(t);

figure;
plot(x, y);
axis equal;
axis tight;
grid on;

theta = 0 : .001 : 2 * pi;
rho = cos(leaf_number * theta);

figure;
polar(theta, rho);

end
