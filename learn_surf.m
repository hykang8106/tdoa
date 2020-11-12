function [] = learn_surf(m)

M = magic(m);
M = repmat(M, 1, 2);
M
size(M)

figure;
surf(M');
axis xy; 
axis tight; 
colormap(jet); view(0, 90);

% plot colorbar
h = colorbar;
set(get(h, 'YLabel'), 'String', 'meter');

end

