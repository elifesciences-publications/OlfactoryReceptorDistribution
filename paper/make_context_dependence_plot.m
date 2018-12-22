% Make plot showing context dependence.

%% Load the data

load(fullfile('save', 'context_dependence.mat'));

%% Make the plot

cmap = linspace(0, 1, 256)' .* hex2color('BC2E2C');

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position(3:4) = [2.4 2];

diff1 = results.K1_pert - results.K1;
diff2 = results.K2_pert - results.K2;
smartscatter(diff1(:), diff2(:), 'size', 10, 'alpha', 1);
% smartscatter(diff1(:), diff2(:), 'size', 5, 'color', hex2color('BC2E2C'), ...
%     'alpha', 0.2);
xlabel('\DeltaK_1');
ylabel('\DeltaK_2');

% data_lim = 1.1*max(abs([diff1(:) ; diff2(:)]));
% xlim([-data_lim, data_lim]);
% ylim([-data_lim, data_lim]);

axis equal;

colormap(cmap);
% colorbar;

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);
preparegraph;

safeprint(fullfile('figs', 'context_dependence'));

%%

cmap = linspace(0, 1, 256)' .* hex2color('BC2E2C');

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position(3:4) = [2.4 2];

diff1 = results.K1_pert - results.K1;
diff2 = results.K2_pert - results.K2;

klim = quantile([abs(diff1(:)) ; abs(diff2(:))], 0.95);

% choose only those diffs that are within the high-density center
sub_mask = (abs(diff1(:)) < klim) & (abs(diff2(:)) < klim);
sub_diff1 = diff1(sub_mask);
sub_diff2 = diff2(sub_mask);

smartscatter(sub_diff1, sub_diff2, 'size', 10, 'alpha', 1);
xlabel('\DeltaK_1');
ylabel('\DeltaK_2');

axis equal;

% xlim([-0.0019 0.0019]);
% ylim([-0.0019 0.0019]);

colormap(cmap);
% colorbar;

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);
preparegraph;

safeprint(fullfile('figs', 'context_dependence_zoom'));
