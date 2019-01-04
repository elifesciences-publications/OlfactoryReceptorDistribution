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

% emphasize regions that have the same sign, and regions with opposite sign
max1 = 1.10*max([diff1(:) ; diff2(:)]);
min1 = 1.10*min([diff1(:) ; diff2(:)]);

max2 = max1;
min2 = min1;

hold on;
fill([0, max1, max1, 0, 0], [0, 0, max2, max2, 0], hex2color('EBF2F8'), ...
    'edgecolor', 'none');
fill([0, min1, min1, 0, 0], [0, 0, min2, min2, 0], hex2color('EBF2F8'), ...
    'edgecolor', 'none');
fill([0, max1, max1, 0, 0], [0, 0, min2, min2, 0], hex2color('FAF0F0'), ...
    'edgecolor', 'none');
fill([0, min1, min1, 0, 0], [0, 0, max2, max2, 0], hex2color('FAF0F0'), ...
    'edgecolor', 'none');

smartscatter(diff1(:), diff2(:), 'size', 10, 'alpha', 1);
% smartscatter(diff1(:), diff2(:), 'size', 5, 'color', hex2color('BC2E2C'), ...
%     'alpha', 0.2);
xlabel('\DeltaK_1');
ylabel('\DeltaK_2');

% data_lim = 1.1*max(abs([diff1(:) ; diff2(:)]));
% xlim([-data_lim, data_lim]);
% ylim([-data_lim, data_lim]);

axis equal;

xlim([min1 max1]);

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

% emphasize regions that have the same sign, and regions with opposite sign
sub_max1 = 1.10*max([sub_diff1(:) ; sub_diff2(:)]);
sub_min1 = 1.10*min([sub_diff1(:) ; sub_diff2(:)]);

sub_max2 = sub_max1;
sub_min2 = sub_min1;

hold on;
fill([0, sub_max1, sub_max1, 0, 0], ...
     [0, 0, sub_max2, sub_max2, 0], hex2color('EBF2F8'), ...
     'edgecolor', 'none');
fill([0, sub_min1, sub_min1, 0, 0], ...
     [0, 0, sub_min2, sub_min2, 0], hex2color('EBF2F8'), ...
     'edgecolor', 'none');
fill([0, sub_max1, sub_max1, 0, 0], ...
     [0, 0, sub_min2, sub_min2, 0], hex2color('FAF0F0'), ...
     'edgecolor', 'none');
fill([0, sub_min1, sub_min1, 0, 0], ...
     [0, 0, sub_max2, sub_max2, 0], hex2color('FAF0F0'), ...
     'edgecolor', 'none');

smartscatter(sub_diff1, sub_diff2, 'size', 10, 'alpha', 1);
xlabel('\DeltaK_1');
ylabel('\DeltaK_2');

axis equal;

xlim([sub_min1 sub_max1]);

% xlim([-0.0019 0.0019]);
% ylim([-0.0019 0.0019]);

colormap(cmap);
% colorbar;

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);
preparegraph;

safeprint(fullfile('figs', 'context_dependence_zoom'));
