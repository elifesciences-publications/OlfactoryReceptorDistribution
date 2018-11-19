% Make plot of Ibarra-Soria-like simulations.

%% Load the data

load(fullfile('save', 'ibarra_soria_like.mat'));

%% Make the plot

K1_mean = mean(K1, 2);
K1_std = std(K1, [], 2);

K2_mean = mean(K2, 2);
K2_std = std(K2, [], 2);

logK_ratio = log2(K2_mean ./ K1_mean);
logK_std = (K1_std ./ K1_mean) + (K2_std ./ K2_mean);

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2.9 2.3];

errorbar(K1_mean, logK_ratio, logK_std, logK_std, K1_std, K1_std, ...
    'marker', 'none', 'linestyle', 'none', 'color', [0.7 0.7 0.7], ...
    'capsize', 0);
hold on;
plot(K1_mean, logK_ratio, '.', 'markersize', 10);

% ylim([-3.5, 3.5]);
ylim([-2.5, 2.5]);
set(gca, 'xscale', 'log');
% xlim([5, 350]);
% xlim([1 max(K1p_mean)*1.2]);
xlim([1 100]);
h = plot(xlim, [0 0], 'r');
set(h, 'color', [1 0 0 0.3], 'linewidth', 1);

xlabel('K^{env1}');
ylabel('log_2 K^{env2}/K^{env1}');

beautifygraph('fontscale', 0.667, 'ticksize', 18);
% beautifygraph('fontscale', 0.667);

set(gca, 'box', 'on', 'xminortick', 'off', 'yminortick', 'off', 'TickLength', [0 0]);

preparegraph('edge', 0);

safeprint(fullfile('figs', 'example_ibarra_soria_style'));
