% Make plots showing robustness to subsampling of receptor and odorants.

%% Load the data

receptor_data = open('save/receptor_robustness.mat');
odorant_data = open('save/odorant_robustness.mat');

%% Set up colors

% set some colors
color_uncertainty = [0.9 0.9 0.9];
color_mean = [0.737 0.180 0.172];

%% Calculate some summary statistics

% for subsampling receptors, focusing on absolute receptor numbers
receptor_data.corr_low = quantile(receptor_data.corr_values, 0.2, 2);
receptor_data.corr_high = quantile(receptor_data.corr_values, 0.8, 2);
receptor_data.corr_mean = nanmean(receptor_data.corr_values, 2);

% for subsampling odorants, focusing on absolute receptor numbers
odorant_data.corr_low = quantile(odorant_data.corr_values, 0.2, 2);
odorant_data.corr_high = quantile(odorant_data.corr_values, 0.8, 2);
odorant_data.corr_mean = nanmean(odorant_data.corr_values, 2);

% for subsampling receptors, focusing on changes in receptor numbers
receptor_data.corr_diff_low = quantile(receptor_data.corr_diff_values, 0.2, 2);
receptor_data.corr_diff_high = quantile(receptor_data.corr_diff_values, 0.8, 2);
receptor_data.corr_diff_mean = nanmean(receptor_data.corr_diff_values, 2);

% for subsampling odorants, focusing on changes in receptor numbers
odorant_data.corr_diff_low = quantile(odorant_data.corr_diff_values, 0.2, 2);
odorant_data.corr_diff_high = quantile(odorant_data.corr_diff_values, 0.8, 2);
odorant_data.corr_diff_mean = nanmean(odorant_data.corr_diff_values, 2);

%% Make the plots

fig = figure;
fig.Units = 'inches';
fig.Color = [1 1 1];
fig.Position(3:4) = [5.8 2.9];

ax = axes;
ax.OuterPosition = [0 0.5 0.5 0.5];

% draw the uncertainty area
receptor_data.actual_fractions = round(...
    receptor_data.receptor_fractions*size(receptor_data.S, 1)) / ...
        size(receptor_data.S, 1);
fill([flipud(receptor_data.actual_fractions(:)) ; receptor_data.actual_fractions(:)], ...
     [flipud(receptor_data.corr_low(:)) ; receptor_data.corr_high(:)], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(receptor_data.actual_fractions, receptor_data.corr_mean, ...
    'color', color_mean, 'linewidth', 1);

% label the axes
xlabel('fraction of receptors removed');
ylabel('robustness of K_i');

% plot the 0 correlation line as a visual guide
xlim([0 0.8333]);
% plot(xlim, [0 0], 'k:');
ylim([0 1]);

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5);

ax = axes;
ax.OuterPosition = [0.5 0.5 0.5 0.5];

% draw the uncertainty area
odorant_data.actual_fractions = round(...
    odorant_data.odorant_fractions*size(odorant_data.S, 1)) / ...
        size(odorant_data.S, 1);
fill([flipud(odorant_data.actual_fractions(:)) ; odorant_data.actual_fractions(:)], ...
     [flipud(odorant_data.corr_low(:)) ; odorant_data.corr_high(:)], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(odorant_data.actual_fractions, odorant_data.corr_mean, ...
    'color', color_mean, 'linewidth', 1);

% label the axes
xlabel('fraction of odorants removed');
ylabel('robustness of K_i');

% plot the 0 correlation line as a visual guide
xlim([0 0.8333]);
% plot(xlim, [0 0], 'k:');
ylim([0 1]);

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5);

ax = axes;
ax.OuterPosition = [0 0 0.5 0.5];

% draw the uncertainty area
receptor_data.actual_fractions = round(...
    receptor_data.receptor_fractions*size(receptor_data.S, 1)) / ...
        size(receptor_data.S, 1);
fill([flipud(receptor_data.actual_fractions(:)) ; receptor_data.actual_fractions(:)], ...
     [flipud(receptor_data.corr_diff_low(:)) ; receptor_data.corr_diff_high(:)], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(receptor_data.actual_fractions, receptor_data.corr_diff_mean, ...
    'color', color_mean, 'linewidth', 1);

% label the axes
xlabel('fraction of receptors removed');
ylabel('robustness of \DeltaK_i');

% plot the 0 correlation line as a visual guide
xlim([0 0.8333]);
% plot(xlim, [0 0], 'k:');
ylim([0 1]);

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5);

ax = axes;
ax.OuterPosition = [0.5 0 0.5 0.5];

% draw the uncertainty area
odorant_data.actual_fractions = round(...
    odorant_data.odorant_fractions*size(odorant_data.S, 1)) / ...
        size(odorant_data.S, 1);
fill([flipud(odorant_data.actual_fractions(:)) ; odorant_data.actual_fractions(:)], ...
     [flipud(odorant_data.corr_diff_low(:)) ; odorant_data.corr_diff_high(:)], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(odorant_data.actual_fractions, odorant_data.corr_diff_mean, ...
    'color', color_mean, 'linewidth', 1);

% label the axes
xlabel('fraction of odorants removed');
ylabel('robustness of \DeltaK_i');

% plot the 0 correlation line as a visual guide
xlim([0 0.8333]);
% plot(xlim, [0 0], 'k:');
ylim([0 1]);

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5);

preparegraph;

safeprint(fullfile('figs', 'robustness_receptor_and_odorant'));
