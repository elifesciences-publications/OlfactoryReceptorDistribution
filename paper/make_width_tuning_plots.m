% Make plots to show how the overlap and inverse overlap matrices relate to
% the optimal receptor distibution at various values of the tuning width in
% an artificial array.

%% Load the tuning results

load(fullfile('save', 'tuning_width_sweep.mat'));

%% Make the plots

% choose the same axis limits for both plots
tuning_min_max = [min(tuning_values) max(tuning_values)];
fixed_y_lim = [-0.1 1];

% set some colors
color_uncertainty = [0.9 0.9 0.9];
color_mean = [0.737 0.180 0.172];

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 3.5 1.3];

% make the axes for the first plot -- K vs. log diag(Q)
ax = axes;
ax.OuterPosition = [0 0 0.5 1];

% draw the uncertainty area
fill([flipud(tuning_values(:)) ; tuning_values(:)], ...
     [flipud(results.summaries.log_diagQ.cp_low(:, :, 1)') ; ...
             results.summaries.log_diagQ.cp_high(:, :, 1)'], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(tuning_values, results.summaries.log_diagQ.cp_mean(:, :, 1), ...
    'color', color_mean, 'linewidth', 1);

% adjust the axes limits
xlim(tuning_min_max);
ylim(fixed_y_lim);

% switch to log space on x axis and adjust the ticks
set(ax, 'xtick', tuning_min_max, ...
    'xticklabel', {[num2str(tuning_min_max(1)) ' (narrow)'], ...
                   [num2str(tuning_min_max(2)) ' (wide)']}, ...
    'xscale', 'log');

% don't waste ink on the axes
ax.LineWidth = 0.5;

% label the axes
xlabel('receptor tuning width');
ylabel('corr(log Q_{aa}/\sigma_a^2, K_a)');

% plot the 0 correlation line as a visual guide
plot(xlim, [0 0], 'k:');

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5);

% make the ticks a bit bigger
ax.TickLength = 2*ax.TickLength;

% make the axes for the first plot -- K vs. -diag(inv(Q))
ax = axes;
ax.OuterPosition = [0.5 0 0.5 1];

% draw the uncertainty area
fill([flipud(tuning_values(:)) ; tuning_values(:)], ...
    -[flipud(results.summaries.diag_invQ.cp_low(:, :, 1)') ; ...
             results.summaries.diag_invQ.cp_high(:, :, 1)'], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(tuning_values, -results.summaries.diag_invQ.cp_mean(:, :, 1), ...
    'color', color_mean, 'linewidth', 1);

% adjust the axes limits
xlim(tuning_min_max);
ylim(fixed_y_lim);

% switch to log space on x axis and adjust the ticks
set(ax, 'xtick', tuning_min_max, ...
    'xticklabel', {[num2str(tuning_min_max(1)) ' (narrow)'], ...
                   [num2str(tuning_min_max(2)) ' (wide)']}, ...
    'xscale', 'log');

% label the axes
xlabel('receptor tuning width');
ylabel('corr(-\sigma_a^2 Q^{-1}_{aa}, K_a)');

% plot the 0 correlation line as a visual guide
plot(xlim, [0 0], 'k:');

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5);

% make the ticks a bit bigger
ax.TickLength = 2*ax.TickLength;

% adjust figure for printing
preparegraph;

safeprint(fullfile('figs', 'logic_vs_tuning.pdf'));

%% Show how well diag(inv(Q)) predicts K at high SNR

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2.4 1.4];

fig.Color = [1 1 1];

crt_Ktot = 15*Ktot;
[crt_K, ~, crt_Q] = calculate_optimal_dist(results.inputs.S{1, end, 1}, ...
    Gamma1, crt_Ktot);
% calculate the diagonal of the inverse overlap matrix
crt_invQ = inv(crt_Q);
crt_diag_invQ = diag(crt_invQ);

% calculate the constant term from the large-SNR approximation
approx_const = crt_Ktot/n_receptors + mean(crt_diag_invQ);
% show the dependence between the diagonal of the inverse overlap matrix
% and the optimal receptor counts
hold on;
drawfitline(crt_diag_invQ, crt_K, 'line', [-1, approx_const], ...
    'style', {'k--', 'linewidth', 1}, 'legend', false);
h = smartscatter(crt_diag_invQ, crt_K, 'color', [0.176, 0.459, 0.733], ...
    'size', 60, 'filled', true);
% [~, ~, h] = scatterfit(crt_diag_invQ, crt_K, ...
%     'scatteropts', {'color', [0.176, 0.459, 0.733], 'size', 60, 'filled', true}, ...
%     'fitopts', {'line', [-1, approx_const], 'style', {'k--', 'linewidth', 1}, ...
%     'legend', false});
% add a thin edge around the data points
h.hscatter.MarkerEdgeColor = [1 1 1];

% ylim([0, max(crt_K)*1.2]);
hold on;
plot([0 max(ylim)], [0 max(ylim)], ':k');

% axis equal;

xlim([0 max(ylim)]);
ylim([0 max(ylim)]);

% label the axes
xlabel('(Q^{-1})_{aa}');
ylabel('K_a');

% beautify, making sure fonts aren't too big
beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'minorticks', 'off', ...
    'ticksize', 10);

% adjust figure for printing
preparegraph;

safeprint(fullfile('figs', 'K_vs_invQ'));

%% How well does diag(inv(Q)) predict K at wide tuning

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2.4 1.4];

fig.Color = [1 1 1];

% choose a particular result
idx_rep = 1;
idx_tuning = length(tuning_values);
crt_Q = results.outputs.Q{idx_rep, idx_tuning, 1};
crt_K = results.outputs.K{idx_rep, idx_tuning, 1};

% calculate the diagonal of the inverse overlap matrix
crt_invQ = inv(crt_Q);
crt_diag_invQ = diag(crt_invQ);

% calculate the constant term from the large-SNR approximation
approx_const = Ktot/n_receptors + mean(crt_diag_invQ);
% show the dependence between the diagonal of the inverse overlap matrix
% and the optimal receptor counts
% scatterfit(crt_diag_invQ, crt_K, ...
%     'scatteropts', {'color', [0.176, 0.459, 0.733], 'size', 60, 'filled', true}, ...
%     'fitopts', {'line', [-1, approx_const], 'style', {'k--', 'linewidth', 1}, ...
%     'legend', false});
h = smartscatter(crt_diag_invQ, crt_K, 'color', [0.176, 0.459, 0.733], ...
    'size', 60, 'filled', true, 'alpha', 0.3);
% add a thin edge around the data points
h.hscatter.MarkerEdgeColor = [1 1 1];

ylim([0, max(crt_K)*1.2]);

% label the axes
xlabel('(Q^{-1})_{aa}');
ylabel('K_a');

% beautify, making sure fonts aren't too big
beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'ticksize', 10);

% adjust figure for printing
preparegraph;

%% Make plots showing how well we can predict \Delta K

% choose the same axis limits for both plots
tuning_min_max = [min(tuning_values) max(tuning_values)];
fixed_y_lim = [-0.1 1];

% set some colors
color_uncertainty = [0.9 0.9 0.9];
color_mean = [0.737 0.180 0.172];

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 3.5 1.3];

% make the axes for the first plot -- K vs. log diag(Q)
ax = axes;
ax.OuterPosition = [0 0 0.5 1];

% draw the uncertainty area
fill([flipud(tuning_values(:)) ; tuning_values(:)], ...
     [flipud(results.summaries.diff_log_diagQ.cp_low(:, :, 1)') ; ...
             results.summaries.diff_log_diagQ.cp_high(:, :, 1)'], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(tuning_values, results.summaries.diff_log_diagQ.cp_mean(:, :, 1), ...
    'color', color_mean, 'linewidth', 1);

% adjust the axes limits
xlim(tuning_min_max);
ylim(fixed_y_lim);

% switch to log space on x axis and adjust the ticks
set(ax, 'xtick', tuning_min_max, ...
    'xticklabel', {[num2str(tuning_min_max(1)) ' (narrow)'], ...
                   [num2str(tuning_min_max(2)) ' (wide)']}, ...
    'xscale', 'log');

% don't waste ink on the axes
ax.LineWidth = 0.5;

% label the axes
xlabel('receptor tuning width');
ylabel('corr(\Deltalog Q_{aa}/\sigma_a^2, \DeltaK_a)');

% plot the 0 correlation line as a visual guide
plot(xlim, [0 0], 'k:');

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5, 'labelsize', 10);

% make the ticks a bit bigger
ax.TickLength = 2*ax.TickLength;

% make the axes for the first plot -- K vs. -diag(inv(Q))
ax = axes;
ax.OuterPosition = [0.5 0 0.5 1];

% draw the uncertainty area
fill([flipud(tuning_values(:)) ; tuning_values(:)], ...
    -[flipud(results.summaries.diff_log_diag_invQ.cp_low(:, :, 1)') ; ...
             results.summaries.diff_log_diag_invQ.cp_high(:, :, 1)'], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(tuning_values, -results.summaries.diff_log_diag_invQ.cp_mean(:, :, 1), ...
    'color', color_mean, 'linewidth', 1);

% adjust the axes limits
xlim(tuning_min_max);
ylim(fixed_y_lim);

% switch to log space on x axis and adjust the ticks
set(ax, 'xtick', tuning_min_max, ...
    'xticklabel', {[num2str(tuning_min_max(1)) ' (narrow)'], ...
                   [num2str(tuning_min_max(2)) ' (wide)']}, ...
    'xscale', 'log');

% label the axes
xlabel('receptor tuning width');
ylabel('corr(-\Deltalog \sigma_a^2Q^{-1}_{aa}, \DeltaK_a)');

% plot the 0 correlation line as a visual guide
plot(xlim, [0 0], 'k:');

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5, 'labelsize', 10);

% make the ticks a bit bigger
ax.TickLength = 2*ax.TickLength;

% adjust figure for printing
preparegraph;

safeprint(fullfile('figs', 'diff_logic_vs_tuning.pdf'));

%% What happens with the SNR as we increase tuning width?

% calculate SNR
% all_snr = cellfun(@(Q) sqrt(Ktot*trace(Q)), results.outputs.Q);
all_snr = cellfun(@(Q) sqrt(trace(Q)/size(Q, 1)), results.outputs.Q);
mean_snr = mean(all_snr, 1);
lo_snr = quantile(all_snr, 0.2, 1);
hi_snr = quantile(all_snr, 0.8, 1);

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [3 2];

fig.Color = [1 1 1];

% set some colors
color_uncertainty = [0.9 0.9 0.9];
color_mean = [0.737 0.180 0.172];

hold on;
% draw the uncertainty area
fill([flipud(tuning_values(:)) ; tuning_values(:)], ...
     [flipud(lo_snr(:, :, 1)') ; hi_snr(:, :, 1)'], ...
     color_uncertainty, 'linestyle', 'none');
% draw the mean
plot(tuning_values, mean_snr(:, :, 1), 'color', color_mean, 'linewidth', 1);

all_n_expressed = cellfun(@(K) sum(K > 1e-3), results.outputs.K);
mean_n_expressed = mean(all_n_expressed, 1);
lo_n_expressed = quantile(all_n_expressed, 0.2, 1);
hi_n_expressed = quantile(all_n_expressed, 0.8, 1);

xlabel('Tuning width');
ylabel('SNR');

beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'ticksize', 10);
preparegraph;

figure;
scatter(all_snr(:), all_n_expressed(:)/n_receptors, 'filled', 'markeredgealpha', 0.2, ...
    'markerfacealpha', 0.2);
% plot(mean_snr(:, :, 1), mean_n_expressed(:, :, 1), 'color', color_mean, 'linewidth', 1);
xlabel('SNR');
ylabel('Fraction of expressed receptors');

beautifygraph;
preparegraph;
