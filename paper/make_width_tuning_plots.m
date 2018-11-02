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
set(ax, 'xtick', tuning_min_max, 'xticklabel', {'narrow', 'wide'}, ...
    'xscale', 'log');

% label the axes
xlabel('receptor tuning');
ylabel('corr(log Q_{aa}, K_a)');

% plot the 0 correlation line as a visual guide
plot(xlim, [0 0], 'k:');

% beautify, making sure fonts aren't too big
beautifygraph('fontscale', 0.667);

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
set(ax, 'xtick', tuning_min_max, 'xticklabel', {'narrow', 'wide'}, ...
    'xscale', 'log');

% label the axes
xlabel('receptor tuning');
ylabel('corr(-Q^{-1}_{aa}, K_a)');

% plot the 0 correlation line as a visual guide
plot(xlim, [0 0], 'k:');

% beautify, making sure fonts aren't too big
beautifygraph('fontscale', 0.667);

% adjust figure for printing
preparegraph;

safeprint(fullfile('figs', 'logic_vs_tuning.pdf'));

%% Show how well diag(inv(Q)) predicts K at high SNR

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2.4 1.4];

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
scatterfit(crt_diag_invQ, crt_K, ...
    'scatteropts', {'color', [0.176, 0.459, 0.733], 'size', 60, 'filled', true}, ...
    'fitopts', {'line', [-1, approx_const], 'style', {'k--', 'linewidth', 1}, ...
    'legend', false});

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
beautifygraph('fontscale', 0.667);

% adjust figure for printing
preparegraph;

%% How well does diag(inv(Q)) predict K at wide tuning

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2.4 1.4];

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
smartscatter(crt_diag_invQ, crt_K, 'color', [0.176, 0.459, 0.733], ...
    'size', 60, 'filled', true, 'alpha', 0.3);

ylim([0, max(crt_K)*1.2]);

% label the axes
xlabel('(Q^{-1})_{aa}');
ylabel('K_a');

% beautify, making sure fonts aren't too big
beautifygraph('fontscale', 0.667);

% adjust figure for printing
preparegraph;

% safe_print(fullfile('figs', 'K_vs_invQ.pdf'));