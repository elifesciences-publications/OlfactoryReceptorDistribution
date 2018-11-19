% Make plots analyzing the dynamical model.

%% Load the results

load(fullfile('save', 'dynmodel_results.mat'));

%% Make the plots

% plot the dynamical evolution of receptor abundances
yl = max([max(dyn_history_from_random.K(:)) max(dyn_history_from_env1.K(:))]);
    
col1 = hex2color('2D75BB');
col2 = hex2color('EE8434');

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.8 2];

% subsample the time points we plot, to avoid enlarging the file
% unnecessarily
subsample_locs = unique([1 round(logspace(0, log10(n_steps), 300)) n_steps]);

% make plot starting with random initial conditions
ax = axes;
ax.OuterPosition = [0 0 0.5 1];
hold on;
% plot final initial conditions as predicted from infomax
plot(n_steps*ones(size(K2_infomax)), K2_infomax, ...
    'd', 'markerfacecolor', col2, 'markeredgecolor', 'none', ...
    'markersize', 3);
% plot the actual trajectory
plot(subsample_locs, dyn_history_from_random.K(:, subsample_locs)', ...
    'k', 'linewidth', 0.25);
% log scale to see all the timescales involved
set(gca, 'xscale', 'log');

xlabel('Time (a.u.)');
ylabel('K_a');

xlim([0 round(n_steps*1.1)]);
ylim([0 yl]);

beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'ticksize', 12);

% sanity check: how far are we from maximum info
disp('Random starting point');
disp(['Maximum info found with calculate_optimal_dist: ' num2str(info2(K2_infomax), '%.4f') ...
    ' run_dyn_model: ' num2str(info2(dyn_history_from_random.K(:, end)), '%.4f')]);
disp(['(Info at initial point: ' num2str(info2(Kini_random), '%.4f') ')']);
disp(' ');

% make plot starting with distribution optimized for environment 1
ax = axes;
ax.OuterPosition = [0.5 0 0.5 1];
hold on;
% plot optimal distribution in environment 1, as predicted from infomax
% (this is the same as initial conditions, up to the perturbation we added)
plot(ones(size(K1_infomax)), K1_infomax, ...
    'd', 'markerfacecolor', col1, 'markeredgecolor', 'none', ...
    'markersize', 3);
% plot final initial conditions as predicted from infomax
plot(n_steps*ones(size(K2_infomax)), K2_infomax, ...
    'd', 'markerfacecolor', col2, 'markeredgecolor', 'none', ...
    'markersize', 3);
% plot the actual trajectory
plot(subsample_locs, dyn_history_from_env1.K(:, subsample_locs)', 'k', 'linewidth', 0.25);
% log scale to see all the timescales involved
set(gca, 'xscale', 'log');

xlabel('Time (a.u.)');
ylabel('K_a');

xlim([0 round(n_steps*1.1)]);
ylim([0 yl]);

beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'ticksize', 12);

% sanity check: how far are we from maximum info
disp('Starting from optimum for different environment');
disp(['Maximum info found with calculate_optimal_dist: ' num2str(info2(K2_infomax), '%.4f') ...
    ' run_dyn_model: ' num2str(info2(dyn_history_from_env1.K(:, end)), '%.4f')]);
disp(['(Info at initial point: ' num2str(info2(Kini_env1), '%.4f') ')']);

preparegraph('edge', 0);

% make sure we're saving vector graphics, not rasterized
fig.Renderer = 'painters';

safeprint(fullfile('figs', 'dynamics_example'));
