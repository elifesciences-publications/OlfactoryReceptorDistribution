% Look at how the tuning of receptors and the kind of environment changes
% affect the optimal receptor distribution.

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   Ktot:
%       Total number of neurons.
%   cov_factor:
%       The concentrations in the artificial environment are chosen
%       arbitrarily. We roughly fit their scale to obtain intermediate SNR
%       at reasonable values of the total number of neurons, Ktot.
%   n_env_samples:
%       Number of pairs of environments to generate for each setting in
%       which we're checking how the type of environmental change affects
%       the magnitude of the change in receptor distribution.

setdefault('Ktot', 25000);
setdefault('cov_factor', 1e-4);
setdefault('n_env_samples', 500);

%% Load Hallem&Carlson sensing data

sensing_fly = open('data/flyResponsesWithNames.mat');
S_fly = sensing_fly.relRates';

% normalize by standard deviations of background rates
S_fly_normalized = bsxfun(@rdivide, S_fly, sensing_fly.bkgStd');

optim_args = {'optimopts', ...
        {'MaxFunctionEvaluations', 50000, 'Display', 'notify-detailed'}, ...
    'method', 'lagsearch'};

% generate colormap for covariance matrix plots
cmap_covmat = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

%% Generate two random environment covariances

% random environments, but make them reproducible
rng(2334);

Gamma1_generic =  cell(1, n_env_samples);
Gamma2_generic =  cell(1, n_env_samples);

K1_generic = cell(1, n_env_samples);
K2_generic = cell(1, n_env_samples);

progress = TextProgress('generating generic environments');
for i = 1:n_env_samples 
    Gamma1_generic{i} = cov_factor * generate_environment('rnd_corr', size(S_fly, 2));
    Gamma2_generic{i} = cov_factor * generate_environment('rnd_corr', size(S_fly, 2));
    
    K1_generic{i} = calculate_optimal_dist(S_fly_normalized, Gamma1_generic{i}, ...
        Ktot, optim_args{:});
    K2_generic{i} = calculate_optimal_dist(S_fly_normalized, Gamma2_generic{i}, ...
        Ktot, optim_args{:});
    
    progress.update(i*100/n_env_samples);
end
progress.done('done!');

%% Generate two environments with non-overlapping dominant odors

% random environments, but make them reproducible
rng(9854);

Gamma1_nonoverlapping =  cell(1, n_env_samples);
Gamma2_nonoverlapping =  cell(1, n_env_samples);

K1_nonoverlapping = cell(1, n_env_samples);
K2_nonoverlapping = cell(1, n_env_samples);

masks = cell(1, n_env_samples);

progress = TextProgress('generating non-overlapping environments');
for i = 1:n_env_samples
    Gamma1_nonoverlapping{i} = cov_factor * generate_environment('rnd_corr', size(S_fly, 2));
    Gamma2_nonoverlapping{i} = cov_factor * generate_environment('rnd_corr', size(S_fly, 2));
    
    % ensure little overlap between odors in env1 and env2
    mask = true(size(S_fly, 2), 1);
%     mask(1:floor(end/2)) = false;
    mask(randperm(length(mask), floor(length(mask)/2))) = false;
    masks{i} = mask;
    % to ensure we maintain positive-definiteness, work on the square roots of
    % the environment matrices
    sqrtGamma1_nonoverlapping = sqrtm(Gamma1_nonoverlapping{i});
    sqrtGamma2_nonoverlapping = sqrtm(Gamma2_nonoverlapping{i});
    % reduce the covariances of half of the odors by a factor of 4; first half
    % for environment 1, second half for environment 2
    sqrtGamma1_nonoverlapping(:, mask) = sqrtGamma1_nonoverlapping(:, mask) / 4;
    sqrtGamma2_nonoverlapping(:, ~mask) = sqrtGamma2_nonoverlapping(:, ~mask) / 4;
    % reconstruct covariance matrices from square roots
    Gamma1_nonoverlapping{i} = sqrtGamma1_nonoverlapping'*sqrtGamma1_nonoverlapping;
    Gamma2_nonoverlapping{i} = sqrtGamma2_nonoverlapping'*sqrtGamma2_nonoverlapping;
    
    % using same value for Ktot as for the generic environments
    K1_nonoverlapping{i} = calculate_optimal_dist(S_fly_normalized, Gamma1_nonoverlapping{i}, ...
        Ktot, optim_args{:});
    K2_nonoverlapping{i} = calculate_optimal_dist(S_fly_normalized, Gamma2_nonoverlapping{i}, ...
        Ktot, optim_args{:});
    
    progress.update(i*100/n_env_samples);
end
progress.done('done!');

%% Compare the change in environments

fig = figure;
fig.Color = [1 1 1];

diff_generic = arrayfun(@(i) Gamma1_generic{i}(:) - Gamma2_generic{i}(:), ...
    1:n_env_samples, 'uniform', false);
diff_nonoverlapping = arrayfun(@(i) Gamma1_nonoverlapping{i}(:) - Gamma2_nonoverlapping{i}(:), ...
    1:n_env_samples, 'uniform', false);

pooled_range = quantile(flatten(cell2mat(diff_generic)), [0.025 0.975]);

cov_bins = linspace(pooled_range(1), pooled_range(2), 100);

hold on;
for i = 1:n_env_samples
    histogram(diff_generic{i}(:), cov_bins, ...
        'edgecolor', [0.933 0.518 0.204], 'edgealpha', 0.1, 'linewidth', 0.5, ...
        'normalization', 'pdf', 'displaystyle', 'stairs');
    histogram(diff_nonoverlapping{i}(:), cov_bins, ...
        'edgecolor', [0.177, 0.459, 0.733], 'edgealpha', 0.1, 'linewidth', 0.5, ...
        'normalization', 'pdf', 'displaystyle', 'stairs');
end

h1 = histogram(diff_generic{1}(:), cov_bins, ...
    'edgecolor', [0.933 0.518 0.204], 'linewidth', 1, ...
    'normalization', 'pdf', 'displaystyle', 'stairs');
h2 = histogram(diff_nonoverlapping{1}(:), cov_bins, ...
    'edgecolor', [0.177, 0.459, 0.733], 'linewidth', 1, ...
    'normalization', 'pdf', 'displaystyle', 'stairs');

xlabel('\Delta\Gamma_{ij}');
ylabel('PDF');

legend([h1 h2], {'generic', 'non-overlapping'});
legend('boxoff');

beautifygraph('linewidth', 0.5);

preparegraph;

%% compare \Delta K in aggregate

fig = figure;
fig.Color = [1 1 1];

diff_K_generic = arrayfun(@(i) K1_generic{i}(:) - K2_generic{i}(:), ...
    1:n_env_samples, 'uniform', false);
diff_K_nonoverlapping = arrayfun(@(i) K1_nonoverlapping{i}(:) - K2_nonoverlapping{i}(:), ...
    1:n_env_samples, 'uniform', false);

pooled_diff_K_generic = abs(flatten(cell2mat(diff_K_generic)));
pooled_diff_K_nonoverlapping = abs(flatten(cell2mat(diff_K_nonoverlapping)));

pooled_K_range = quantile(pooled_diff_K_generic, [0.025 0.975]);

K_bins = linspace(pooled_K_range(1), pooled_K_range(2), 100);

hold on;
h1 = histogram(pooled_diff_K_generic(:), K_bins, ...
    'edgecolor', [0.933 0.518 0.204], 'linewidth', 1, ...
    'normalization', 'pdf', 'displaystyle', 'stairs');
h2 = histogram(pooled_diff_K_nonoverlapping(:), K_bins, ...
    'edgecolor', [0.177, 0.459, 0.733], 'linewidth', 1, ...
    'normalization', 'pdf', 'displaystyle', 'stairs');

xlabel('abs(\DeltaK_i)');
ylabel('PDF');

legend([h1 h2], {'generic', 'non-overlapping'});
legend('boxoff');

beautifygraph('linewidth', 0.5);

preparegraph;

[~, ks_p] = kstest2(pooled_diff_K_generic, pooled_diff_K_nonoverlapping, ...
    'tail', 'larger');
rs_p = ranksum(pooled_diff_K_generic, pooled_diff_K_nonoverlapping, ...
    'tail', 'left');

disp(['KS test p = ' num2str(ks_p, '%g') '.']);
disp(['Wilcoxon rank sum test p = ' num2str(rs_p, '%g') '.']);

%% Save the data

save(fullfile('save', 'change_nonoverlapping.mat'), 'Gamma1_generic', 'Gamma2_generic', ...
    'Gamma1_nonoverlapping', 'Gamma2_nonoverlapping', 'masks', ...
    'K1_generic', 'K2_generic', 'K1_nonoverlapping', 'K2_nonoverlapping', ...
    'Ktot', 'S_fly', 'S_fly_normalized', ...
    'n_env_samples', 'optim_args', 'sensing_fly');

%% SCRATCH

%% compare an example \Delta K

dist_plots = {{'env_change_with_overlap', K1_generic{1}, K2_generic{1}}, ...
    {'env_change_no_overlap', K1_nonoverlapping{1}, K2_nonoverlapping{1}}};
% , ...
%     {'env_change_narrow_tuning', K1_generic_narrow, K2_generic_narrow}, ...
%     {'env_change_wide_tuning', K1_generic_wide, K2_generic_wide}};

dist_lim = 0.06;

for i = 1:length(dist_plots)
    crt_K1 = dist_plots{i}{2};
    crt_K2 = dist_plots{i}{3};
    
    fig = figure;
    fig.Units = 'inches';
    fig.Position = [fig.Position(1:2) 3 0.8];
    
    fig.Color = [1 1 1];
    
    ax = axes;
    
    plotDistChange(crt_K1, crt_K2, 'method', 'bar', 'receptornames', sensing_fly.orNames, ...
        'labelrot', 60, 'colors', {[0.177, 0.459, 0.733], [0.933 0.518 0.204]});
    ylim([-dist_lim dist_lim]);
    hold on;
    plot(xlim, [0, 0], 'k', 'linewidth', 0.25);
    set(ax, 'ycolor', 'none');

    text(-1.8, 0, '\Delta K_a', 'horizontalalignment', 'center', 'fontsize', 8);
    
    beautifygraph('fontscale', 0.667, 'linewidth', 0.5);
    
    preparegraph('edge', 0);
    
%     safeprint(fullfile('figs', dist_plots{i}{1}));
end
