% check the effects of a change in environment

%% Load fly sensing matrix

% fly sensing matrix
datafile = open('flyResponsesWithNames.mat');
data = datafile.relRates;

S_fly = data';
sigma_fly = 1000;

%% Load mouse sensing matrix

% load response curves
mouse_responses = open('data/mouse_receptor_curves.mat');
S_mouse = mouse_responses.sensing_matrix(:, :, 11);
% set missing data to 0
S_mouse(isnan(S_mouse)) = 0;
sigma_mouse = 5;

%% General set up

% set up generic optimization options
optim_options_nodisp = {'MaxFunctionEvaluations', 50000};
optim_options = [optim_options_nodisp {'Display', 'notify-detailed'}];
optim_args_nodisp = {'optimopts', optim_options_nodisp, 'method', 'lagsearch'};
optim_args = {'optimopts', optim_options, 'method', 'lagsearch'};

% generate colormap for covariance matrix plots
cmap_covmat = divergent([0.21 0.17 0.53], [0.98 0.40 0.17], 256);

%% Generating artificial sensing matrices

% random narrowly-tuned sensing matrix
% make this reproducible
rng(3874);

% random widely-tuned sensing matrix
% this is only used for its norm later in the code
% XXX I should get rid of this and use something based on S_fly instead
S_wide_legacy = generate_random_sensing(24, 50, [0.2 1.0], 200);

%% Example using fly sensing matrix

% make this reproducible
rng(9843857);

% start from a random environment covariance matrix
Gamma0_example = generate_environment('rnd_corr', size(S_fly, 2));
% add variance to two different sets of odorants
Gamma1_example = generate_environment('delta_rnd_diag', Gamma0_example, ...
    'deltasize', 100, 'deltapos', [16, 21, 23, 25, 43, 46, 81, 91, 99, 105]);
Gamma2_example = generate_environment('delta_rnd_diag', Gamma0_example, ...
    'deltasize', 100, 'deltapos', [5, 18, 33, 36, 46, 53, 66, 71, 84, 101, 107]);

% values of Ktot for low and high SNR
Ktot_example = [300 30000];
% get optimal repertoires
[K1_example, ~, Q1_example] = calculate_optimal_dist(S_fly/sigma_fly, ...
    Gamma1_example, Ktot_example, optim_args{:});
[K2_example, ~, Q2_example] = calculate_optimal_dist(S_fly/sigma_fly, ...
    Gamma2_example, Ktot_example, optim_args{:});

%% ... plot the change

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.8 2];

% show change at low SNR
% subplot(1, 2, 2);
ax1 = axes;
ax1.Units = 'normalized';
ax1.OuterPosition = [1/2 0 1/2 1];
plotDistChange(K1_example(:, 1), K2_example(:, 1), 'beautifyopts', {'fontscale', 0.833});
ylim([0 0.08]);
title('Low SNR');

% subplot(1, 2, 1);
% show change at high SNR
ax2 = axes;
ax2.Units = 'normalized';
ax2.OuterPosition = [0 0 1/2 1];
plotDistChange(K1_example(:, 2), K2_example(:, 2), 'beautifyopts', {'fontscale', 0.833});
ylim([0 0.08]);
title('High SNR');

preparegraph('edge', 0);

safe_print(fullfile('figs', 'natural_env_change_example.pdf'));

%% ... plot covariance matrices used for the example

% we have diagrams for these odorants
odor_pic_names = {'1-hexanol', 'E2-hexenol', 'E3-hexenol', 'Z3-hexenol', ...
    '2-heptanone', 'butyl acetate', 'pentyl acetate', 'methyl hexanoate', ...
    'ethyl hexanoate', 'methyl benzoate', 'ethyl benzoate'};
odor_pic_idxs = zeros(size(odor_pic_names));
for i = 1:length(odor_pic_names)
    odor_pic_idxs(i) = find(strcmp(datafile.odorNames, odor_pic_names{i}));
end;

% choose a scale for the colormap, and locations for ticks
cticks = [-1.5, -1, -0.5, 0, 0.5, 1, 1.5];
cticklabels = arrayfun(@num2str, cticks, 'uniform', false);
clims = [-1.5 1.5];

% draw the zoomed-in detail of the covariance matrix for environment 1
fig1 = figure;
fig1.Units = 'inches';
fig1.Position = [fig1.Position(1:2) 1.9 1.9];

imagesc(Gamma1_example(odor_pic_idxs, odor_pic_idxs), clims);
colormap(cmap_covmat);

% overlay a white grid, so that each element is separated
hold on;
n_odor = length(odor_pic_idxs);
for i = 0:n_odor
    plot([0.5, n_odor + 0.5], i + [0.5, 0.5], 'w', 'linewidth', 2);
    plot(i + [0.5, 0.5], [0.5, n_odor + 0.5], 'w', 'linewidth', 2);
end

beautifygraph;

axis equal;
axis off;

preparegraph;

safe_print(fullfile('figs', 'odor_cov_example_env1_zoom.pdf'));

% draw the zoomed-in detail of the covariance matrix for environment 1
fig2 = figure;
fig2.Units = 'inches';
fig2.Position = [fig2.Position(1:2) 1.9 1.9];

imagesc(Gamma2_example(odor_pic_idxs, odor_pic_idxs), clims);
colormap(cmap_covmat);

% overlay a white grid, so that each element is separated
hold on;
n_odor = length(odor_pic_idxs);
for i = 0:n_odor
    plot([0.5, n_odor + 0.5], i + [0.5, 0.5], 'w', 'linewidth', 2);
    plot(i + [0.5, 0.5], [0.5, n_odor + 0.5], 'w', 'linewidth', 2);
end

beautifygraph;

axis equal;
axis off;

preparegraph;

safe_print(fullfile('figs', 'odor_cov_example_env2_zoom.pdf'));

% draw the color bar
fig3 = figure;
fig3.Units = 'inches';
fig3.Position = [fig3.Position(1:2) 1.5 1.9];

beautifygraph;

axis off;
caxis(clims);
colormap(cmap_covmat);
hcb = colorbar('ytick', cticks, 'yticklabel', cticklabels);

preparegraph;

safe_print(fullfile('figs', 'odor_cov_example_zoom_colorbar.pdf'));

% draw the full background covariance matrix, Gamma0_example
fig4 = figure;
fig4.Units = 'inches';
fig4.Position = [fig4.Position(1:2) 1.9 1.9];

% reorder the indices so that the ones that we have diagrams for are
% bunched together in the middle (that was it's easier to show how we zoom
% in)
reorder = [1:min(odor_pic_idxs)-1 odor_pic_idxs ...
    setdiff(setdiff(1:size(Gamma0_example, 1), odor_pic_idxs), 1:min(odor_pic_idxs)-1)];
box_start = min(odor_pic_idxs) - 0.5;
box_end = min(odor_pic_idxs)+length(odor_pic_idxs)-0.5;

% plotMatrix(Gamma0_ex(reorder, reorder), [-0.5 0.5]);
plotMatrix(Gamma0_example(reorder, reorder), clims);
colormap(cmap_covmat);

% draw a box showing where the zoomed-in matrices are located
hold on;
% draw a thick outline in white, and a thin one in black, to improve
% visibility
plot([box_start box_start box_end box_end box_start], ...
     [box_start box_end box_end box_start box_start], 'w', 'linewidth', 2);
plot([box_start box_start box_end box_end box_start], ...
     [box_start box_end box_end box_start box_start], 'k');

beautifygraph('fontscale', 0.833);

axis equal;
set(gca, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');

preparegraph;

safe_print(fullfile('figs', 'odor_cov_example_env0.pdf'));

%% Effect of removing weak responses

S_fly_no_weak = S_fly;
% mask = abs(S_fly_nw) > quantile(abs(S_fly_nw(:)), 0.5);
mask = abs(S_fly_no_weak) > 50;
S_fly_no_weak(~mask) = 0;

[K1_example_no_weak, ~, Q1_example_no_weak] = calculate_optimal_dist(...
    S_fly_no_weak/sigma_fly, Gamma1_example, Ktot_example, optim_args{:});

plotDistChange(K1_example(:, 1), K1_example_no_weak(:, 1));

%% Distribution change upon environment change, as a function of tuning width

% random environments and sensing matrices, but make them reproducible
rng(135753);

sensitivity_vs_tuning.Gamma1 = generate_environment(...
    'rnd_corr', 50, 'corrbeta', 50);
sensitivity_vs_tuning.Gamma2 = generate_environment(...
    'rnd_corr', 50, 'corrbeta', 50);

% fix the total number of neurons
sensitivity_vs_tuning.Ktot = 4000;

% choose a range of tuning widths
%tuning_range = linspace(0.05, 2.0, 32);
sensitivity_vs_tuning.tuning_range = logspace(log10(0.01), log10(0.2), 48);
% interpolate sigmas to stay in a regime of intermediate SNR
sensitivity_vs_tuning.sigmas = logspace(log10(50), log10(0.6), ...
    length(sensitivity_vs_tuning.tuning_range));
% number of samples for each tuning width
sensitivity_vs_tuning.n_per_tuning = 24;
cell_empty = cell(sensitivity_vs_tuning.n_per_tuning, ...
    length(sensitivity_vs_tuning.tuning_range));
sensitivity_vs_tuning.S = cell_empty;
sensitivity_vs_tuning.K1 = cell_empty;
sensitivity_vs_tuning.Q1 = cell_empty;
sensitivity_vs_tuning.K2 = cell_empty;
sensitivity_vs_tuning.Q2 = cell_empty;

tic;
% normalizing to S_wide_legacy for historical reasons (that's where we
% found the values for sigma_narrow, sigma_wide)
normalize_to_S_wide_legacy = @(S) S*norm(S_wide_legacy)/norm(S);
for i = 1:length(sensitivity_vs_tuning.tuning_range)
    disp([int2str(i) '/' int2str(length(sensitivity_vs_tuning.tuning_range)) '...']);
    crt_sigma = sensitivity_vs_tuning.sigmas(i);
    for j = 1:sensitivity_vs_tuning.n_per_tuning
        % generate a random sensing matrix and normalize
        crt_S = normalize_to_S_wide_legacy(generate_random_sensing(24, 50, ...
            sensitivity_vs_tuning.tuning_range(i), 200));
        sensitivity_vs_tuning.S{j, i} = crt_S;
        
        % get the implied optimal distributions
        [sensitivity_vs_tuning.K1{j, i}, ~, sensitivity_vs_tuning.Q1{j, i}] = ...
            calculate_optimal_dist(crt_S/crt_sigma, ...
            sensitivity_vs_tuning.Gamma1, sensitivity_vs_tuning.Ktot, ...
            optim_args_nodisp{:});
        [sensitivity_vs_tuning.K2{j, i}, ~, sensitivity_vs_tuning.Q2{j, i}] = ...
            calculate_optimal_dist(crt_S/crt_sigma, ...
            sensitivity_vs_tuning.Gamma2, sensitivity_vs_tuning.Ktot, ...
            optim_args_nodisp{:});
    end
end
disp(['Generating tuning width sweep took ' num2str(toc, '%.2f') ' seconds.']);

%% ... calculate correlation coefficients

zero_empty = zeros(sensitivity_vs_tuning.n_per_tuning, ...
    length(sensitivity_vs_tuning.tuning_range));
sensitivity_vs_tuning.corrs = struct;
sensitivity_vs_tuning.corrs.cs_Q_abs = zero_empty;
sensitivity_vs_tuning.corrs.cp_Q_abs = zero_empty;
sensitivity_vs_tuning.corrs.cs_Q_diff = zero_empty;
sensitivity_vs_tuning.corrs.cp_Q_diff = zero_empty;

sensitivity_vs_tuning.corrs.cs_logQ_abs = zero_empty;
sensitivity_vs_tuning.corrs.cp_logQ_abs = zero_empty;
sensitivity_vs_tuning.corrs.cs_logQ_diff = zero_empty;
sensitivity_vs_tuning.corrs.cp_logQ_diff = zero_empty;

sensitivity_vs_tuning.corrs.cs_invQ_abs = zero_empty;
sensitivity_vs_tuning.corrs.cp_invQ_abs = zero_empty;
sensitivity_vs_tuning.corrs.cs_invQ_diff = zero_empty;
sensitivity_vs_tuning.corrs.cp_invQ_diff = zero_empty;

sensitivity_vs_tuning.corrs.cs_invQpred_abs = zero_empty;
sensitivity_vs_tuning.corrs.cp_invQpred_abs = zero_empty;
sensitivity_vs_tuning.corrs.cs_invQpred_diff = zero_empty;
sensitivity_vs_tuning.corrs.cp_invQpred_diff = zero_empty;

for i = 1:length(sensitivity_vs_tuning.tuning_range)
    disp([int2str(i) '/' int2str(length(sensitivity_vs_tuning.tuning_range)) '...']);
    crt_sigma = sensitivity_vs_tuning.sigmas(i);
    for j = 1:sensitivity_vs_tuning.n_per_tuning
        crt_Q1 = sensitivity_vs_tuning.Q1{j, i};
        crt_Q2 = sensitivity_vs_tuning.Q2{j, i};
        
        crt_K1 = sensitivity_vs_tuning.K1{j, i};
        crt_K2 = sensitivity_vs_tuning.K2{j, i};
        
        % correlation with diagonal elements of Q
        sensitivity_vs_tuning.corrs.cs_Q_abs(j, i) = corr(...
            diag(crt_Q1), crt_K1, 'type', 'spearman');
        sensitivity_vs_tuning.corrs.cp_Q_abs(j, i) = corr(...
            diag(crt_Q1), crt_K1);
        
        % after taking the log
        sensitivity_vs_tuning.corrs.cs_logQ_abs(j, i) = corr(...
            log(diag(crt_Q1)), crt_K1, 'type', 'spearman');
        sensitivity_vs_tuning.corrs.cp_logQ_abs(j, i) = corr(...
            log(diag(crt_Q1)), crt_K1);
        
        crt_invQ1 = inv(crt_Q1);
        crt_invQ2 = inv(crt_Q2);
        
        % correlation with diagonal elements of inv(Q)
        sensitivity_vs_tuning.corrs.cs_invQ_abs(j, i) = corr(...
            diag(crt_invQ1), crt_K1, 'type', 'spearman');
        sensitivity_vs_tuning.corrs.cp_invQ_abs(j, i) = corr(...
            diag(crt_invQ1), crt_K1);
        
        crt_invQ1pred = max(...
            sensitivity_vs_tuning.Ktot / size(sensitivity_vs_tuning.S{j, i}, 1) + ...
            mean(diag(crt_invQ1)) - diag(crt_invQ1), 0);
        crt_invQ2pred = max(...
            sensitivity_vs_tuning.Ktot / size(sensitivity_vs_tuning.S{j, i}, 1) + ...
            mean(diag(crt_invQ2)) - diag(crt_invQ2), 0);
        
        % correlation with clipped version of diagonal elements of inv(Q)
        sensitivity_vs_tuning.corrs.cs_invQpred_abs(j, i) = corr(...
            crt_invQ1pred, crt_K1, 'type', 'spearman');
        sensitivity_vs_tuning.corrs.cp_invQpred_abs(j, i) = corr(...
            crt_invQ1pred, crt_K1);

        % correlations between differences
        crt_K_diff = crt_K2 - crt_K1;
        
        crt_Qd_diff = diag(crt_Q2) - diag(crt_Q1);
        crt_logQd_diff = log(diag(crt_Q2)) - log(diag(crt_Q1));
        crt_invQd_diff = diag(crt_invQ2) - diag(crt_invQ1);
        crt_invQdpred_diff = crt_invQ2pred - crt_invQ1pred;
        
        sensitivity_vs_tuning.corrs.cs_Q_diff(j ,i) = corr(...
            crt_Qd_diff, crt_K_diff, 'type', 'spearman');
        sensitivity_vs_tuning.corrs.cp_Q_diff(j, i) = corr(...
            crt_Qd_diff, crt_K_diff);
        
        sensitivity_vs_tuning.corrs.cs_logQ_diff(j, i) = corr(...
            crt_logQd_diff, crt_K_diff, 'type', 'spearman');
        sensitivity_vs_tuning.corrs.cp_logQ_diff(j, i) = corr(...
            crt_logQd_diff, crt_K_diff);
        
        sensitivity_vs_tuning.corrs.cs_invQ_diff(j, i) = corr(...
            crt_invQd_diff, crt_K_diff, 'type', 'spearman');
        sensitivity_vs_tuning.corrs.cp_invQ_diff(j, i) = corr(...
            crt_invQd_diff, crt_K_diff);

        sensitivity_vs_tuning.corrs.cs_invQpred_diff(j, i) = corr(...
            crt_invQdpred_diff, crt_K_diff, 'type', 'spearman');
        sensitivity_vs_tuning.corrs.cp_invQpred_diff(j, i) = corr(...
            crt_invQdpred_diff, crt_K_diff);
    end
end

% summarize (calculate means, standard deviations, etc.)
fields = fieldnames(sensitivity_vs_tuning.corrs);
for i = 1:length(fields)
    crt_corrs = sensitivity_vs_tuning.corrs.(fields{i});
    
    sensitivity_vs_tuning.corrs.summary.([fields{i} '_mean']) = ...
        mean(crt_corrs, 1);
    sensitivity_vs_tuning.corrs.summary.([fields{i} '_median']) = ...
        median(crt_corrs, 1);

    sensitivity_vs_tuning.corrs.summary.([fields{i} '_std']) = ...
        std(crt_corrs, [], 1);
    
    sensitivity_vs_tuning.corrs.summary.([fields{i} '_low']) = ...
        quantile(crt_corrs, 0.2, 1);
    sensitivity_vs_tuning.corrs.summary.([fields{i} '_high']) = ...
        quantile(crt_corrs, 0.8, 1);
end

%% ... show how well log(diag(Q)) and diag(inv(Q)) predict receptor distribution

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 3 1.8];

ax = axes;
ax.OuterPosition = [0 0 0.5 1];

fill([flipud(sensitivity_vs_tuning.tuning_range(:)) ; ...
        sensitivity_vs_tuning.tuning_range(:)], ...
    [flipud(sensitivity_vs_tuning.corrs.summary.cp_logQ_abs_low(:)) ; ...
        sensitivity_vs_tuning.corrs.summary.cp_logQ_abs_high(:)], ...
    [0.9 0.9 0.9], 'linestyle', 'none');
hold on;
plot(sensitivity_vs_tuning.tuning_range, ...
    sensitivity_vs_tuning.corrs.summary.cp_logQ_abs_mean, 'linewidth', 1);

tuning_min_max = [min(sensitivity_vs_tuning.tuning_range) max(sensitivity_vs_tuning.tuning_range)];
xlim(tuning_min_max);
fixed_y_lim = [-0.1 1];
ylim(fixed_y_lim);

set(gca, 'xtick', tuning_min_max, 'xticklabel', {'narrow', 'wide'}, ...
    'xscale', 'log');
xlabel('receptor tuning');
ylabel('corr(log Q_{aa}, K_a)');

plot(xlim, [0 0], 'k:');

beautifygraph('fontscale', 0.833);

ax = axes;
ax.OuterPosition = [0.5 0 0.5 1];

fill([flipud(sensitivity_vs_tuning.tuning_range(:)) ; ...
        sensitivity_vs_tuning.tuning_range(:)], ...
    -[flipud(sensitivity_vs_tuning.corrs.summary.cp_invQ_abs_low(:)) ; ...
        sensitivity_vs_tuning.corrs.summary.cp_invQ_abs_high(:)], ...
    [0.9 0.9 0.9], 'linestyle', 'none');
hold on;
plot(sensitivity_vs_tuning.tuning_range, ...
    -sensitivity_vs_tuning.corrs.summary.cp_invQ_abs_mean, 'linewidth', 1);

tuning_min_max = [min(sensitivity_vs_tuning.tuning_range) max(sensitivity_vs_tuning.tuning_range)];
xlim(tuning_min_max);
ylim(fixed_y_lim);

set(gca, 'xtick', tuning_min_max, 'xticklabel', {'narrow', 'wide'}, ...
    'xscale', 'log');
xlabel('receptor tuning');
ylabel('corr(-Q^{-1}_{aa}, K_a)');

plot(xlim, [0 0], 'k:');

beautifygraph('fontscale', 0.833);
preparegraph;

safe_print(fullfile('figs', 'logic_vs_tuning.pdf'));

%% ... how well do log(diag(Q)), diag(inv(Q)) predict distribution changes

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 3 1.8];

ax = axes;
ax.OuterPosition = [0 0 0.5 1];

fill([flipud(sensitivity_vs_tuning.tuning_range(:)) ; ...
        sensitivity_vs_tuning.tuning_range(:)], ...
    [flipud(sensitivity_vs_tuning.corrs.summary.cp_logQ_diff_low(:)) ; ...
        sensitivity_vs_tuning.corrs.summary.cp_logQ_diff_high(:)], ...
    [0.9 0.9 0.9], 'linestyle', 'none');
hold on;
plot(sensitivity_vs_tuning.tuning_range, ...
    sensitivity_vs_tuning.corrs.summary.cp_logQ_diff_mean, 'linewidth', 1);

tuning_min_max = [min(sensitivity_vs_tuning.tuning_range) max(sensitivity_vs_tuning.tuning_range)];
xlim(tuning_min_max);
fixed_y_lim = [-0.1 1];
ylim(fixed_y_lim);

set(gca, 'xtick', tuning_min_max, 'xticklabel', {'narrow', 'wide'}, ...
    'xscale', 'log');
xlabel('receptor tuning');
ylabel('corr(\Delta log Q_{aa}, \Delta K_a)');

plot(xlim, [0 0], 'k:');

beautifygraph('fontscale', 0.833);

ax = axes;
ax.OuterPosition = [0.5 0 0.5 1];

fill([flipud(sensitivity_vs_tuning.tuning_range(:)) ; ...
        sensitivity_vs_tuning.tuning_range(:)], ...
    -[flipud(sensitivity_vs_tuning.corrs.summary.cp_invQ_diff_low(:)) ; ...
        sensitivity_vs_tuning.corrs.summary.cp_invQ_diff_high(:)], ...
    [0.9 0.9 0.9], 'linestyle', 'none');
hold on;
plot(sensitivity_vs_tuning.tuning_range, ...
    -sensitivity_vs_tuning.corrs.summary.cp_invQ_diff_mean, 'linewidth', 1);

tuning_min_max = [min(sensitivity_vs_tuning.tuning_range) max(sensitivity_vs_tuning.tuning_range)];
xlim(tuning_min_max);
ylim(fixed_y_lim);

set(gca, 'xtick', tuning_min_max, 'xticklabel', {'narrow', 'wide'}, ...
    'xscale', 'log');
xlabel('receptor tuning');
ylabel('corr(-\Delta Q^{-1}_{aa}, \Delta K_a)');

plot(xlim, [0 0], 'k:');

beautifygraph('fontscale', 0.833);
preparegraph;

safe_print(fullfile('figs', 'logic_vs_tuning_envchange.pdf'));

%% ... K vs. diag(inv(Q)) at wide tuning

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2.9 2];

idx_rep = 1;
idx_tuning = 48;
% idx_tuning = 10;
crt_Q = sensitivity_vs_tuning.Q1{idx_rep, idx_tuning};
crt_K = sensitivity_vs_tuning.K1{idx_rep, idx_tuning};
crt_invQ = inv(crt_Q);
crt_invQd = diag(crt_invQ);
approx_const = sensitivity_vs_tuning.Ktot/size(sensitivity_vs_tuning.S{1, 1}, 1) ...
    + mean(crt_invQd);
% crt_invQd_clipped = min(crt_invQd, approx_const);
scatterfit(crt_invQd, crt_K, 'scatteropts', {60, 'filled'}, 'fitopts', ...
    {'line', [-1, approx_const], ...
    'style', {'k--', 'linewidth', 1}, 'legend', false});

% ylim([0 max(crtK)]);

xlabel('(Q^{-1})_{aa}');
ylabel('K_a');

beautifygraph('fontscale', 0.833);
preparegraph;

safe_print(fullfile('figs', 'K_vs_invQ.pdf'));

%% Look at effects of generic environment change

% random environments, but make them reproducible
rng(2334);
 
Gamma1_generic = generate_environment('rnd_corr', size(S_fly, 2));
Gamma2_generic = generate_environment('rnd_corr', size(S_fly, 2));

%% ... search for Ktot such that most receptors in both environments are "on"

test_Ktot_generic = logspace(log10(0.1), log10(50000), 32);

% using the fly data here
tic;
test_K1_generic = calculate_optimal_dist(S_fly/sigma_fly, Gamma1_generic, ...
    test_Ktot_generic, optim_args{:});
test_K2_generic = calculate_optimal_dist(S_fly/sigma_fly, Gamma2_generic, ...
    test_Ktot_generic, optim_args{:});
disp(['Calculating optimal distributions for two environments at ' ...
    int2str(length(test_Ktot_generic)) ' Ktot values took ' num2str(toc, '%.1f') ...
    ' seconds.']);

% find numbers of receptor types that are populated above the 0.1% mark
pop_threshold = 0.001;
test_count1_generic = sum(bsxfun(@ge, test_K1_generic, pop_threshold*test_Ktot_generic), 1);
test_count2_generic = sum(bsxfun(@ge, test_K2_generic, pop_threshold*test_Ktot_generic), 1);

%% ... get the optimal repertoires at a chosen Ktot

% choosing a Ktot manually based on results from the previous cell
Ktot_generic = 25000;

% show where we are on the plot of number of receptors expressed vs.
% number of OSN
semilogx(test_Ktot_generic, test_count1_generic, 'linewidth', 1);
hold all;
semilogx(test_Ktot_generic, test_count2_generic, 'linewidth', 1);
plot([25000 25000], [0 25], 'k--');
text(25000, 10, 'chosen K_{tot}  ', 'horizontalalignment', 'right');

xlabel('Total number of OSN');
ylabel('Number of expressed receptor types');

% now calculate the optimal repertoires at the chosen K_tot
K1_generic = calculate_optimal_dist(S_fly/sigma_fly, Gamma1_generic, ...
    Ktot_generic, optim_args{:});
K2_generic = calculate_optimal_dist(S_fly/sigma_fly, Gamma2_generic, ...
    Ktot_generic, optim_args{:});

%% ... compare the repertoires

figure;
scatterfit(K1_generic, K2_generic, 'scatteropts', {300, '.r'});
lims = [min(min(K1_generic), min(K2_generic)) - 50 max(max(K1_generic), max(K2_generic)) + 50];
xlim(lims);
ylim(lims);
hold on;
drawfitline(lims, lims, 'line', [1 0], 'style', {':k'}, 'legend', false);
title('Environments from the same ensemble');

xlabel('Receptor abundances in environment 1');
ylabel('Receptor abundances in environment 2');

beautifygraph;
preparegraph;

figure;
plotDistChange(K1_generic, K2_generic, 'beautifyopts', {'fontscale', 0.833});
ylim([0 max([K1_generic(:) ; K2_generic(:)])/Ktot_generic]);

%% Look at effects of changing the set of odors that are present

% random environments, but make them reproducible
rng(9854);

Gamma1_nonoverlapping = generate_environment('rnd_corr', size(S_fly, 2));
Gamma2_nonoverlapping = generate_environment('rnd_corr', size(S_fly, 2));

% ensure little overlap between odors in env1 and env2
mask = true(size(Gamma1_nonoverlapping, 1), 1);
mask(1:floor(end/2)) = false;
% to ensure we maintain positive-definiteness, work on the square roots of
% the environment matrices
sqrtGamma1_nonoverlapping = sqrtm(Gamma1_nonoverlapping);
sqrtGamma2_nonoverlapping = sqrtm(Gamma2_nonoverlapping);
% reduce the covariances of half of the odors by a factor of 4; first half
% for environment 1, second half for environment 2
sqrtGamma1_nonoverlapping(:, mask) = sqrtGamma1_nonoverlapping(:, mask) / 4;
sqrtGamma2_nonoverlapping(:, ~mask) = sqrtGamma2_nonoverlapping(:, ~mask) / 4;
% reconstruct covariance matrices from square roots
Gamma1_nonoverlapping = sqrtGamma1_nonoverlapping'*sqrtGamma1_nonoverlapping;
Gamma2_nonoverlapping = sqrtGamma2_nonoverlapping'*sqrtGamma2_nonoverlapping;

% using same value for Ktot as for the generic environments
K1_nonoverlapping = calculate_optimal_dist(S_fly/sigma_fly, Gamma1_nonoverlapping, ...
    Ktot_generic, optim_args{:});
K2_nonoverlapping = calculate_optimal_dist(S_fly/sigma_fly, Gamma2_nonoverlapping, ...
    Ktot_generic, optim_args{:});

%% ... and compare
figure;
scatterfit(K1_nonoverlapping, K2_nonoverlapping, 'scatteropts', {300, '.r'});
lims = [min(min(K1_nonoverlapping), min(K2_nonoverlapping)) - 50 max(max(K1_nonoverlapping), max(K2_nonoverlapping)) + 50];
xlim(lims);
ylim(lims);
hold on;
drawfitline(lims, lims, 'line', [1 0], 'style', {':k'}, 'legend', false);
title('Environments with non-overlapping odorants');

xlabel('Receptor abundances in environment 1');
ylabel('Receptor abundances in environment 2');

beautifygraph;
preparegraph;

figure;

plotDistChange(K1_nonoverlapping, K2_nonoverlapping, 'beautifyopts', {'fontscale', 0.833});
ylim([0 max([K1_nonoverlapping(:) ; K2_nonoverlapping(:)])/Ktot_generic]);

%% ... make plots for random change vs. non-overlapping environments

% rng_perc = 0.95;
rng_perc = 0.9965;
gRng = quantile([abs(Gamma1_generic(:)) ; abs(Gamma2_generic(:))], rng_perc);
gRng_no = quantile([abs(Gamma1_nonoverlapping(:)) ; abs(Gamma2_nonoverlapping(:))], rng_perc);
gRng_all = max(gRng, gRng_no);

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.8 2.6];

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.05 0.35 0.2 0.65];
plotMatrix(Gamma1_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 1', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.25 0.35 0.2 0.65];
plotMatrix(Gamma2_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
% ylabel('odorant 2');
ylabel(' ');
title('Env. 2', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.55 0.35 0.2 0.65];
plotMatrix(Gamma1_nonoverlapping, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 1', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.75 0.35 0.2 0.65];
plotMatrix(Gamma2_nonoverlapping, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
% ylabel('odorant 2');
ylabel(' ');
title('Env. 2', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0 0 1/2 0.4];
% yl = max(abs([K1(:) - K2(:) ; K1b(:) - K2b(:)])) / Ktot;
yl = 0.06;
plotDistChange(K1_generic, K2_generic, 'method', 'bar', 'receptornames', datafile.orNames);
ylim([-yl yl]);
hold on;
plot(xlim, [0 0], 'k', 'linewidth', 0.25);
set(ax, 'ycolor', 'none');
title('Generic environments', 'fontsize', 10);
text(-1.5, 0, '\Delta K_a', 'rotation', 90, 'horizontalalignment', 'center');

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [1/2 0 1/2 0.4];
plotDistChange(K1_nonoverlapping, K2_nonoverlapping, 'method', 'bar', 'receptornames', datafile.orNames);
ylim([-yl yl]);
hold on;
plot(xlim, [0 0], 'k', 'linewidth', 0.25);
set(ax, 'ycolor', 'none');
title('Non-overlapping environments', 'fontsize', 10);
text(-1.5, 0, '\Delta K_a', 'rotation', 90, 'horizontalalignment', 'center');

preparegraph('edge', 0);

safe_print(fullfile('figs', 'env_change_vs_overlap.pdf'));

%% ... different viz: plots for random change vs. non-overlapping environments

gRng = quantile([abs(Gamma1_generic(:)) ; abs(Gamma2_generic(:))], 0.95);
gRng_no = quantile([abs(Gamma1_nonoverlapping(:)) ; abs(Gamma2_nonoverlapping(:))], 0.95);
gRng_all = max(gRng, gRng_no);

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.8 3.5];

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.05 0.65 0.2 0.35];
plotMatrix(Gamma1_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 1', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.25 0.65 0.2 0.35];
plotMatrix(Gamma2_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 2', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.55 0.65 0.2 0.35];
plotMatrix(Gamma1_nonoverlapping, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 1', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.75 0.65 0.2 0.35];
plotMatrix(Gamma2_nonoverlapping, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 2', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0 0 1/2 0.6];
plotDistChange(K1_generic, K2_generic, 'beautifyopts', {'fontscale', 0.833});
ylim([0 0.07]);
title('Generic environments', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [1/2 0 1/2 0.6];
plotDistChange(K1_nonoverlapping, K2_nonoverlapping, 'beautifyopts', {'fontscale', 0.833});
ylim([0 0.07]);
title('Non-overlapping environments', 'fontsize', 10);

preparegraph('edge', 0);

% safe_print(fullfile('figs', 'env_change_vs_overlap_different.pdf'));

%% Look at generic environments with narrow vs. wide tuning

% keep things reproducible
rng(123);

% generate narrowly- and broadly-tuned sensing matrices
S_narrow = generate_random_sensing(size(S_fly, 1), size(S_fly, 2), 0.05, 200);
S_wide = generate_random_sensing(size(S_fly, 1), size(S_fly, 2), 0.8, 200);

% normalize them in a similar way
S_narrow = S_narrow*max(S_fly(:))/max(S_narrow(:));
S_wide = S_wide*max(S_fly(:))/max(S_wide(:));

% they still require different sigmas (which are just used as
% normalization; see below) to fall in the intermediate-SNR regime
sigma_narrow = 3500;
sigma_wide = 350;

% calculate the optimal distributions at both narrow and wide tuning, in
% both environments
K1_generic_narrow = calculate_optimal_dist(S_narrow/sigma_narrow, ...
    Gamma1_generic, Ktot_generic, optim_args{:}, 'lagstart', size(S_fly, 1)/Ktot_generic);
K2_generic_narrow = calculate_optimal_dist(S_narrow/sigma_narrow, ...
    Gamma2_generic, Ktot_generic, optim_args{:}, 'lagstart', size(S_fly, 1)/Ktot_generic);

K1_generic_wide = calculate_optimal_dist(S_wide/sigma_wide, Gamma1_generic, ...
    Ktot_generic, optim_args{:}, 'lagstart', size(S_fly, 1)/Ktot_generic);
K2_generic_wide = calculate_optimal_dist(S_wide/sigma_wide, Gamma2_generic, ...
    Ktot_generic, optim_args{:}, 'lagstart', size(S_fly, 1)/Ktot_generic);

% sanity check on convergence: make sure the total number of neurons is as
% prescribed
if abs(sum(K1_generic_narrow) / Ktot_generic - 1) > 1e-4
    error('K1_generic_narrow doesn''t sum up to Ktot_generic.');
end
if abs(sum(K2_generic_narrow) / Ktot_generic - 1) > 1e-4
    error('K2_generic_narrow doesn''t sum up to Ktot_generic.');
end
if abs(sum(K1_generic_wide) / Ktot_generic - 1) > 1e-4
    error('K1_generic_wide doesn''t sum up to Ktot_generic.');
end
if abs(sum(K2_generic_wide) / Ktot_generic - 1) > 1e-4
    error('K2_generic_wide doesn''t sum up to Ktot_generic.');
end

%% ... save sensing matrix plots

% make plots of the sensing matrices used here, and save them
% (this is used in some talks)
all_mats = containers.Map;

all_mats('S_narrow') = S_narrow;
all_mats('S_wide') = S_wide;

map_keys = all_mats.keys;
for i = 1:length(map_keys)
    crt_name = map_keys{i};
    crt_S = all_mats(crt_name);

    fig = figure;
    fig.Units = 'inches';
    
    fig_x = 3.5;
    fig_y = 1;
    
    fig.Position = [fig.Position(1:2) fig_x fig_y];
    
    ax = axes;
    ax.Units = 'inches';
    
    edge_x = 0.25;
    edge_y = 0.25;
    
    margin_x = 0.1;
    
    ax_x = (fig_x - edge_x - margin_x);
    ax_y = ax_x*size(crt_S, 1)/size(crt_S, 2);
    
    ax.Position = [edge_x edge_y ax_x ax_y];
    
    imagesc(crt_S, [0 250]);
    
    beautifygraph;
    axis equal;
    box on;
    set(gca, 'xminortick', 'off', 'yminortick', 'off', 'xtick', [], 'ytick', []);
    xlabel('Odorants');
    ylabel('Receptors');
    
    preparegraph;
    
    safe_print(fullfile('figs', 'tuning_figs', [crt_name '.pdf']));
end

%% ... make plots for random change vs. non-overlapping

% rng_perc = 0.95;
rng_perc = 0.9965;
gRng_all = quantile([abs(Gamma1_generic(:)) ; abs(Gamma2_generic(:))], rng_perc);

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.8 2.6];

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.05 0.35 0.2 0.65];
plotMatrix(Gamma1_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 1', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.25 0.35 0.2 0.65];
plotMatrix(Gamma2_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
% ylabel('odorant 2');
ylabel(' ');
title('Env. 2', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.55 0.35 0.2 0.65];
plotMatrix(Gamma1_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 1', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.75 0.35 0.2 0.65];
plotMatrix(Gamma2_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
% ylabel('odorant 2');
ylabel(' ');
title('Env. 2', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0 0 1/2 0.4];
% yl = max(abs([K1g_narrow(:) - K2g_narrow(:) ; K1g_wide(:) - K2g_wide(:)])) / Ktot;
yl = 0.06;
plotDistChange(K1_generic_narrow, K2_generic_narrow, 'method', 'bar', 'receptornames', datafile.orNames);
ylim([-yl yl]);
hold on;
plot(xlim, [0 0], 'k', 'linewidth', 0.25);
set(ax, 'ycolor', 'none');
title('Narrow tuning', 'fontsize', 10);
text(-1.5, 0, '\Delta K_a', 'rotation', 90, 'horizontalalignment', 'center');

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [1/2 0 1/2 0.4];
plotDistChange(K1_generic_wide, K2_generic_wide, 'method', 'bar', 'receptornames', datafile.orNames);
ylim([-yl yl]);
hold on;
plot(xlim, [0 0], 'k', 'linewidth', 0.25);
set(ax, 'ycolor', 'none');
title('Wide tuning', 'fontsize', 10);
text(-1.5, 0, '\Delta K_a', 'rotation', 90, 'horizontalalignment', 'center');

preparegraph('edge', 0);

safe_print(fullfile('figs', 'env_change_vs_tuning.pdf'));

%% ... different viz: plots for distribution change as a function of tuning width

gRng_all = quantile([abs(Gamma1_generic(:)) ; abs(Gamma2_generic(:))], 0.95);
% gbRng = quantile([abs(Gamma1b(:)) ; abs(Gamma2b(:))], 0.98);
% gRng_all = max(gRng, gbRng);

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.8 3.5];
% fig.Position = [fig.Position(1:2) 5.8 2.1];

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.05 0.65 0.2 0.35];
plotMatrix(Gamma1_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 1', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.25 0.65 0.2 0.35];
plotMatrix(Gamma2_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 2', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.55 0.65 0.2 0.35];
plotMatrix(Gamma1_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 1', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
ax.OuterPosition = [0.75 0.65 0.2 0.35];
plotMatrix(Gamma2_generic, [-gRng_all gRng_all]);
colormap(cmap_covmat);
beautifygraph('fontscale', 0.833);
%colorbar;
axis equal;
set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
xlabel('odorant 1');
ylabel('odorant 2');
title('Env. 2', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
% ax.OuterPosition = [0 0 1/2 1];
ax.OuterPosition = [0 0 1/2 0.6];
plotDistChange(K1_generic_narrow, K2_generic_narrow, 'beautifyopts', {'fontscale', 0.833});
ylim([0 0.07]);
title('Narrow tuning', 'fontsize', 10);

ax = axes;
ax.Units = 'normalized';
% ax.OuterPosition = [1/2 0 1/2 1];
ax.OuterPosition = [1/2 0 1/2 0.6];
plotDistChange(K1_generic_wide, K2_generic_wide, 'beautifyopts', {'fontscale', 0.833});
ylim([0 0.07]);
title('Wide tuning', 'fontsize', 10);

preparegraph('edge', 0);

% safe_print(fullfile('figs', 'env_change_vs_tuning_different.pdf'));

%% Putting error bars on abundance change, mouse sensing matrix

% random environments, but make them reproducible
rng(335632);

% how many perturbations to calculate
mouse_with_errorbars.n_perturbations = 24;

% setting the unperturbed values of the parameters
% using the mouse sensing matrix
mouse_with_errorbars.S_0 = S_mouse;
% ... with a little bit of noise to remove the degeneracy from the many
% elements of S that are exactly zero
mouse_with_errorbars.S_0(mouse_with_errorbars.S_0 == 0) = ...
    0.1*randn(1, sum(mouse_with_errorbars.S_0(:) == 0));
% mouse_with_errorbars.sigma0 = 1;

% generate random environments
mouse_with_errorbars.Gamma1_0 = generate_environment('rnd_corr', ...
    size(mouse_with_errorbars.S_0, 2));
mouse_with_errorbars.Gamma2_0 = generate_environment('rnd_corr', ...
    size(mouse_with_errorbars.S_0, 2));

% set the unperturbed total number of OSN
% (remember that this is in fact divided by the receptor noise squared)
mouse_with_errorbars.Ktot_0 = 2000;

% initialize storage for the perturbed parameters
cell_empty = cell(1, mouse_with_errorbars.n_perturbations);
mouse_with_errorbars.Gamma1 = cell_empty;
mouse_with_errorbars.Gamma2 = cell_empty;

zero_empty = zeros(mouse_with_errorbars.n_perturbations, 1);
mouse_with_errorbars.Ktot1 = zero_empty;
mouse_with_errorbars.Ktot2 = zero_empty;

mouse_with_errorbars.S1 = cell_empty;
mouse_with_errorbars.S2 = cell_empty;

mouse_with_errorbars.K1 = zeros(size(mouse_with_errorbars.S_0, 1), ...
    mouse_with_errorbars.n_perturbations);
mouse_with_errorbars.K2 = zeros(size(mouse_with_errorbars.S_0, 1), ...
    mouse_with_errorbars.n_perturbations);

mouse_with_errorbars.Q1 = cell_empty;
mouse_with_errorbars.Q2 = cell_empty;

% set perturbation amounts
mouse_with_errorbars.perturbation_sizes.Gamma = 0.01;
mouse_with_errorbars.perturbation_sizes.S = 0.05;
mouse_with_errorbars.perturbation_sizes.Ktot = 0;

% start calculating the perturbations and how they affect the optimal
% receptor distribution
for i = 1:mouse_with_errorbars.n_perturbations
    disp([int2str(i) '/' int2str(mouse_with_errorbars.n_perturbations) '...']);
    
    % perturb environments
    if mouse_with_errorbars.perturbation_sizes.Gamma > 0
        crt_Gamma1 = generate_environment(...
            'delta_rnd_unif', mouse_with_errorbars.Gamma1_0, ...
            'factorsize', mouse_with_errorbars.perturbation_sizes.Gamma);
        crt_Gamma2 = generate_environment(...
            'delta_rnd_unif', mouse_with_errorbars.Gamma2_0, ...
            'factorsize', mouse_with_errorbars.perturbation_sizes.Gamma);
    else
        crt_Gamma1 = mouse_with_errorbars.Gamma1_0;
        crt_Gamma2 = mouse_with_errorbars.Gamma2_0;
    end
    
    mouse_with_errorbars.Gamma1{i} = crt_Gamma1;
    mouse_with_errorbars.Gamma2{i} = crt_Gamma2;
    
    % sanity check: symmetry is preserved
    if max(abs(flatten(crt_Gamma1 - crt_Gamma1'))) > 1e-6
        error('Perturbed Gamma1 not symmetric!');
    end
    if max(abs(flatten(crt_Gamma2 - crt_Gamma2'))) > 1e-6
        error('Perturbed Gamma2 not symmetric!');
    end
    
    % sanity check: positive-definiteness is preserved
    if any(eig(crt_Gamma1) < -1e-6)
        error('Perturbed Gamma1 not positive definite!');
    end
    if any(eig(crt_Gamma2) < -1e-6)
        error('Perturbed Gamma2 not positive definite!');
    end
        
    % perturb Ktot
    crt_Ktot1 = mouse_with_errorbars.Ktot_0*(1 + ...
        mouse_with_errorbars.perturbation_sizes.Ktot*randn);
    crt_Ktot2 = mouse_with_errorbars.Ktot_0*(1 + ...
        mouse_with_errorbars.perturbation_sizes.Ktot*randn);
    
    mouse_with_errorbars.Ktot1(i) = crt_Ktot1;
    mouse_with_errorbars.Ktot2(i) = crt_Ktot2;
    
    % sanity check: Ktots are still positive
    if crt_Ktot1 <= 0
        error('Perturbed Ktot1 <= 0!');
    end
    if crt_Ktot2 <= 0
        error('Perturbed Ktot2 <= 0!');
    end
    
    % perturb S
    crt_S1 = mouse_with_errorbars.S_0 + ...
        mouse_with_errorbars.perturbation_sizes.S *...
        randn(size(mouse_with_errorbars.S_0));
    crt_S2 = mouse_with_errorbars.S_0 + ...
        mouse_with_errorbars.perturbation_sizes.S *...
        randn(size(mouse_with_errorbars.S_0));
    
    mouse_with_errorbars.S1{i} = crt_S1;
    mouse_with_errorbars.S2{i} = crt_S2;
    
    % calculate optimal repertoires
    [crt_K1, ~, crt_Q1] = calculate_optimal_dist(crt_S1, crt_Gamma1, crt_Ktot1, ...
        optim_args{:}, 'lagstart', size(crt_S1, 1)/crt_Ktot1);
    [crt_K2, ~, crt_Q2] = calculate_optimal_dist(crt_S2, crt_Gamma2, crt_Ktot2, ...
        optim_args{:}, 'lagstart', size(crt_S2, 1)/crt_Ktot2);
    
    mouse_with_errorbars.K1(:, i) = crt_K1;
    mouse_with_errorbars.K2(:, i) = crt_K2;
    
    mouse_with_errorbars.Q1{i} = crt_Q1;
    mouse_with_errorbars.Q2{i} = crt_Q2;
    
    % sanity check: total number of neurons is as prescribed
    if abs(sum(crt_K1)/crt_Ktot1 - 1) > 1e-4
        error('K1 doesn''t sum to Ktot1!');
    end
    if abs(sum(crt_K2)/crt_Ktot2 - 1) > 1e-4
        error('K2 doesn''t sum to Ktot2!');
    end
end

%% ... make Ibarra-Soria-like plot

mouse_with_errorbars.summary = struct;

mouse_with_errorbars.summary.K1_mean = mean(mouse_with_errorbars.K1, 2);
mouse_with_errorbars.summary.K1_std = std(mouse_with_errorbars.K1, [], 2);

mouse_with_errorbars.summary.K2_mean = mean(mouse_with_errorbars.K2, 2);
mouse_with_errorbars.summary.K2_std = std(mouse_with_errorbars.K2, [], 2);

logK_ratio = log2(mouse_with_errorbars.summary.K2_mean ./ ...
    mouse_with_errorbars.summary.K1_mean);
logK_std = (mouse_with_errorbars.summary.K1_std ./ mouse_with_errorbars.summary.K1_mean) ...
    + (mouse_with_errorbars.summary.K2_std ./ mouse_with_errorbars.summary.K2_mean);

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2.9 2.3];

errorbar(mouse_with_errorbars.summary.K1_mean, logK_ratio, logK_std, logK_std, ...
    mouse_with_errorbars.summary.K1_std, mouse_with_errorbars.summary.K1_std, ...
    'marker', 'none', 'linestyle', 'none', 'color', [0.7 0.7 0.7], ...
    'capsize', 0);
hold on;
plot(mouse_with_errorbars.summary.K1_mean, logK_ratio, '.', 'markersize', 10);

ylim([-3.5, 3.5]);
set(gca, 'xscale', 'log');
% xlim([5, 350]);
% xlim([1 max(K1p_mean)*1.2]);
xlim([1 100]);
h = plot(xlim, [0 0], 'r');
set(h, 'color', [1 0 0 0.3], 'linewidth', 1);

xlabel('K^{env1}');
ylabel('log_2 K^{env2}/K^{env1}');

beautifygraph('fontscale', 0.833);

set(gca, 'box', 'on', 'xminortick', 'off', 'yminortick', 'off', 'TickLength', [0 0]);

preparegraph('edge', 0);

safe_print(fullfile('figs', 'example_ibarra_soria_style.pdf'));

%% Run tests of population dynamics model

% NOTE: this uses some of the results from the previous section
% (simulations with error bars using mouse sensing matrix)

rng(3488247);

% number of steps to run the simulation for
n_steps = 5000;

% choose a total number of OSN
Ktot_dynmodel = mouse_with_errorbars.Ktot_0;

% calculate optimal receptor distributions, overlap matrices, and mutual
% information functions for one of the samples from the error bars
% simulations
% (K and Q are recalculated; we're mostly just using this to get the mutual
% information function)
[K1_infomax_for_dynmodel, ~, Q1_dynmodel, info1_dynmodel] = ...
    calculate_optimal_dist(mouse_with_errorbars.S1{1}, mouse_with_errorbars.Gamma1{1}, ...
        Ktot_dynmodel, optim_args{:}, ...
        'lagstart', size(mouse_with_errorbars.S1{1}, 1)/Ktot_dynmodel);
[K2_infomax_for_dynmodel, ~, Q2_dynmodel, info2_dynmodel] = ...
    calculate_optimal_dist(mouse_with_errorbars.S2{1}, mouse_with_errorbars.Gamma2{1}, ...
        Ktot_dynmodel, optim_args{:}, ...
        'lagstart', size(mouse_with_errorbars.S2{1}, 1)/Ktot_dynmodel);

% random initial receptor abundances
Kini_dynmodel_random = rand(size(K1_infomax_for_dynmodel));
Kini_dynmodel_random = Kini_dynmodel_random*Ktot_dynmodel/sum(Kini_dynmodel_random);
% run the population dynamics model in the second environment, starting
% from random distribution
% learning rate chosen so that exponential growth has doubling time 2
[~, dyn_history_from_random] = run_dyn_model(Q2_dynmodel, Ktot_dynmodel, ...
    'guessK', Kini_dynmodel_random, 'tolinfo', 0, 'tolK', 0, ...
    'maxsteps', n_steps, 'guesslbd', 0.01679, 'ratelbd', 1e-8, ...
    'rate', log(2)/2);

% initial receptor abundances set to optimum for environment 1, plus random
% jitter
Kini_dynmodel_env1 = abs(K1_infomax_for_dynmodel + ...
    0.05*Ktot_dynmodel/length(K1_infomax_for_dynmodel)*...
    (2*rand(size(K1_infomax_for_dynmodel)) - 1));
% make sure we're still normalized to the correct total number of neurons
Kini_dynmodel_env1 = Kini_dynmodel_env1*Ktot_dynmodel/sum(Kini_dynmodel_env1);
[~, dyn_history_from_env1] = run_dyn_model(Q2_dynmodel, Ktot_dynmodel, ...
    'guessK', Kini_dynmodel_env1, 'tolinfo', 0, 'tolK', 0, ...
    'maxsteps', n_steps, 'guesslbd', 0.01679, 'ratelbd', 1e-8, ...
    'rate', log(2)/2);

% plot the dynamical evolution of receptor abundances
yl = max([max(dyn_history_from_random.K(:)) max(dyn_history_from_env1.K(:))]);
    
col1 = [0.21 0.17 0.53];
col2 = [0.98 0.80 0.17];

fig = figure;
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
plot(n_steps*ones(size(K2_infomax_for_dynmodel)), K2_infomax_for_dynmodel, ...
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

beautifygraph('fontscale', 0.833);

% sanity check: how far are we from maximum info
disp('Random starting point');
disp(['Maximum info found with calculate_optimal_dist: ' num2str(info2_dynmodel(K2_infomax_for_dynmodel), '%.4f') ...
    ' run_dyn_model: ' num2str(info2_dynmodel(dyn_history_from_random.K(:, end)), '%.4f')]);
disp(['(Info at initial point: ' num2str(info2_dynmodel(Kini_dynmodel_random), '%.4f') ')']);
disp(' ');

% make plot starting with distribution optimized for environment 1
ax = axes;
ax.OuterPosition = [0.5 0 0.5 1];
hold on;
% plot optimal distribution in environment 1, as predicted from infomax
% (this is the same as initial conditions, up to the perturbation we added)
plot(ones(size(K1_infomax_for_dynmodel)), K1_infomax_for_dynmodel, ...
    'd', 'markerfacecolor', col1, 'markeredgecolor', 'none', ...
    'markersize', 3);
% plot final initial conditions as predicted from infomax
plot(n_steps*ones(size(K2_infomax_for_dynmodel)), K2_infomax_for_dynmodel, ...
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

beautifygraph('fontscale', 0.833);

% sanity check: how far are we from maximum info
disp('Starting from optimum for different environment');
disp(['Maximum info found with calculate_optimal_dist: ' num2str(info2_dynmodel(K2_infomax_for_dynmodel), '%.4f') ...
    ' run_dyn_model: ' num2str(info2_dynmodel(dyn_history_from_env1.K(:, end)), '%.4f')]);
disp(['(Info at initial point: ' num2str(info2_dynmodel(Kini_dynmodel_env1), '%.4f') ')']);

preparegraph('edge', 0);

% make sure we're saving vector graphics, not rasterized
fig.Renderer = 'painters';

safe_print(fullfile('figs', 'dynamics_example.pdf'));

%% Look at convergence times

% trying out different definition for "convergence time"

% going 90% of the way from initial to final value
t_convergence_90p = zeros(size(K1_infomax_for_dynmodel));
% convergence to within 5% of final value
t_convergence_5p_error = zeros(size(K1_infomax_for_dynmodel));
% convergence to within 5% of Ktot
t_convergence_01p_Ktot = zeros(size(K1_infomax_for_dynmodel));

for i = 1:length(K1_infomax_for_dynmodel)
    crt_trajectory = dyn_history_from_random.K(i, :);
    
    deltaK = crt_trajectory(end) - crt_trajectory(1);
    [~, t_convergence_90p(i)] = min(abs(crt_trajectory - crt_trajectory(1) - 0.9*deltaK));
    
    [~, t_convergence_5p_error(i)] = min(abs(abs(crt_trajectory - crt_trajectory(end)) ...
        - 0.05*crt_trajectory(end)));
    
    [~, t_convergence_01p_Ktot(i)] = min(abs(abs(crt_trajectory - crt_trajectory(end)) ...
        - 0.001*Ktot_dynmodel));
end

%% Try to mimic perturbation from Ibarra-Soria

% random starting environment, but make it reproducible
rng(923967);

Gamma1_ibarra_soria = generate_environment('rnd_corr', size(S_mouse, 2));
details1_ibarra_soria.factor = sqrtm(Gamma1_ibarra_soria);

% [Gamma1_ibarra_soria, details1_ibarra_soria] = generate_environment(...
%     'rnd_product', size(S_mouse, 2), 'factorsize', 0.4/sqrt(size(S_mouse, 2)));

% find odor mixture from Ibarra-Soria
% XXX Eugenol unfortunately is not part of this sensing matrix...
% XXX r-carvone is the same as (-)-carvone
mixture_ibarra_soria = {'Acetophenone', 'Heptanal', '(-)-Carvone'};
idx_odors_ibarra_soria = cellfun(@(name) find(strcmp(mouse_responses.odorants, name)), ...
    mixture_ibarra_soria);
n_odors_ibarra_soria = length(idx_odors_ibarra_soria);

% we're adding a large amount of variation to the Ibarra-Soria odorants
ibarra_soria_amount = 0.5;
sqrtGamma2_ibarra_soria = details1_ibarra_soria.factor;
sqrtGamma2_ibarra_soria(:, idx_odors_ibarra_soria) = ...
    sqrtGamma2_ibarra_soria(:, idx_odors_ibarra_soria) + ...
    ibarra_soria_amount*randn(size(sqrtGamma2_ibarra_soria, 1), n_odors_ibarra_soria);

Gamma2_ibarra_soria = sqrtGamma2_ibarra_soria'*sqrtGamma2_ibarra_soria;

% visualize the environments

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 12 6];
subplot(1, 2, 1);
g1Rng = max(abs(Gamma1_ibarra_soria(:)));
g2Rng = max(abs(Gamma2_ibarra_soria(:)));
% gRng = max(g1Rng, g2Rng);
gRng = g1Rng;
plotMatrix(Gamma1_ibarra_soria, [-gRng gRng]);
colormap(cmap_covmat);
%colorbar;
xlabel('odorant 1');
ylabel('odorant 2');
title('Control environment');

subplot(1, 2, 2);
plotMatrix(Gamma2_ibarra_soria, [-gRng gRng]);
colormap(cmap_covmat);
%colorbar;
xlabel('odorant 1');
ylabel('odorant 2');
title('Perturbed environment');

preparegraph;

%% ... run the simulation

Ktot_ibarra_soria = 15000;

% using same value for Ktot
[K1_ibarra_soria, ~, Q1_ibarra_soria] = calculate_optimal_dist(S_mouse/sigma_mouse, ...
    Gamma1_ibarra_soria, Ktot_ibarra_soria, optim_args{:});
[K2_ibarra_soria, ~, Q2_ibarra_soria] = calculate_optimal_dist(S_mouse/sigma_mouse, ...
    Gamma2_ibarra_soria, Ktot_ibarra_soria, optim_args{:});

%% ... plot the difference

logK_ratio_ibarra_soria = log2(K2_ibarra_soria ./ K1_ibarra_soria);

figure;
% fig.Units = 'inches';
% fig.Position = [fig.Position(1:2) 2.9 2.3];

hold on;
plot(K1_ibarra_soria, logK_ratio_ibarra_soria, '.', 'markersize', 10);

ylim([-3.5, 3.5]);
set(gca, 'xscale', 'log');
xlim([1 700]);
h = plot(xlim, [0 0], 'r');
set(h, 'color', [1 0 0 0.3], 'linewidth', 1);

xlabel('K^{env1}');
ylabel('log_2 K^{env2}/K^{env1}');

beautifygraph('fontscale', 0.833);

set(gca, 'box', 'on', 'xminortick', 'off', 'yminortick', 'off', 'TickLength', [0 0]);

preparegraph('edge', 0);

%% Sensitivity of abundance change to starting environment

% random starting environment, but make it reproducible
rng(9834578);

start_sensitivity.n_samples = 30;
start_sensitivity.S = S_mouse;
start_sensitivity.sigma = sigma_mouse;

start_sensitivity.Ktot = 50000;

start_sensitivity.beta_corr = [4 15 50];
% start_sensitivity.beta_corr = 50;

start_sensitivity.diff_amount = 0.5;
start_sensitivity.idxs = randi(size(S_mouse, 2), 5, 1);

start_sensitivity.Kdiff = zeros(size(S_mouse, 1), start_sensitivity.n_samples);

start_sensitivity.Gamma1 = cell(1, start_sensitivity.n_samples);
start_sensitivity.Gamma2 = cell(1, start_sensitivity.n_samples);

start_sensitivity.K1 = zeros(size(S_mouse, 1), start_sensitivity.n_samples);
start_sensitivity.K2 = zeros(size(S_mouse, 1), start_sensitivity.n_samples);

start_sensitivity.Q1 = cell(1, start_sensitivity.n_samples);
start_sensitivity.Q2 = cell(1, start_sensitivity.n_samples);

start_sensitivity.all_beta_corr = zeros(1, start_sensitivity.n_samples);

for i = 1:start_sensitivity.n_samples
    disp(i);
    
    crt_beta_corr = start_sensitivity.beta_corr(floor(...
            (i-1)/start_sensitivity.n_samples*length(start_sensitivity.beta_corr)) ...
        + 1);
    start_sensitivity.all_beta_corr(i) = crt_beta_corr;
    
    start_sensitivity.Gamma1{i} = generate_environment('rnd_corr', size(S_mouse, 2), ...
        'corrbeta', crt_beta_corr);
    
    % we're adding a large amount of variation to the Ibarra-Soria odorants
    sqrtGamma2_start_sensitivity = sqrtm(start_sensitivity.Gamma1{i});
    sqrtGamma2_start_sensitivity(:, start_sensitivity.idxs) = ...
        sqrtGamma2_start_sensitivity(:, start_sensitivity.idxs) + ...
        start_sensitivity.diff_amount*...
        randn(size(sqrtGamma2_start_sensitivity, 1), length(start_sensitivity.idxs));
    start_sensitivity.Gamma2{i} = sqrtGamma2_start_sensitivity'*sqrtGamma2_start_sensitivity;
    
    [start_sensitivity.K1(:, i), ~, start_sensitivity.Q1{i}] = ...
        calculate_optimal_dist(start_sensitivity.S/start_sensitivity.sigma, ...
        start_sensitivity.Gamma1{i}, ...
        start_sensitivity.Ktot, optim_args{:});
    [start_sensitivity.K2(:, i), ~, start_sensitivity.Q2{i}] = ...
        calculate_optimal_dist(start_sensitivity.S/start_sensitivity.sigma, ...
        start_sensitivity.Gamma2{i}, ...
        start_sensitivity.Ktot, optim_args{:});
    
    % sanity check: total number of neurons is as prescribed
    if abs(sum(start_sensitivity.K1(:, i))/start_sensitivity.Ktot - 1) > 1e-4
        error('K1 doesn''t sum to Ktot!');
    end
    if abs(sum(start_sensitivity.K2(:, i))/start_sensitivity.Ktot - 1) > 1e-4
        error('K2 doesn''t sum to Ktot!');
    end
    
    start_sensitivity.Kdiff(:, i) = start_sensitivity.K2(:, i) - ...
        start_sensitivity.K1(:, i);
end

%% Calculate some summaries

start_sensitivity.summary.Kdiff_mean = mean(start_sensitivity.Kdiff, 2);
start_sensitivity.summary.Kdiff_std = std(start_sensitivity.Kdiff, [], 2);






%% SCRATCH

%% Histogram of deltaK in wide vs. narrow tuning regimes

rng(2334);

n_samples = 100;

K1_all_narrow = zeros(size(S_fly, 1), n_samples);
K2_all_narrow = zeros(size(S_fly, 1), n_samples);

K1_all_wide = zeros(size(S_fly, 1), n_samples);
K2_all_wide = zeros(size(S_fly, 1), n_samples);

for i = 1:n_samples
    disp(i);
    
    crt_Gamma1 = generate_environment('rnd_corr', size(S_fly, 2));
    crt_Gamma2 = generate_environment('rnd_corr', size(S_fly, 2));

    % calculate the optimal distributions at both narrow and wide tuning, in
    % both environments
    K1_all_narrow(:, i) = calculate_optimal_dist(S_narrow/sigma_narrow, ...
        crt_Gamma1, Ktot_generic, optim_args{:}, 'lagstart', size(S_fly, 1)/Ktot_generic);
    K2_all_narrow(:, i) = calculate_optimal_dist(S_narrow/sigma_narrow, ...
        crt_Gamma2, Ktot_generic, optim_args{:}, 'lagstart', size(S_fly, 1)/Ktot_generic);

    K1_all_wide(:, i) = calculate_optimal_dist(S_wide/sigma_wide, crt_Gamma1, ...
        Ktot_generic, optim_args{:}, 'lagstart', size(S_fly, 1)/Ktot_generic);
    K2_all_wide(:, i) = calculate_optimal_dist(S_wide/sigma_wide, crt_Gamma2, ...
        Ktot_generic, optim_args{:}, 'lagstart', size(S_fly, 1)/Ktot_generic);

    % sanity check on convergence: make sure the total number of neurons is as
    % prescribed
    if abs(sum(K1_all_narrow(:, i)) / Ktot_generic - 1) > 1e-4
        error('K1_narrow doesn''t sum up to Ktot_generic.');
    end
    if abs(sum(K2_all_narrow(:, i)) / Ktot_generic - 1) > 1e-4
        error('K2_narrow doesn''t sum up to Ktot_generic.');
    end
    if abs(sum(K1_all_wide(:, i)) / Ktot_generic - 1) > 1e-4
        error('K1_wide doesn''t sum up to Ktot_generic.');
    end
    if abs(sum(K2_all_wide(:, i)) / Ktot_generic - 1) > 1e-4
        error('K2_wide doesn''t sum up to Ktot_generic.');
    end
end