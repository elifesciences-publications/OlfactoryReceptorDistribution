% Make plots showing how the tuning of receptors and the kind of environment
% changes affect the optimal receptor distribution.

%% Load the results

load(fullfile('save', 'change_nonoverlapping.mat'));

%% Define a nice colormap

cmap_covmat = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

%% Plot example covariance matrices

env_matrices = {...
    {'env_generic1', Gamma1_generic{1}, []}, ...
    {'env_generic2', Gamma2_generic{1}, []}, ...
    {'env_noverlap1', Gamma1_nonoverlapping{1}, masks{1}}, ...
    {'env_noverlap2', Gamma2_nonoverlapping{1}, masks{1}}};

% figure out a good color range for the matrices
pooled_entries = flatten(cell2mat(cellfun(@(c) c{2}, env_matrices, 'uniform', false)));

% env_color_lim = 0.435;
env_color_lim = quantile(abs(pooled_entries), 0.95);

for i = 1:length(env_matrices)
    crt_mat = env_matrices{i}{2};
    
    % order rows and columns of matrix to make obvious the block structure
    crt_mask = env_matrices{i}{3};
    if ~isempty(crt_mask)
        idxs1 = find(crt_mask);
        idxs2 = find(~crt_mask);
        idxs = [idxs1(:) idxs2(:)];
        crt_mat = crt_mat(idxs, idxs);
    end
    
    % draw the matrix
    fig = figure;
    
    ax = axes;
    ax.Units = 'pixels';
    ax.Position = [ax.Position(1:2) fliplr(size(crt_mat))];
    
    plotheat(crt_mat, [-env_color_lim env_color_lim]);
    colormap(cmap_covmat);
    
    beautifygraph('fontscale', 0.667);
    preparegraph;

    axis equal;
    set(ax, 'xtick', [], 'ytick', [], 'box', 'on');
        
    fig.Units = 'pixels';
    fig.Position = [fig.Position(1:2) ax.OuterPosition(3:4)];
    
    crt_frame = getframe;
    
     imwrite(crt_frame.cdata, fullfile('figs', [env_matrices{i}{1} '.png']));
end

%% Plot distribution of \Delta \Gamma_{ij}

fig = figure;
fig.Color = [1 1 1];

fig.Units = 'inches';
fig.Position(3:4) = [3 2];

% calculate change in entries of covariance matrices
diff_generic = arrayfun(@(i) Gamma1_generic{i}(:) - Gamma2_generic{i}(:), ...
    1:n_env_samples, 'uniform', false);
diff_nonoverlapping = arrayfun(@(i) Gamma1_nonoverlapping{i}(:) - Gamma2_nonoverlapping{i}(:), ...
    1:n_env_samples, 'uniform', false);

% find a reasonable range for the histogram plots
pooled_range = quantile(flatten(cell2mat(diff_generic)), [0.025 0.975]);

cov_bins = linspace(pooled_range(1), pooled_range(2), 100);

% draw traces from all the samples
hold on;
for i = 1:n_env_samples
    histogram(diff_generic{i}(:), cov_bins, ...
        'edgecolor', cmap_covmat(end, :), 'edgealpha', 0.02, 'linewidth', 0.5, ...
        'normalization', 'pdf', 'displaystyle', 'stairs');
    histogram(diff_nonoverlapping{i}(:), cov_bins, ...
        'edgecolor', cmap_covmat(1, :), 'edgealpha', 0.02, 'linewidth', 0.5, ...
        'normalization', 'pdf', 'displaystyle', 'stairs');
end

% draw the pooled histogram with a thicker line
pooled_diff_generic = flatten(cell2mat(diff_generic));
pooled_diff_nonoverlapping = flatten(cell2mat(diff_nonoverlapping));
h1 = histogram(pooled_diff_generic, cov_bins, ...
    'edgecolor', cmap_covmat(end, :), 'linewidth', 1, ...
    'normalization', 'pdf', 'displaystyle', 'stairs');
h2 = histogram(pooled_diff_nonoverlapping, cov_bins, ...
    'edgecolor', cmap_covmat(1, :), 'linewidth', 1, ...
    'normalization', 'pdf', 'displaystyle', 'stairs');

% fix the labels and ranges
xlabel('\Delta\Gamma_{ij}');
ylabel('PDF');

xlim(pooled_range);

% add a legend
hl = legend([h1 h2], {'generic', 'non-overlapping'}, 'fontsize', 6, ...
    'location', 'northeast');
hl.Position(1:2) = hl.Position(1:2) + [0.1 0.1];
legend('boxoff');

% beautify, prepare for printing, and print
beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'ticksize', 12);

preparegraph;

safeprint(fullfile('figs', 'deltacov_dist_generic_cs_nonoverlapping'), ...
    'type', 'png');

%% Plot \Delta K distributions

fig = figure;
fig.Color = [1 1 1];

fig.Units = 'inches';
fig.Position(3:4) = [3 2];

% calculate change in optimal receptor distributions
diff_K_generic = arrayfun(@(i) K1_generic{i}(:) - K2_generic{i}(:), ...
    1:n_env_samples, 'uniform', false);
diff_K_nonoverlapping = arrayfun(@(i) K1_nonoverlapping{i}(:) - K2_nonoverlapping{i}(:), ...
    1:n_env_samples, 'uniform', false);

% pool all the changes in receptor numbers together
% use absolute differences since by definition the mean change is 0
% (because Ktot is fixed)
pooled_diff_K_generic = abs(flatten(cell2mat(diff_K_generic)));
pooled_diff_K_nonoverlapping = abs(flatten(cell2mat(diff_K_nonoverlapping)));

% find a suitable range for the histograms
pooled_K_range = quantile(pooled_diff_K_generic, 0.95);
% pooled_K_range = max([pooled_diff_K_generic(:) ; pooled_diff_K_nonoverlapping(:)]);

K_bins = linspace(0, pooled_K_range, 100);

% draw the histograms
hold on;
% h1 = histogram(pooled_diff_K_generic(:), K_bins, ...
%     'edgecolor', cmap_covmat(end, :), 'linewidth', 1, ...
%     'normalization', 'cdf', 'displaystyle', 'stairs');
% h2 = histogram(pooled_diff_K_nonoverlapping(:), K_bins, ...
%     'edgecolor', cmap_covmat(1, :), 'linewidth', 1, ...
%     'normalization', 'cdf', 'displaystyle', 'stairs');
% n1 = histcounts(pooled_diff_K_generic(:), K_bins, 'normalization', 'cdf');
% n2 = histcounts(pooled_diff_K_nonoverlapping(:), K_bins, 'normalization', 'cdf');
n1 = histcounts(pooled_diff_K_generic(:), K_bins, 'normalization', 'pdf');
n2 = histcounts(pooled_diff_K_nonoverlapping(:), K_bins, 'normalization', 'pdf');

h1 = plot(K_bins(1:end-1), n1, 'color', cmap_covmat(end, :), 'linewidth', 1);
h2 = plot(K_bins(1:end-1), n2, 'color', cmap_covmat(1, :), 'linewidth', 1);

threshold = 50;
x_shift = 3;
plot([threshold threshold], [0 0.03], 'color', [0.5 0.5 0.5], 'linewidth', 0.5);

generic_y = 0.028;
text(threshold - 0.5*x_shift, generic_y, [int2str(100*mean(pooled_diff_K_generic <= threshold)) '%<'], ...
    'color', cmap_covmat(end, :), 'horizontalalignment', 'right');
text(threshold + x_shift, generic_y, ['>' int2str(100*mean(pooled_diff_K_generic > threshold)) '%'], ...
    'color', cmap_covmat(end, :), 'horizontalalignment', 'left');

nonovl_y = 0.022;
text(threshold - 0.5*x_shift, nonovl_y, [int2str(100*mean(pooled_diff_K_nonoverlapping <= threshold)) '%<'], ...
    'color', cmap_covmat(1, :), 'horizontalalignment', 'right');
text(threshold + x_shift, nonovl_y, ['>' int2str(100*mean(pooled_diff_K_nonoverlapping > threshold)) '%'], ...
    'color', cmap_covmat(1, :), 'horizontalalignment', 'left');

% fix axes labels an ranges
xlabel('|\DeltaK_i|');
% ylabel('CDF');
ylabel('PDF');

% xlim([0 pooled_K_range]);
% ylim([0 1]);
xlim([0 200]);
ylim([0 0.03]);

% add a legend
hl = legend([h1 h2], {'generic', 'non-overlapping'}, 'location', 'southeast');
hl.Position(1:2) = hl.Position(1:2) + [0.05 0.1];
legend('boxoff');

% beautify, prepare for printing, and print
beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'ticksize', 12);

preparegraph;

% calculate p-values for differences between distributions
[~, ks_p] = kstest2(pooled_diff_K_generic, pooled_diff_K_nonoverlapping);
% rank-sum test is one-tailed with the alternative hypothesis being that
% the \Delta Ks tend to be larger in the non-overlapping regime
rs_p = ranksum(pooled_diff_K_generic, pooled_diff_K_nonoverlapping, ...
    'tail', 'left');

disp(['KS test p = ' num2str(ks_p, '%g') '.']);
disp(['Wilcoxon rank sum test p = ' num2str(rs_p, '%g') '.']);

% safeprint(fullfile('figs', 'deltaK_cdf_generic_vs_nonoverlapping'));
safeprint(fullfile('figs', 'deltaK_pdf_generic_vs_nonoverlapping'));

%% Load tuning results

load(fullfile('save', 'change_vs_tuning.mat'));

%% Plot \Delta K vs. tuning per receptor

fig = figure;
fig.Color = [1 1 1];

fig.Units = 'inches';
fig.Position(3:4) = [3 2];

diffK_varied = flatten(cell2mat(arrayfun(@(i) K2_varied{i} - K1_varied{i}, ...
    1:n_sensing_samples, 'uniform', false)));
sigmas_varied_flat = flatten(cell2mat(sigmas_varied));
S_varied_ipr = flatten(cell2mat(cellfun(@(v) ipr(abs(v)), S_varied, 'uniform', false)));

% bin_edges = linspace(min(S_varied_ipr), max(S_varied_ipr), 5);
bin_edges = linspace(min(sigmas_varied_flat), max(sigmas_varied_flat) + eps, 4);

ipr_step = diff(bin_edges(1:2));
% discretized_ipr = bin_edges(1) + ...
%     floor((S_varied_ipr - bin_edges(1)) / ipr_step) * ipr_step;
discretized_sigmas = bin_edges(1) + ...
    floor((sigmas_varied_flat - bin_edges(1)) / ipr_step) * ipr_step;

abs_diffK_varied = abs(diffK_varied);
pooled_diffK = arrayfun(@(i) median(abs_diffK_varied(S_varied_ipr >= bin_edges(i) & ...
    S_varied_ipr < bin_edges(i+1))), 1:length(bin_edges)-1);
pooled_diffK_std = arrayfun(@(i) std(abs_diffK_varied(S_varied_ipr >= bin_edges(i) & ...
    S_varied_ipr < bin_edges(i+1))), 1:length(bin_edges)-1);

% mask = abs(diffK_varied) > 1e-3;
% stripplot(discretized_ipr, diffK_varied, 'jitter', ipr_step*0.5, 'kde', true, ...
%     'boxes', true);
hold on;
% zl = plot([min(discretized_sigmas) max(discretized_sigmas)], [0, 0], ':k');
plot([0.1 0.7], [0 0], ':k', 'linewidth', 0.5);
stripplot(discretized_sigmas, diffK_varied, 'jitter', ipr_step*0.5, 'kde', true, ...
    'boxes', true, 'colors', [0.788, 0.365, 0.388], 'marker', '.', ...
    'boxOpts', {'linewidth', 1, 'color', 'k'});
% stripplot(discretized_ipr(mask), diffK_varied(mask), 'jitter', ipr_step*0.5, 'kde', true, ...
%     'boxes', true);

% hold on;
% bin_centers = 0.5*(bin_edges(1:end-1) + bin_edges(2:end));
% bar(bin_centers, pooled_diffK);
% errorbar(bin_centers, pooled_diffK, pooled_diffK_std);

% scatterfit(S_varied_ipr, abs(diffK_varied));

yl = max(abs(ylim));
ylim(1.1*yl*[-1, 1]);

xlabel('Receptor tuning width');
ylabel('\DeltaK');

set(gca, 'xtick', [0.2, 0.4, 0.6], 'xticklabel', {...
    '[0.2, 0.4]', '[0.4, 0.6]', '[0.6, 0.8]'});

beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'ticksize', 12, ...
    'minorticks', 'off');
preparegraph;

safeprint(fullfile('figs', 'deltaK_vs_tuning'));

%% SCRATCH BELOW

%% Plot distribution of log relative change in K

fig = figure;
fig.Color = [1 1 1];

fig.Units = 'inches';
fig.Position(3:4) = [3 2];

% calculate change in optimal receptor distributions
logratio_K_generic = arrayfun(@(i) log2(K2_generic{i}(:) ./ K1_generic{i}(:)), ...
    1:n_env_samples, 'uniform', false);
logratio_K_nonoverlapping = arrayfun(@(i) log2(K2_nonoverlapping{i}(:) ./ K1_nonoverlapping{i}(:)), ...
    1:n_env_samples, 'uniform', false);

% pool all the changes in receptor numbers together
% use absolute differences since by definition the mean change is 0
% (because Ktot is fixed)
pooled_logratio_K_generic = abs(flatten(cell2mat(logratio_K_generic)));
pooled_logratio_K_nonoverlapping = abs(flatten(cell2mat(logratio_K_nonoverlapping)));

% find a suitable range for the histograms
% pooled_K_logratio_range = quantile(pooled_logratio_K_generic, [0.025 0.975]);
pooled_K_logratio_range = [0 0.5];

K_logratio_bins = linspace(pooled_K_logratio_range(1), pooled_K_logratio_range(2), 100);

% draw the histograms
hold on;
h1 = histogram(pooled_logratio_K_generic(:), K_logratio_bins, ...
    'edgecolor', cmap_covmat(end, :), 'linewidth', 1, ...
    'normalization', 'pdf', 'displaystyle', 'stairs');
h2 = histogram(pooled_logratio_K_nonoverlapping(:), K_logratio_bins, ...
    'edgecolor', cmap_covmat(1, :), 'linewidth', 1, ...
    'normalization', 'pdf', 'displaystyle', 'stairs');

% fix axes labels an ranges
xlabel('|log_2 K_i^{(2)} / K_i^{(1)}|');
ylabel('PDF');

xlim([0 pooled_K_logratio_range(2)]);

% add a legend
hl = legend([h1 h2], {'generic', 'non-overlapping'});
hl.Position(1:2) = hl.Position(1:2) + [0.1 0.1];
legend('boxoff');

% beautify, prepare for printing, and print
beautifygraph('fontscale', 0.667, 'linewidth', 0.5, 'ticksize', 12);

preparegraph;

% calculate p-values for differences between distributions
[~, ks_p] = kstest2(pooled_logratio_K_generic, pooled_logratio_K_nonoverlapping);
% rank-sum test is one-tailed with the alternative hypothesis being that
% the \Delta Ks tend to be larger in the non-overlapping regime
rs_p = ranksum(pooled_logratio_K_generic, pooled_logratio_K_nonoverlapping, ...
    'tail', 'left');

disp(['KS test p = ' num2str(ks_p, '%g') '.']);
disp(['Wilcoxon rank sum test p = ' num2str(rs_p, '%g') '.']);
