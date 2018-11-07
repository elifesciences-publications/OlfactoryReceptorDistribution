% Make plots meant to show generic properties of the optimization.

%% Load Hallem&Carlson sensing data

sensing_fly = open('data/flyResponsesWithNames.mat');
S_fly = sensing_fly.relRates';

%% Plot a sample of the Hallem&Carlson sensing data

odor_choices = {'1-hexanol', 'E2-hexenol', 'E3-hexenol', 'Z3-hexenol', ...
    '2-heptanone', 'butyl acetate', 'pentyl acetate', 'methyl hexanoate', ...
    'ethyl hexanoate', 'methyl benzoate', 'ethyl benzoate'};
or_choices = {'Or2a', 'Or19a', 'Or35a', 'Or47b', 'Or67a', 'Or85b'};

odor_idxs = cellfun(@(s) find(strcmp(sensing_fly.odorNames, s)), odor_choices);
or_idxs = cellfun(@(s) find(strcmp(sensing_fly.orNames, s)), or_choices);

fig = figure;
fig.Position(3:4) = fig.Position(3)*[1 length(or_choices)/length(odor_choices)];

ax = axes;
ax.Position = [0 0 1 1];

% cmap = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64);
cmap = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

cl = quantile(flatten(S_fly(or_idxs, odor_idxs)), 0.9);

plotheat(S_fly(or_idxs, odor_idxs), [-cl cl]);
colormap(cmap);

axis off;

preparegraph;

safeprint(fullfile('figs', 'example_sensing'));

%% Make environment

% random environment, but make it reproducible
rng(92376);

% * the concentrations used in Hallem&Carlson were high; here we're using
%   concentrations 1000 times smaller (so that the variance is about 1e-6)
% * the actual variances of odor concentrations in natural scenes are
%   unknown, so we effectively fit them to obtain neuron numbers in the
%   right ballpark
Gamma = generate_environment('rnd_diag_const', size(S_fly, 2), ...
    'diag_mu', 1e-6, 'diag_size', 1e-6, 'offdiag_mu', 1e-7);

%% Get distribution at various values for Ktot

% Ktot_values = [0 logspace(log10(0.1), log10(15000), 100)];
Ktot_values = [0 logspace(log10(5), log10(7.5e5), 100)];

% normalize by standard deviations of background rates
S_fly_normalized = bsxfun(@rdivide, S_fly, sensing_fly.bkgStd');

[K, info_values, Q, info_fct] = calculate_optimal_dist(...
    S_fly_normalized, Gamma, Ktot_values);

[M, N] = size(S_fly);

% find where each receptor kicks in
minK_by_type = ones(M, 1);
K_nz = bsxfun(@rdivide, K, Ktot_values) > 0.004;
for i = 1:M
    idx = find(~K_nz(i, :), 1, 'last');
    if ~isempty(idx)
        minK_by_type(i) = idx+1;
    end
end
[minK_ordered, rec_order] = sort(minK_by_type);

%% Plot "natural" receptor usage

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 6.0 1.8];

% make the background white so the plot looks the same in the window as it
% looks on paper
fig.Color = [1 1 1];

% plot the dependence on Ktot for small Ktot
ax1 = axes;
ax1.Units = 'normalized';
ax1.OuterPosition = [1/3 0 1/3 1];

% approximate maximum value of Ktot for this plot
smallKtot_limit = 130;
[~, idx_max] = min(abs(Ktot_values - smallKtot_limit));

% show the OSN numbers as an area plot
h_area = area(Ktot_values(1:idx_max), K(:, 1:idx_max)', ...
    'linewidth', 0.5, 'facecolor', 'flat');
hold on;

% mark transitions where new receptors kick in
for i = 1:M-1
    if minK_ordered(i+1) <= length(Ktot_values) && minK_ordered(i+1) <= idx_max
        line(repmat(Ktot_values(minK_ordered(i+1)), [2 1]), ylim, ...
            'linestyle', ':', 'color', [0.7 0.7 0.7], 'linewidth', 1);
    end
end

xlim([0 Ktot_values(idx_max)]);
ylim([0 Ktot_values(idx_max)]);

xlh = xlabel('Number of OSNs');
xlh.Units = 'characters';
xlh.Position(2) = xlh.Position(2) - 0.15;
ylabel('OSNs by receptor type');

% find out which receptors we've displayed, and in what order they appeared
% so we can create a legend
rec_dispd = rec_order(minK_ordered <= idx_max);
h_leg = legend(h_area(rec_dispd), sensing_fly.orNames(rec_dispd), 'location', 'northwest', ...
    'fontsize', 6);
h_leg.Position = [h_leg.Position(1) - 0.005 h_leg.Position(2) + 0.09 h_leg.Position(3:4)];
legend('boxoff');
beautifygraph(ax1, 'fontscale', 0.667, 'ticksize', 12);

ax1.XMinorTick = 'off';
ax1.YMinorTick = 'off';
ax1.TickLength = 3*ax1.TickLength;

ax1.LineWidth = 0.5;

% plot the dependence on Ktot for large Ktot
ax2 = axes;
ax2.Units = 'normalized';
ax2.OuterPosition = [0 0 1/3 1];

% show the OSN numbers as an area plot
area(Ktot_values, K', 'linewidth', 0.5, 'facecolor', 'flat');
hold on;

xlim([0 Ktot_values(end)]);
ylim([0 Ktot_values(end)]);

ylabel('OSNs by receptor type');
beautifygraph('fontscale', 0.667, 'ticksize', 12);

% format the ticks in a nicer way
ax2.XTick = [0 4*10^5];
ax2.XTickLabel = {0, '4 \times 10^5'};

xlh = xlabel('Number of OSNs');
xlh.Units = 'characters';
xlh.Position(2) = xlh.Position(2) - 0.13;
beautifygraph(ax2, 'fontscale', 0.667);

ax2.XMinorTick = 'off';
ax2.YMinorTick = 'off';
ax2.TickLength = 3*ax2.TickLength;

ax2.LineWidth = 0.5;

% make plot showing how the number of receptor expressed in non-zero
% numbers of neurons grows with the total number of neurons
ax3 = axes;
ax3.Units = 'normalized';
ax3.OuterPosition = [2/3 0 1/3 1];
rectype_count = zeros(size(Ktot_values));
for i = 1:length(minK_ordered)
    rectype_count(minK_ordered(i):end) = rectype_count(minK_ordered(i):end) + 1;
end
semilogx(Ktot_values, rectype_count, 'linewidth', 1);
xlabel('Number of OSNs');
ylabel('Active receptor types');
xlim([Ktot_values(2), Ktot_values(end)]);
ylim([0 M+1]);
beautifygraph(ax3, 'fontscale', 0.667, 'ticksize', 12);

% select a good crop of ticks
ax3.XTick = [10^1 10^2 10^3 10^4 10^5];

ax3.XMinorTick = 'off';
ax3.YMinorTick = 'off';
ax3.TickLength = 3*ax3.TickLength;

ax3.LineWidth = 0.5;

preparegraph('edge', 0);

% try to ensure that all the plots are aligned properly
ymax = max([ax1.Position(2) ax2.Position(2) ax3.Position(2)]);
hmin = min([ax1.Position(4) ax2.Position(4) ax3.Position(4)]);
ax1.Position(2) = ymax;
ax2.Position(2) = ymax;
ax3.Position(2) = ymax;
ax1.Position(4) = hmin;
ax2.Position(4) = hmin;
ax3.Position(4) = hmin;

safeprint(fullfile('figs', 'dependence_on_Ktot'));
