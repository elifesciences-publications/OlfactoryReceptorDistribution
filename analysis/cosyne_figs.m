% Making some figures for Cosyne 2017 abstract.

%% Example of receptor usage

% make an artificial sensing matrix and an artificial odor environment
rng('default');

M = 4;
N = 10;
S = 0.1*randn(M, N);

sigma = 10.0;
Gamma = 0.001*ones(N) + eye(N);

% calculate the receptor distribution
Ktot_values = logspace(log10(1), log10(250), 100);
[K, info_values, ~, ~] = calculate_optimal_dist(S/sigma, Gamma, Ktot_values);

% find where each receptor kicks in
minK_by_type = ones(M, 1);
K_nz = bsxfun(@rdivide, K, Ktot_values) > 0.004;
for i = 1:M
    idx = find(~K_nz(i, :), 1, 'last');
    if ~isempty(idx)
        minK_by_type(i) = idx+1;
    end
end
% calculate information values with subsets of receptors
% 1st column: only using the one with smallest minK
% 2nd column: using the two receptors with smallest minK
% etc.
[minK_ordered, rec_order] = sort(minK_by_type);
info_values_sub = repmat(info_values(:), [1 M-1]);
for i = 1:M-1
    subS = S(rec_order(1:i), :);
    [~, info_values_sub(:, i)] = calculate_optimal_dist(subS/sigma, Gamma, Ktot_values);
end

%% Plot example of receptor usage

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 6 2.8];

% plot the receptor distribution
subplot(1, 2, 1);
area(Ktot_values, K', 'linewidth', 1);
hold on;
for i = 1:M-1
    if minK_ordered(i+1) <= length(Ktot_values)
        line(repmat(Ktot_values(minK_ordered(i+1)), [2 1]), ylim, ...
            'linestyle', ':', 'color', [0.7 0.7 0.7], 'linewidth', 1);
    end
end
axis equal;
xlim([0 max(Ktot_values)]);
ylim([0 max(Ktot_values)]);

xlabel('Total OSN number');
ylabel('Number of OSN by receptor type');
legend({'OR 1', 'OR 2', 'OR 3', 'OR 4'}, 'location', 'northwest');
beautifygraph;

ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;

% plot the information per OSN
subplot(1, 2, 2);
hold on;
%color_choice = parula(4);
for i = 1:M-1
    plot(Ktot_values, info_values_sub(:, i)./Ktot_values(:), '--', ...
        'color', [0.4 0.4 0.4]);
end
plot(Ktot_values, info_values./Ktot_values, 'k', 'linewidth', 2);
axis square;
xlabel('Total OSN number');
ylabel('Information per OSN (nats)');
xlim([0 max(Ktot_values)]);
tmp = ylim;
for i = 1:M-1
    if minK_ordered(i+1) <= length(Ktot_values)
        line(repmat(Ktot_values(minK_ordered(i+1)), [2 1]), tmp, ...
            'linestyle', ':', 'color', [0.7 0.7 0.7], 'linewidth', 1);
    end
end
ylim(tmp);

beautifygraph;

ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;

preparegraph('edge', 0);

safe_print(fullfile('figs', 'example_receptor_dist.pdf'));

%% Effect of changing environment

% for reproducibility
rng('default');

M = 4;  % # receptor types
N = 10; % # odorants

% make a response matrix in which odorants 1 and 2 only excite the first
% and second receptor, respectively
% the first and second receptor otherwise have identical responses to all
% other odorants
% other receptors responds randomly to all other odorants, but not to 1 & 2
bkg_response = 1;
special_response = 1;
S = bkg_response*[repmat(rand(1, N), 2, 1) ; rand(M-2, N)];

S(:, 1) = 0;
S(:, 2) = 0;
S(1, 1) = special_response;
S(2, 2) = special_response;

sigma = 10;

Gamma = 0.1*diag(1 + rand(1, N));
Gamma1 = Gamma + diag([1 0 zeros(1, N-2)]);
Gamma2 = Gamma + diag([0 1 zeros(1, N-2)]);

Ktot_values = logspace(log10(1), log10(10000), 50);
[K1, info_values_1, Q1, ~] = calculate_optimal_dist(S/sigma, Gamma1, Ktot_values); %#ok<ASGLU>
[K2, info_values_2, Q2, ~] = calculate_optimal_dist(S/sigma, Gamma2, Ktot_values); %#ok<ASGLU>

%% Make plot for changing environment

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 6 2];

lowSnrIdx = 36;
highSnrIdx = 50;

% show change at low SNR
subplot(1, 2, 1);
bar([K1(:, lowSnrIdx) K2(:, lowSnrIdx)]/Ktot_values(lowSnrIdx));
ylabel('Fraction of OSN');
legend({'env 1', 'env 2'});
legend('boxoff');
title('Low SNR');
xlim([0 4.5]);
ylim([0 0.6]);
beautifygraph;
ax = gca;
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;
ax.XTick = 1:M;
ax.XTickLabel = arrayfun(@(i) ['OR ' int2str(i)], 1:M, 'uniform', false);

subplot(1, 2, 2);
bar([K1(:, highSnrIdx) K2(:, highSnrIdx)]/Ktot_values(highSnrIdx));
xlim([0 4.5]);
ylim([0 0.6]);
ylabel('Fraction of OSN');
title('High SNR');
legend({'env 1', 'env 2'});
legend('boxoff');
beautifygraph;
ax = gca;
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;
ax.XTick = 1:M;
ax.XTickLabel = arrayfun(@(i) ['OR ' int2str(i)], 1:M, 'uniform', false);

preparegraph('edge', 0);

safe_print(fullfile('figs', 'example_env_change.pdf'));

%% Make plots using data from the fly

datafile = open('flyResponsesWithNames.mat');
data = datafile.relRates;

%fnet = open('flavornet_odor_composition.mat');

S = data';
[M, N] = size(S);

sigma = 10.0;

%Gamma = 0.00001*ones(N) + 0.01*eye(N);
Gamma = 5e-6*ones(N) + 1.5e-4*eye(N);

Ktot_values = [0 logspace(log10(0.1), log10(15000), 250)];
[K, info_values, A, info_fct] = calculate_optimal_dist(S/sigma, Gamma, Ktot_values); %#ok<ASGLU>

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
% calculate information values with subsets of receptors
% 1st column: only using the one with smallest minK
% 2nd column: using the two receptors with smallest minK
% etc.
% info_values_sub = repmat(info_values(:), [1 M-1]);
% for i = 1:M-1
%     subS = S(rec_order(1:i), :);
%     [~, info_values_sub(:, i)] = calculate_optimal_dist(subS/sigma, Gamma, Ktot_values);
% end

%% Plot "natural" receptor usage

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 7.5 2];

% plot the receptor distribution
subplot(1, 3, 1);
idx_max = 62;
% scale noise to bring K values to reasonable range
K_scaling = 50;
h_area = area(K_scaling*Ktot_values(1:idx_max), K_scaling*K(:, 1:idx_max)', 'linewidth', 1);
hold on;
for i = 1:M-1
    if minK_ordered(i+1) <= length(Ktot_values) && minK_ordered(i+1) <= idx_max
        line(repmat(K_scaling*Ktot_values(minK_ordered(i+1)), [2 1]), ylim, ...
            'linestyle', ':', 'color', [0.7 0.7 0.7], 'linewidth', 1);
    end
end
axis equal;
xlim([0 K_scaling*Ktot_values(idx_max)]);
ylim([0 K_scaling*Ktot_values(idx_max)]);

xlabel('OSN number');
ylabel('OSNs by receptor type');
% find out which receptors we've displayed, and in what order they appeared
rec_dispd = rec_order(minK_ordered <= idx_max);
h_leg = legend(h_area(rec_dispd), datafile.orNames(rec_dispd), 'location', 'northwest', ...
    'fontsize', 6);
h_leg.Position = [h_leg.Position(1) + 0.005 h_leg.Position(2) + 0.05 h_leg.Position(3:4)];
legend('boxoff');
beautifygraph;

ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;

% plot the information per OSN
subplot(1, 3, 2);
h_area = area(K_scaling*Ktot_values, K_scaling*K', 'linewidth', 0.5);
hold on;
% for i = 1:M-1
%     if minK_ordered(i+1) <= length(Ktot_values)
%         line(repmat(K_scaling*Ktot_values(minK_ordered(i+1)), [2 1]), ylim, ...
%             'linestyle', ':', 'color', [0.7 0.7 0.7], 'linewidth', 1);
%     end
% end
axis equal;
xlim([0 K_scaling*Ktot_values(end)]);
ylim([0 K_scaling*Ktot_values(end)]);

xlabel('OSN number');
ylabel('OSNs by receptor type');
% find out which receptors we've displayed, and in what order they appeared
% rec_dispd = rec_order(minK_ordered <= idx_max);
% legend(h_area(rec_dispd), datafile.orNames(rec_dispd), 'location', 'northwest', ...
%     'fontsize', 8);
% legend('boxoff');
beautifygraph;

min_nz_value = min(ax.XTick(ax.XTick > 0));
scale_exp = floor(log10(min_nz_value));
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
ticks_to_str = @(x) iif(abs(x) > eps, ...
    [num2str(x/10^scale_exp) '\times10^{' int2str(scale_exp) '}'], ...
    true, '0');

ax = gca;
ax.XTick = [0 4*10^5];
ax.XTickLabel = arrayfun(ticks_to_str, ax.XTick, 'uniform', false);
%ax.YTickLabel = arrayfun(ticks_to_str, ax_high_snr.YTick, 'uniform', false);

ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;

subplot(1, 3, 3);
rectype_count = zeros(size(Ktot_values));
for i = 1:length(minK_ordered)
    rectype_count(minK_ordered(i):end) = rectype_count(minK_ordered(i):end) + 1;
end
semilogx(K_scaling*Ktot_values, rectype_count, 'linewidth', 1);
xlabel('OSN number');
ylabel('Active receptor types');
xlim([K_scaling*Ktot_values(2), K_scaling*Ktot_values(end)]);
ylim([0 M+1]);
beautifygraph;

ax = gca;
ax.XTick = [10^1 10^2 10^3 10^4 10^5];

ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;

preparegraph('edge', 0);

safe_print(fullfile('figs', 'natural_receptor_dist.pdf'));

%% Plot showing general form of solution

range = 1:150;

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 3 2];

plot(K_scaling*Ktot_values(range), K(:, range), 'linewidth', 1);

%title('Receptor numbers');
xlabel('K');
ylabel('K_\alpha');
legend([datafile.orNames(1:3) ; {'...'}], 'location', 'northwest');

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'distribution_example.pdf'));

%% Fly: effect of changing environment 

rng('default');
%odor1 = [1 0 zeros(1, N-2)];
%odor2 = [0 1 zeros(1, N-2)];
odor1 = zeros(1, N);
odor2 = zeros(1, N);
odor1(rand(N, 1) < 0.1) = 1;
odor2(rand(N, 1) < 0.1) = 1;
Gamma1 = Gamma + 0.1*diag(odor1);
Gamma2 = Gamma + 0.1*diag(odor2);

Ktot_values = logspace(log10(0.1), log10(15000), 50);
optim_options = {'Algorithm', 'sqp', 'MaxFunctionEvaluations', 50000, ...
    'Display', 'notify-detailed', 'StepTolerance', 1e-14};
[K1, info_values_1, A1, ~] = calculate_optimal_dist(S/sigma, Gamma1, Ktot_values, ...
    optim_options{:});
[K2, info_values_2, A2, ~] = calculate_optimal_dist(S/sigma, Gamma2, Ktot_values, ...
    optim_options{:});

%% Make plot for changing environment

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 6 2];

lowSnrIdx = 24;
highSnrIdx = 50;

% show change at low SNR
subplot(1, 2, 1);
% bar([K1(:, lowSnrIdx) K2(:, lowSnrIdx)]/Ktot_values(lowSnrIdx));
bar_width = 0.8;
bar_col1 = [0.21 0.17 0.53];
bar_col2 = [0.98 0.80 0.17];
bar((1:M) - bar_width/4, K1(:, lowSnrIdx)/Ktot_values(lowSnrIdx), bar_width/2, ...
    'facecolor', bar_col1);
hold on;
bar((1:M) + bar_width/4, K2(:, lowSnrIdx)/Ktot_values(lowSnrIdx), bar_width/2, ...
    'facecolor', bar_col2);
ylabel('Fraction of OSN');
legend({'env 1', 'env 2'});
legend('boxoff');
title('Low SNR');
xlim([0 M+bar_width]);
ylim([0 0.2]);
beautifygraph;
ax = gca;
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;
ax.XTick = 1:M;
ax.XTickLabel = datafile.orNames;
ax.XTickLabelRotation = 90;
ax.FontSize = 7;

subplot(1, 2, 2);
% bar([K1(:, highSnrIdx) K2(:, highSnrIdx)]/Ktot_values(highSnrIdx));
bar((1:M) - bar_width/4, K1(:, highSnrIdx)/Ktot_values(highSnrIdx), bar_width/2, ...
    'facecolor', bar_col1);
hold on;
bar((1:M) + bar_width/4, K2(:, highSnrIdx)/Ktot_values(highSnrIdx), bar_width/2, ...
    'facecolor', bar_col2);
xlim([0 M+bar_width]);
ylim([0 0.2]);
ylabel('Fraction of OSN');
title('High SNR');
legend({'env 1', 'env 2'});
legend('boxoff');
beautifygraph;
ax = gca;
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;
ax.XTick = 1:M;
ax.XTickLabel = datafile.orNames;
ax.XTickLabelRotation = 90;
ax.FontSize = 7;

preparegraph('edge', 0);

safe_print(fullfile('figs', 'natural_env_change.pdf'));

%% Different kind of plot

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 6 2];

lowSnrIdx = 24;
highSnrIdx = 50;

% show change at low SNR
subplot(1, 2, 1);
plot([zeros(1, M) ; ones(1, M)], ...
    [K1(:, lowSnrIdx)'/Ktot_values(lowSnrIdx) ; K2(:, lowSnrIdx)'/Ktot_values(lowSnrIdx)], ...
    'k-');
hold on;
plot(zeros(1, M), K1(:, lowSnrIdx)'/Ktot_values(lowSnrIdx), 'd', ...
    'markerfacecolor', [0.21 0.17 0.53], 'markeredgecolor', 'none');
plot(ones(1, M), K2(:, lowSnrIdx)'/Ktot_values(lowSnrIdx), 'd', ...
    'markerfacecolor', [0.98 0.80 0.17], 'markeredgecolor', 'none');
ylabel('Fraction of OSN');
title('Low SNR');
ylim([0 0.1]);
beautifygraph;
ax = gca;
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = [0 0];
ax.XTick = [0 1];
ax.XTickLabel = {'env1', 'env2'};

subplot(1, 2, 2);
plot([zeros(1, M) ; ones(1, M)], ...
    [K1(:, highSnrIdx)'/Ktot_values(highSnrIdx) ; K2(:, highSnrIdx)'/Ktot_values(highSnrIdx)], ...
    'k-');
hold on;
plot(zeros(1, M), K1(:, highSnrIdx)'/Ktot_values(highSnrIdx), 'd', ...
    'markerfacecolor', [0.21 0.17 0.53], 'markeredgecolor', 'none');
plot(ones(1, M), K2(:, highSnrIdx)'/Ktot_values(highSnrIdx), 'd', ...
    'markerfacecolor', [0.98 0.80 0.17], 'markeredgecolor', 'none');
ylim([0 0.1]);
ylabel('Fraction of OSN');
title('High SNR');
beautifygraph;
ax = gca;
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = [0 0];
ax.XTick = [0 1];
ax.XTickLabel = {'env1', 'env2'};

preparegraph('edge', 0);

safe_print(fullfile('figs', 'natural_env_change_alt.pdf'));

%% Differenter kind of plot

diff_vec = K2 - K1;

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 3 2];

K_scaling = 50;
plot(K_scaling*Ktot_values, K_scaling*diff_vec, 'linewidth', 1);
xlabel('K_{tot}');
ylabel('\Delta K_\alpha');

xlim([0, max(K_scaling*Ktot_values)]);

ax = gca;
ax.Box = 'off';
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = [0 0];

beautifygraph;
preparegraph;

safe_print(fullfile('figs', 'abs_diff_env_change.pdf'));

%% Show the shape of the mutual information surface

datafile = open('flyResponsesWithNames.mat');
data = datafile.relRates;

%fnet = open('flavornet_odor_composition.mat');

S = data';
[M, N] = size(S);

sigma = 10.0;

%Gamma = 0.00001*ones(N) + 0.01*eye(N);
Gamma = 5e-6*ones(N) + 1.5e-4*eye(N);

% calculate the receptor distribution
Ktot = 1500;
[K, info_max, Q, info_fct] = calculate_optimal_dist(S/sigma, Gamma, Ktot);
%A = inv(Q);

% sample nearby points
% n_pts = 1000;
% K_pts = zeros(M, n_pts);
% info_pts = zeros(n_pts, 1);
% for i = 1:n_pts
%      K_trial = max(0, K + 0.01*max(K)*rand(size(K)));
%      K_trial = K_trial/sum(K_trial)*sum(K);
%      K_pts(:, i) = K_trial;
%      info_pts(i) = info_fct(K_pts(:, i));
% end

% calculate hessian
L = inv(eye(size(Q)) + diag(K)*Q);
hessian = -0.5*L.*L'; %#ok<MINV>

% find directions of maximum variation
[pcs, ~] = eig(hessian);
K_pc_rng = 0.01*max(K);
x = linspace(-K_pc_rng, K_pc_rng, 100);
y = linspace(-K_pc_rng, K_pc_rng, 100);
info_surf = zeros(length(x), length(y));
v1 = pcs(:, 1) - ones(M, 1)*(ones(1, M)*pcs(:, 1))/M;
v1 = v1/norm(v1);
v2 = pcs(:, 2) - ones(M, 1)*(ones(1, M)*pcs(:, 2))/M;
v2 = v2/norm(v2);
for i = 1:length(x)
    for j = 1:length(y)
        crtK = K + v1*x(i) + v2*y(j);
        if any(crtK < 0)
            info_surf(i, j) = nan;
        else
            info_surf(i, j) = info_fct(crtK);
        end
    end
end

%% Covariance matrices plots

odor_pic_names = {'1-hexanol', 'E2-hexenol', 'E3-hexenol', 'Z3-hexenol', ...
    '2-heptanone', 'butyl acetate', 'pentyl acetate', 'methyl hexanoate', ...
    'ethyl hexanoate', 'methyl benzoate', 'ethyl benzoate'};
odor_pic_idxs = zeros(size(odor_pic_names));
for i = 1:length(odor_pic_names)
    odor_pic_idxs(i) = find(strcmp(datafile.odorNames, odor_pic_names{i}));
end

%cticks = [1e-5/3, 1e-5, 1e-4/3, 1e-4, 1e-3/3, 1e-3, 1e-2/3, 1e-2, 1e-1/3, 1e-1];
cticks = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1];
cticklabels = {'10^{-5}', '10^{-4}', '0.001', '0.01', '0.1'};
clims = [1e-5/3, 1e-1];

fig1 = figure;
fig1.Units = 'inches';
fig1.Position = [fig1.Position(1:2) 2 2];

imagesc(log10(Gamma1(odor_pic_idxs, odor_pic_idxs)), log10(clims));

hold on;
n_odor = length(odor_pic_idxs);
for i = 0:n_odor
    line([0.5, n_odor + 0.5], i + [0.5, 0.5], 'color', 'w', 'linewidth', 2);
    line(i + [0.5, 0.5], [0.5, n_odor + 0.5], 'color', 'w', 'linewidth', 2);
end

beautifygraph;

axis equal;
axis off;

preparegraph;

safe_print(fullfile('figs', 'odor_cov_env1.pdf'));

fig2 = figure;
fig2.Units = 'inches';
fig2.Position = [fig2.Position(1:2) 2 2];

h = imagesc(log10(Gamma2(odor_pic_idxs, odor_pic_idxs)), log10(clims));

hold on;
n_odor = length(odor_pic_idxs);
for i = 0:n_odor
    line([0.5, n_odor + 0.5], i + [0.5, 0.5], 'color', 'w', 'linewidth', 2);
    line(i + [0.5, 0.5], [0.5, n_odor + 0.5], 'color', 'w', 'linewidth', 2);
end

beautifygraph;

axis equal;
axis off;

preparegraph;

safe_print(fullfile('figs', 'odor_cov_env2.pdf'));

fig3 = figure;
fig3.Units = 'inches';
fig3.Position = [fig3.Position(1:2) 1.5 2];

beautifygraph;

axis off;
caxis(log10(clims));
hcb = colorbar('ytick', log10(cticks), 'yticklabel', cticklabels);

preparegraph;

safe_print(fullfile('figs', 'odor_cov_colorbar.pdf'));
