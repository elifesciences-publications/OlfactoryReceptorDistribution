% Make plots about the nonlinear model.

%% Load data

load(fullfile('save', 'nonlinear_example.mat'));

%% Visualize sensing matrix and environments

cmap = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [5.5 1.5];
fig.Color = [1 1 1];

sensing_frac = 0.58;

ax = axes;
ax.OuterPosition = [0 0.15 sensing_frac 1];
cl = quantile(abs(S(:)), 0.95);
plotheat(S, [-cl, cl]);
colormap(cmap);
% colorbar('southoutside');
colorbar('northoutside', 'fontsize', 6.7);
xlabel('Odorants');
ylabel('Receptors');
% title('Sensing matrix');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

cl = quantile(abs([Gamma1(:) ; Gamma2(:)]), 0.95);

env_frac = (1 - sensing_frac) / 2;
ax = axes;
ax.OuterPosition = [sensing_frac 0 env_frac 1];
plotheat(Gamma1, [-cl, cl]);
colormap(cmap);
colorbar('fontsize', 6.7);
% title('Environment matrix 1');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

ax = axes;
ax.OuterPosition = [sensing_frac + env_frac 0 env_frac 1];
plotheat(Gamma2, [-cl, cl]);
colormap(cmap);
colorbar('fontsize', 6.7);
% title('Environment matrix 2');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);
preparegraph('edge', 0);

safeprint(fullfile('figs', 'nonlinear_setting'));

%% Show exact response distribution in the linear case

bin_edges = linspace(-1.5, 4.0, 101);
n_bins = length(bin_edges)-1;

% do it in a dumb way, for simplicity
Pr_exact1 = zeros(n_bins);
Pr_exact2 = zeros(n_bins);
mean_exact = S*c0(:);

% components of Pr are actually histogram counts, need to multiply by bin size
factor = diff(bin_edges(1:2))^n_receptors;
for i = 1:n_bins
    col = 0.5*(bin_edges(i) + bin_edges(i+1));
    for j = 1:n_bins
        row = 0.5*(bin_edges(j) + bin_edges(j+1));
        for k = 1:n_bins
            depth = 0.5*(bin_edges(k) + bin_edges(k+1));
            Pr_exact1(i, j, k) = factor*mvnpdf([col row depth], ...
                mean_exact', cov_exact1);
            Pr_exact2(i, j, k) = factor*mvnpdf([col row depth], ...
                mean_exact', cov_exact2);
        end
    end
end

%% Plot exact distribution

fig = figure;
fig.Color = [1 1 1];
fig.Position(3) = 3*fig.Position(3);
fig.Position(4) = 2*fig.Position(4);

for j = 1:2
    switch j
        case 1
            Pr_exact = Pr_exact1;
        case 2
            Pr_exact = Pr_exact2;
    end
    for i = 1:3
        ax = axes;
        ax.OuterPosition = [(i-1)/3 1-j/2 1/3 1/2];
        
        switch i
            case 1
                imagesc(bin_edges([1 end]), bin_edges([1 end]), ...
                    squeeze(Pr_exact(:, :, end/2)));
                lab1 = 'r_1';
                lab2 = 'r_2';
            case 2
                imagesc(bin_edges([1 end]), bin_edges([1 end]), ...
                    squeeze(Pr_exact(:, end/2, :)));
                lab1 = 'r_1';
                lab2 = 'r_3';
            case 3
                imagesc(bin_edges([1 end]), bin_edges([1 end]), ...
                    squeeze(Pr_exact(end/2, :, :)));
                lab1 = 'r_2';
                lab2 = 'r_3';
        end
        colormap(flipud(gray));
        set(gca, 'ydir', 'normal');
%         set(gca, 'xticklabel', arrayfun(@(x) num2str(bin_edges(x), '%.2f'), get(gca, 'xtick'), 'uniform', false));
%         set(gca, 'yticklabel', arrayfun(@(x) num2str(bin_edges(x), '%.2f'), get(gca, 'ytick'), 'uniform', false));
        title('Probability distribution (exact)');
        
        ylabel(lab1);
        xlabel(lab2);
        
        beautifygraph('linewidth', 0.5);
    end
end

preparegraph;

%% Find optimal receptor distribution using linear method, and plot

K1linear = calculate_optimal_dist(S./sigma(:), Gamma1, Ktot);
K2linear = calculate_optimal_dist(S./sigma(:), Gamma2, Ktot);

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position(3:4) = [2 2];

ax = axes;
ax.OuterPosition = [0 1/2 1 1/2];
bar(K1linear, 'edgecolor', 'none');
xlabel('receptor');
ylabel('abundance');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

ax = axes;
ax.OuterPosition = [0 0 1 1/2];
bar(K2linear, 'edgecolor', 'none');
xlabel('receptor');
ylabel('abundance');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

preparegraph('edge', 0);

%% Make information plots for the nonlinear model

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position(3:4) = [2 2];

imagesc(Kbins([1 end]), Kbins([1 end]), info_map_nonlinear1, 'alphadata', ...
    ~isnan(info_map_nonlinear1));
% colormap(flipud(gray));
colormap(cmap(end/2:end, :));
colorbar;
axis image;
xlabel('K_2');
ylabel('K_1');

set(gca, 'ydir', 'normal');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

preparegraph('edge', 0);

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position(3:4) = [2 2];

imagesc(Kbins([1 end]), Kbins([1 end]), info_map_nonlinear2, 'alphadata', ...
    ~isnan(info_map_nonlinear2));
% colormap(flipud(gray));
colormap(cmap(end/2:end, :));
colorbar;
axis image;
xlabel('K_2');
ylabel('K_1');

set(gca, 'ydir', 'normal');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

preparegraph('edge', 0);

%% Combine plots

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position(3:4) = [4 3];

ax = axes;
ax.OuterPosition = [0.01 0.78 0.35 0.2];
bar(K1linear, 'edgecolor', 'none', 'facecolor', hex2color('2D75BB'));
xlabel('receptor');
ylabel('abundance');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

ax = axes;
ax.OuterPosition = [0 0.12 0.5 0.8];
imagesc(Kbins([1 end]), Kbins([1 end]), info_map_nonlinear1, 'alphadata', ...
    ~isnan(info_map_nonlinear1));
% colormap(flipud(gray));
colormap(cmap(end/2:end, :));
colorbar;
axis image;
xlabel('K_2');
ylabel('K_1');

set(gca, 'ydir', 'normal');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

[~, max_idx1] = max(info_map_nonlinear1(:));
[tmpK1idx, tmpK2idx] = ind2sub(size(info_map_nonlinear1), max_idx1);
tmpK1 = Kbins(tmpK1idx);
tmpK2 = Kbins(tmpK2idx);
K1nonlinear = [tmpK1 tmpK2 Ktot - (tmpK1 + tmpK2)];

ax = axes;
ax.OuterPosition = [0.01 0 0.35 0.2];
bar(K1nonlinear, 'edgecolor', 'none', 'facecolor', hex2color('EE8434'));
xlabel('receptor');
ylabel('abundance');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

% second environment

ax = axes;
ax.OuterPosition = [0.56 0.78 0.35 0.2];
bar(K2linear, 'edgecolor', 'none', 'facecolor', hex2color('2D75BB'));
xlabel('receptor');
ylabel('abundance');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

ax = axes;
ax.OuterPosition = [0.55 0.12 0.5 0.8];
imagesc(Kbins([1 end]), Kbins([1 end]), info_map_nonlinear2, 'alphadata', ...
    ~isnan(info_map_nonlinear1));
% colormap(flipud(gray));
colormap(cmap(end/2:end, :));
colorbar;
axis image;
xlabel('K_2');
ylabel('K_1');

set(gca, 'ydir', 'normal');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

[~, max_idx2] = max(info_map_nonlinear2(:));
[tmpK1idx, tmpK2idx] = ind2sub(size(info_map_nonlinear2), max_idx2);
tmpK1 = Kbins(tmpK1idx);
tmpK2 = Kbins(tmpK2idx);
K2nonlinear = [tmpK1 tmpK2 Ktot - (tmpK1 + tmpK2)];

ax = axes;
ax.OuterPosition = [0.56 0 0.35 0.2];
bar(K2nonlinear, 'edgecolor', 'none', 'facecolor', hex2color('EE8434'));
xlabel('receptor');
ylabel('abundance');

beautifygraph('minorticks', 'off', 'linewidth', 0.5, ...
    'fontscale', 0.667, 'ticksize', 10, 'labelsize', 10);

preparegraph('edge', 0);

safeprint(fullfile('figs', 'nonlinear_results'));
