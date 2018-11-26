% Make plots showing an example of the effects of environment change on
% receptor distribution.
%
% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   S_choice:
%       Choice of sensing matrix to use. Options are:
%         'fly'     -- Drosophila sensing matrix from Hallem&Carlson (2006)
%         'mouse'   -- mouse sensing matrix from Saito et al. (2009)
%         'fly_scrambled'
%                   -- fly sensing matrix randomly scrambled across both
%                      rows and columns
%         'mouse_scrambled'
%                   -- mouse sensing matrix randomly scrambled across both
%                      rows and columns
%         'artificial'
%                   -- an artificial sensing matrix using a range of tuning
%                      widths
%         'artificial_narrow'
%                   -- artificial sensing matrix with narrowly-tuned
%                      receptors
%         'artificial_wide'
%                   -- artificial sensing matrix with broadly-tuned
%                      receptors
%         'gaussian'-- artificial sensing matrix with entries drawn from a
%                      zero-centered Gaussian
%         'binary'  -- artificial sensing matrix with binary 0 or 1 entries
%         'signed'  -- artificial sensing matrix with entries that are
%                      either 0 or +1/-1
%   S_size:
%       Size of sensing matrix, [n_receptors, n_odorants]. This is usually
%       [24, 110] (the size of the Drosophila matrix), unless `S_choice` is
%       'mouse' or 'mouse_scrambled'.
%   artif_tuning:
%       Range of tuning width used with `generate_random_sensing` when
%       `S_choice` is one of the artificial settings. Only the first (last)
%       element of `artif_tuning` is used for 'artificial_narrow'
%       ('artifical_wide').
%   artif_snr:
%       SNR setting to use with `generate_random_sensing` when `S_choice`
%       is one of the artificial settings.

%% Setup

setdefault('S_choice', 'fly');
setdefault('S_size', [24, 110]);
setdefault('artif_tuning', [0.2 0.8]);
setdefault('artif_snr', 200);

%% Preprocess some options

if ~strcmp(S_choice, 'fly')
    postfix = ['_' S_choice];
else
    postfix = '';
end

%% Load fly (Hallem&Carlson) sensing data

sensing_fly = open('data/flyResponsesWithNames.mat');
S_fly = sensing_fly.relRates';

% normalize by standard deviations of background rates
S_fly_normalized = bsxfun(@rdivide, S_fly, sensing_fly.bkgStd');

optim_args = {'optimopts', ...
        {'MaxFunctionEvaluations', 50000, 'Display', 'notify-detailed'}, ...
    'method', 'lagsearch'};

% generate colormap for covariance matrix plots
% cmap_covmat = divergent([0.21 0.17 0.53], [0.98 0.40 0.17], 256);
cmap_covmat = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

%% Load mouse (Saito et al.) sensing matrix

% load response curves
mouse_responses = open(fullfile('data', 'mouse_receptor_curves.mat'));
S_mouse = mouse_responses.sensing_matrix(:, :, 11);
% set missing data to 0
S_mouse(isnan(S_mouse)) = 0;
% remove receptors that don't respond to anything
mask = ~all(S_mouse == 0, 2);
S_mouse = S_mouse(mask, :);

%% Select a sensing matrix

% keep random stuff reproducible
rng(8237);

switch S_choice
    case 'fly'
        S = S_fly_normalized;
    case 'mouse'
        S = 1500*S_mouse;
    case 'fly_scrambled'
        S = 0.1*S_fly_normalized;
        S(:) = S(randperm(numel(S)));
    case 'mouse_scrambled'
        S = 1500*S_mouse;
        S(:) = S(randperm(numel(S)));
    case 'artificial'
        S = 100*generate_random_sensing(S_size(1), S_size(2), artif_tuning, artif_snr);
    case 'artificial_narrow'
        S = 100*generate_random_sensing(S_size(1), S_size(2), artif_tuning(1), artif_snr);
    case 'artificial_wide'
        S = 100*generate_random_sensing(S_size(1), S_size(2), artif_tuning(2), artif_snr);
    case 'gaussian'
        S = 2*randn(S_size);
    case 'binary'
        S = 5*(rand(S_size) > 0.7);
    case 'signed'
        S = 5*sign(rand(S_size) - 0.5);
        S = S .* (rand(S_size) > 0.7);
    otherwise
        error('Unrecognized setting for S_choice.');
end

%% Example receptor abundances in two environments

% make this reproducible
rng(9843857);

% start from a random environment covariance matrix
Gamma0 = generate_environment('rnd_corr', size(S, 2));
% add variance to two different sets of odorants
if isequal(size(S), [24 110])
    Gamma1 = generate_environment('delta_rnd_diag', Gamma0, ...
        'delta_size', 100, 'delta_pos', [16, 21, 23, 25, 43, 46, 81, 91, 99, 105]);
    Gamma2 = generate_environment('delta_rnd_diag', Gamma0, ...
        'delta_size', 100, 'delta_pos', [5, 18, 33, 36, 46, 53, 66, 71, 84, 101, 107]);
else
    idxs_all = randperm(size(S, 2), 20);
    idxs1 = idxs_all(1:10);
    idxs2 = idxs_all(11:20);
    Gamma1 = generate_environment('delta_rnd_diag', Gamma0, ...
        'delta_size', 100, 'delta_pos', idxs1);
    Gamma2 = generate_environment('delta_rnd_diag', Gamma0, ...
        'delta_size', 100, 'delta_pos', idxs2);
end

% The concentrations in the artificial environment are chosen arbitrarily.
% We roughly fit their scale to obtain intermediate SNR at reasonable
% values of the total number of neurons, Ktot.
factor = 2e-4;
Gamma0 = factor*Gamma0;
Gamma1 = factor*Gamma1;
Gamma2 = factor*Gamma2;

% values of Ktot for low and high SNR
Ktot = [100 40000];
% get optimal repertoires
[K1, ~, Q1] = calculate_optimal_dist(S, Gamma1, Ktot, optim_args{:});
[K2, ~, Q2] = calculate_optimal_dist(S, Gamma2, Ktot, optim_args{:});

%% ... plot the change

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.8 2];

fig.Color = [1 1 1];

distChangeColors = {hex2color('2D75BB'), hex2color('EE8434')};

frac1 = bsxfun(@rdivide, K1, Ktot);
frac2 = bsxfun(@rdivide, K2, Ktot);
max_frac1 = max([frac1(:) ; frac2(:)]);

yl = 1.1485*max_frac1;

% show change at low SNR
ax1 = axes;
ax1.Units = 'normalized';
ax1.OuterPosition = [1/2 0 1/2 1];
plotDistChange(K1(:, 1), K2(:, 1), 'colors', distChangeColors, ...
    'beautifyopts', {'fontscale', 0.667});
ylim([0 yl]);
title('Low SNR');

beautifygraph('fontscale', 0.667, 'ticksize', 12, 'linewidth', 0.5);

% show change at high SNR
ax2 = axes;
ax2.Units = 'normalized';
ax2.OuterPosition = [0 0 1/2 1];
plotDistChange(K1(:, 2), K2(:, 2), 'colors', distChangeColors, ...
    'beautifyopts', {'fontscale', 0.667});
ylim([0 yl]);
title('High SNR');

beautifygraph('fontscale', 0.667, 'ticksize', 12, 'linewidth', 0.5);

preparegraph('edge', 0);

safeprint(fullfile('figs', ['natural_env_change_example' postfix]));

%% ...plot the sensing matrix

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position(3:4) = 3.5*[1 size(S, 1)/size(S, 2)];

ax = axes;
ax.Position = [0 0 1 1];

cmap = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

cl = quantile(flatten(abs(S)), 0.9);
if cl == 0
    cl = max(abs(S(:)));
end

plotheat(S, [-cl cl], 'edge', 0);
colormap(cmap);

axis off;

preparegraph;

safeprint(fullfile('figs', ['example_full_sensing' postfix]));

%% ... plot covariance matrices used for the example

% NOTE: this doesn't make sense unless we use the fly covariance matrices

if isequal(size(S), [24 110])
    % we have diagrams for these odorants
    odor_pic_names = {'1-hexanol', 'E2-hexenol', 'E3-hexenol', 'Z3-hexenol', ...
        '2-heptanone', 'butyl acetate', 'pentyl acetate', 'methyl hexanoate', ...
        'ethyl hexanoate', 'methyl benzoate', 'ethyl benzoate'};
    odor_pic_idxs = zeros(size(odor_pic_names));
    for i = 1:length(odor_pic_names)
        odor_pic_idxs(i) = find(strcmp(sensing_fly.odorNames, odor_pic_names{i}));
    end
    
    % a useful set of functions below
    iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
    ticks_to_str = @(x, scale_exp) iif(abs(x) > eps, ...
        [num2str(x/10^scale_exp) '\times10^{' int2str(scale_exp) '}'], ...
        true, '0');
    
    % choose a scale for the colormap, and locations for ticks
    cticks = factor*[-1.5, -1, -0.5, 0, 0.5, 1, 1.5];
    cticklabels = arrayfun(@(x) ticks_to_str(x, -4), cticks, 'uniform', false);
    clims = factor*[-1.5 1.5];
    
    % draw the zoomed-in detail of the covariance matrix for environment 1
    fig1 = figure;
    fig1.Units = 'inches';
    fig1.Position = [fig1.Position(1:2) 1.9 1.9];
    
    fig1.Color = [1 1 1];
    
    imagesc(Gamma1(odor_pic_idxs, odor_pic_idxs), clims);
    colormap(cmap_covmat);
    
    % overlay a white grid, so that each element is separated
    hold on;
    n_odor = length(odor_pic_idxs);
    for i = 0:n_odor
        plot([0.5, n_odor + 0.5], i + [0.5, 0.5], 'color', [1 1 1], 'linewidth', 2);
        plot(i + [0.5, 0.5], [0.5, n_odor + 0.5], 'color', [1 1 1], 'linewidth', 2);
    end
    
    beautifygraph;
    
    axis equal;
    axis off;
    
    preparegraph;
    
    safeprint(fullfile('figs', 'odor_cov_example_env1_zoom'));
    
    % draw the zoomed-in detail of the covariance matrix for environment 1
    fig2 = figure;
    fig2.Units = 'inches';
    fig2.Position = [fig2.Position(1:2) 1.9 1.9];
    
    fig2.Color = [1 1 1];
    
    imagesc(Gamma2(odor_pic_idxs, odor_pic_idxs), clims);
    colormap(cmap_covmat);
    
    % overlay a white grid, so that each element is separated
    hold on;
    n_odor = length(odor_pic_idxs);
    for i = 0:n_odor
        plot([0.5, n_odor + 0.5], i + [0.5, 0.5], 'color', [1 1 1], 'linewidth', 2);
        plot(i + [0.5, 0.5], [0.5, n_odor + 0.5], 'color', [1 1 1], 'linewidth', 2);
    end
    
    beautifygraph;
    
    axis equal;
    axis off;
    
    preparegraph;
    
    safeprint(fullfile('figs', 'odor_cov_example_env2_zoom'));
    
    % draw the color bar
    fig3 = figure;
    fig3.Units = 'inches';
    fig3.Position = [fig3.Position(1:2) 1.5 1.9];
    
    fig3.Color = [1 1 1];
    
    beautifygraph;
    
    axis off;
    caxis(clims);
    colormap(cmap_covmat);
    hcb = colorbar('west', 'ytick', cticks, 'yticklabel', cticklabels);
    
    preparegraph;
    
    safeprint(fullfile('figs', 'odor_cov_example_zoom_colorbar'));
    
    % draw the full background covariance matrix, Gamma0_example
    fig4 = figure;
    fig4.Units = 'inches';
    fig4.Position = [fig4.Position(1:2) 1.9 1.9];
    
    fig4.Color = [1 1 1];
    
    % reorder the indices so that the ones that we have diagrams for are
    % bunched together in the middle (that was it's easier to show how we zoom
    % in)
    reorder = [1:min(odor_pic_idxs)-1 odor_pic_idxs ...
        setdiff(setdiff(1:size(Gamma0, 1), odor_pic_idxs), 1:min(odor_pic_idxs)-1)];
    box_start = min(odor_pic_idxs) - 0.5;
    box_end = min(odor_pic_idxs)+length(odor_pic_idxs)-0.5;
    
    plotheat(Gamma0(reorder, reorder), clims);
    colormap(cmap_covmat);
    
    % draw a box showing where the zoomed-in matrices are located
    hold on;
    % draw a thick outline in white, and a thin one in black, to improve
    % visibility
    plot([box_start box_start box_end box_end box_start], ...
        [box_start box_end box_end box_start box_start], 'w', 'linewidth', 2);
    plot([box_start box_start box_end box_end box_start], ...
        [box_start box_end box_end box_start box_start], 'k');
    
    beautifygraph('fontscale', 0.667, 'linewidth', 0.5);
    
    axis equal;
    set(gca, 'xtick', [], 'ytick', [], 'box', 'on');
    xlabel('odorant 1');
    ylabel('odorant 2');
    
    preparegraph;
    
    safeprint(fullfile('figs', 'odor_cov_example_env0'));
end
