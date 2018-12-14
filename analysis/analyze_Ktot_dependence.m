% Check how the number of receptor types with non-zero OSN populations
% depends on the total number of neurons for a variety of sensing matrices.
%
% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   S_choice:
%       Choices of sensing matrix to use. This should be a cell array of
%       options from the following list, or 'all' to use all of them.
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
%   Ktot_range:
%       Range for the total number of neurons to use in optimizations.
%   n_Ktot:
%       Number of Ktot values to consider. The value Ktot = 0 is always
%       added to this.

%% Setup

setdefault('S_choice', 'all');
setdefault('S_size', [24, 110]);
setdefault('artif_tuning', [0.2 0.8]);
setdefault('artif_snr', 200);
setdefault('Ktot_range', [5, 7.5e5]);
setdefault('n_Ktot', 100);

%% Postprocess some options

if strcmp(S_choice, 'all')
    S_choice = {'fly', 'fly_scrambled', 'mouse', 'mouse_scrambled', ...
        'artificial', 'artificial_narrow', 'artificial_wide', ...
        'gaussian', 'binary', 'signed'};
end

%% Setup some non-configurable settings

optim_args = {'optimopts', ...
        {'MaxFunctionEvaluations', 50000, 'Display', 'notify-detailed'}, ...
    'method', 'lagsearch', 'sumtol', 1e-3};

%% Load fly (Hallem&Carlson) sensing data

sensing_fly = open('data/flyResponsesWithNames.mat');
S_fly = sensing_fly.relRates';

% normalize by standard deviations of background rates
S_fly_normalized = bsxfun(@rdivide, S_fly, sensing_fly.bkgStd');

%% Load mouse (Saito et al.) sensing matrix

% load response curves
mouse_responses = open(fullfile('data', 'mouse_receptor_curves.mat'));
S_mouse = mouse_responses.sensing_matrix(:, :, 11);
% set missing data to 0
S_mouse(isnan(S_mouse)) = 0;
% remove receptors that don't respond to anything
mask = ~all(S_mouse == 0, 2);
S_mouse = S_mouse(mask, :);

%% Calculate the receptor abundances at a series of Ktot values

% select the values of Ktot that we'll be considering
Ktot_values = [0 logspace(log10(Ktot_range(1)), log10(Ktot_range(2)), n_Ktot)];

% an extended array to be used for receptors that are never activated
% within the considered range of Ktot
Ktot_values_ext = [Ktot_values inf];

results = containers.Map;
for i = 1:length(S_choice)
    % select a sensing matrix
    
    % keep random stuff reproducible
    rng(8237);
    
    crt_S_choice = S_choice{i};
    
    switch crt_S_choice
        case 'fly'
            S = S_fly_normalized;
        case 'mouse'
            % adding a bit of random noise to remove degeneracies from the very
            % sparse S_mouse (rank(S_mouse) = 50 while size(S_mouse) = [59, 63]).
            S = 50*(S_mouse + 0.01*randn(size(S_mouse)));
        case 'fly_scrambled'
            S = 0.1*S_fly_normalized;
            S(:) = S(randperm(numel(S)));
        case 'mouse_scrambled'
            S = 50*S_mouse;
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

    % make environment
    
    % random environment, but make it reproducible
    rng(92376);
    
    Gamma = generate_environment('rnd_diag_const', size(S, 2), ...
        'diag_mu', 1e-6, 'diag_size', 1e-6, 'offdiag_mu', 1e-7);
    
    % get the optimal distribution at all Ktot values
    disp(['Working on ' crt_S_choice '...']);
    [K, info_values] = calculate_optimal_dist(S, Gamma, Ktot_values);
    
    % find where each receptor kicks in
    M = size(S, 1);
    
    % find where each receptor kicks in
    minK_by_type = ones(M, 1);
    K_nz = bsxfun(@rdivide, K, Ktot_values) > 0.004;
    for j = 1:M
        idx = find(~K_nz(j, :), 1, 'last');
        if ~isempty(idx)
            minK_by_type(j) = idx+1;
        end
    end
    [minK_ordered, rec_order] = sort(minK_by_type);
    
    % make a structure of results for this particular sensing matrix
    crt_results = struct;
    crt_results.S = S;
    crt_results.Gamma = Gamma;
    crt_results.K = K;
    crt_results.info_values = info_values;
    crt_results.Ktot_values = Ktot_values;
    % starting index in Ktot_values array where each receptor type starts
    % getting expressed
    crt_results.starting_Ktot_idx = minK_by_type;
    % the corresponding Ktot value
    crt_results.starting_Ktot = Ktot_values_ext(minK_by_type);
    % indices in Ktot_values where more receptor types are added
    crt_results.steps_Ktot_idx = minK_ordered;
    % Ktot values where the number of expressed receptor types increases
    crt_results.steps_Ktot = Ktot_values_ext(minK_ordered);
    % the order in which receptors get expressed
    crt_results.expression_order = rec_order;
    crt_results.rec_order = rec_order;
    
    results(crt_S_choice) = crt_results;
end

%% Plot number of active receptors vs. Ktot for all sensing matrices

n_plots = length(S_choice);
nx = ceil(sqrt(n_plots));
ny = ceil(n_plots / nx);

ax_x = min(3, 12/nx);
ax_y = ax_x*3/4;

fig_x = ax_x*nx;
fig_y = ax_y*ny;

fig = figure;
fig.Units = 'inches';
fig.Position = [1 1 fig_x fig_y];
fig.Color = [1 1 1];

for i = 1:n_plots
    crt_S_choice = S_choice{i};
    
    crt_x = mod(i-1, nx);
    crt_y = floor((i-1)/nx);
    
    ax = axes;
    ax.OuterPosition = [crt_x/nx 1-(crt_y+1)/ny 1/nx 1/ny];

    % make plot showing how the number of receptor expressed in non-zero
    % numbers of neurons grows with the total number of neurons
    crt_results = results(crt_S_choice);
    crt_Ktot_values = crt_results.Ktot_values;
    crt_minK_ordered = crt_results.steps_Ktot_idx;
    crt_S = crt_results.S;
    crt_K = crt_results.K;
    
%     crt_K_nz = (bsxfun(@rdivide, crt_K, crt_Ktot_values) > 0.004);
    crt_K_nz = (bsxfun(@rdivide, crt_K, crt_Ktot_values) > 0.01);
    
    rectype_count = sum(crt_K_nz, 1);
    semilogx(crt_Ktot_values, rectype_count, 'linewidth', 1);
    xlabel('Number of OSNs');
    ylabel('Active receptor types');
    xlim([crt_Ktot_values(2), crt_Ktot_values(end)]);
    ylim([0 size(crt_S, 1)+1]);
    crt_name = crt_S_choice;
    crt_name(crt_name == '_') = ' ';
    title(crt_name);
    beautifygraph(ax, 'fontscale', 0.667, 'ticksize', 12, 'linewidth', 0.5);

    % select a good crop of ticks
    ax.XTick = [10^1 10^2 10^3 10^4 10^5];

    ax.XMinorTick = 'off';
    ax.YMinorTick = 'off';
    ax.TickLength = 3*ax.TickLength;
end
preparegraph('edge', 0);

%% Save the results

save(fullfile('save', 'Ktot_dependence.mat'), 'results', 'Ktot_values', ...
    'Ktot_range', 'n_Ktot', 'S_choice', 'S_size', 'artif_snr', 'artif_tuning', ...
    'optim_args');
