% Show examples of context dependence, i.e., the same perturbation leads to
% different effects depending on starting environment.
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
%   Ktot:
%       Total number of neurons to use in optimizations.
%   corr_beta:
%       Parameter to use for generating the background environment matrices.
%   n_samples:
%       Number of environment samples to consider.

%% Setup

setdefault('S_choice', 'fly');
setdefault('S_size', [24, 110]);
setdefault('artif_tuning', [0.2 0.8]);
setdefault('artif_snr', 200);
setdefault('Ktot', 2);
setdefault('corr_beta', 8);
setdefault('n_samples', 100);

%% Preprocess some options

if ~strcmp(S_choice, 'fly')
    postfix = ['_' S_choice];
else
    postfix = '';
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

%% Select a sensing matrix

% keep random stuff reproducible
rng(8237);

switch S_choice
    case 'fly'
        S = S_fly_normalized;
    case 'mouse'
        % adding a bit of random noise to remove degeneracies from the very
        % sparse S_mouse (rank(S_mouse) = 50 while size(S_mouse) = [59, 63]).
        S = 5000*(S_mouse + 0.01*randn(size(S_mouse)));
    case 'fly_scrambled'
        S = 0.1*S_fly_normalized;
        S(:) = S(randperm(numel(S)));
    case 'mouse_scrambled'
        S = 5000*S_mouse;
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

S_size = size(S);

%% Generate environments and perturb both in the same way

% random environments, but make them reproducible
rng(32392076);

results.Gamma1 = cell(1, n_samples);
results.Gamma2 = cell(1, n_samples);

results.idx_odorant = zeros(1, n_samples);
results.Gamma1_pert = cell(1, n_samples);
results.Gamma2_pert = cell(1, n_samples);

results.K1 = zeros(size(S, 1), n_samples);
results.K2 = zeros(size(S, 1), n_samples);

results.K1_pert = zeros(size(S, 1), n_samples);
results.K2_pert = zeros(size(S, 1), n_samples);

% sometimes convergence fails; in that case, draw another example
n_tries = 6;

progress = TextProgress('generating samples');
for i = 1:n_samples
    for j = 1:n_tries
        % generate two environments (the "contexts")
        Gamma1 = generate_environment('rnd_corr', size(S, 2), 'corr_beta', corr_beta);
        Gamma2 = generate_environment('rnd_corr', size(S, 2), 'corr_beta', corr_beta);
        
        % make sure matrices are symmetric
        Gamma1 = 0.5*(Gamma1 + Gamma1');
        Gamma2 = 0.5*(Gamma2 + Gamma2');
        
        % generate a perturbation
        % choose an odorant whose variance to increase
        idx_odorant = randi(size(Gamma1, 1));
        
        % choose a perturbation
        perturbation = 5*sqrt(median(Gamma1(:)))*randn(size(Gamma1, 1), 1);
        
        % increase variance for that particular odorant
        % chop away tiny complex component
        factor1 = real(sqrtm(Gamma1));
        factor1(:, idx_odorant) = factor1(:, idx_odorant) + perturbation;
        Gamma1_pert = factor1'*factor1;
        
        factor2 = real(sqrtm(Gamma2));
        factor2(:, idx_odorant) = factor1(:, idx_odorant) + perturbation;
        Gamma2_pert = factor2'*factor2;
        
        try
            % calculate optimal distributions
            K1 = calculate_optimal_dist(S, Gamma1, Ktot, optim_args{:});
            K2 = calculate_optimal_dist(S, Gamma2, Ktot, optim_args{:});
            
            K1_pert = calculate_optimal_dist(S, Gamma1_pert, Ktot, optim_args{:});
            K2_pert = calculate_optimal_dist(S, Gamma2_pert, Ktot, optim_args{:});
            
            worked = true;
        catch
            worked = false;
        end
        if worked
            break;
        end
    end
    if ~worked
        error('Convergence failed %d times at i = %d.', n_tries, i);
    end

    % store the results
    results.Gamma1{i} = Gamma1;
    results.Gamma2{i} = Gamma2;
    
    results.idx_odorant(i) = idx_odorant;
    results.Gamma1_pert{i} = Gamma1_pert;
    results.Gamma2_pert{i} = Gamma2_pert;
    
    results.K1(:, i) = K1;
    results.K2(:, i) = K2;
    
    results.K1_pert(:, i) = K1_pert;
    results.K2_pert(:, i) = K2_pert;
    
    % update progress bar
    progress.update(i/n_samples*100);
end
progress.done('done!');

%% Compare \Delta K in the two contexts

diff1 = results.K1_pert - K1;
diff2 = results.K2_pert - K2;
scatterfit(diff1(:), diff2(:), 'fitopts', {'showci', true, 'legend', 'fp'});
xlabel('\DeltaK_1');
ylabel('\DeltaK_2');

%% Save the results

save(fullfile('save', ['context_dependence' postfix '.mat']), 'S', 'S_choice', ...
    'S_size', 'artif_snr', 'artif_tuning', 'corr_beta', 'n_samples', ...
    'n_tries', 'optim_args', 'postfix', 'results');

%% SCRATCH

%% Generate two environments

% random environments, but make them reproducible
rng(32392076);

% * the concentrations used in Hallem&Carlson were high; here we're using
%   concentrations 1000 times smaller (so that the variance is about 1e-6)
% * the actual variances of odor concentrations in natural scenes are
%   unknown, so we effectively fit them to obtain neuron numbers in the
%   right ballpark
% Gamma1 = generate_environment('rnd_diag_const', size(S_fly, 2), ...
%     'diag_mu', 1e-6, 'diag_size', 1e-6, 'offdiag_mu', 1e-7);
% Gamma2 = generate_environment('rnd_diag_const', size(S_fly, 2), ...
%     'diag_mu', 1e-6, 'diag_size', 1e-6, 'offdiag_mu', 1e-7);
Gamma1 = generate_environment('rnd_corr', size(S, 2), 'corr_beta', corr_beta);
Gamma2 = generate_environment('rnd_corr', size(S, 2), 'corr_beta', corr_beta);

% make sure matrices are symmetric
Gamma1 = 0.5*(Gamma1 + Gamma1');
Gamma2 = 0.5*(Gamma2 + Gamma2');

%% Apply the same perturbation to both environments

% choose an odorant whose variance to increase
idx_odorant = 86;

% choose a perturbation
perturbation = 5*sqrt(median(Gamma1(:)))*randn(size(Gamma1, 1), 1);

% increase variance for that particular odorant
% chop away tiny complex component
factor1 = real(sqrtm(Gamma1));
factor1(:, idx_odorant) = factor1(:, idx_odorant) + perturbation;
Gamma1_pert = factor1'*factor1;

factor2 = real(sqrtm(Gamma2));
factor2(:, idx_odorant) = factor1(:, idx_odorant) + perturbation;
Gamma2_pert = factor2'*factor2;

%% Show the environments

cmap = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

fig = figure;
fig.Position(3:4) = 2*fig.Position(3:4);
fig.Color = [1 1 1];

cl = quantile(abs([Gamma1(:) ; Gamma2(:)]), 0.95);

ax = axes;
ax.OuterPosition = [0 1/2 1/2 1/2];
plotheat(Gamma1, [-cl, cl]);
colormap(cmap);
colorbar;
title('Environment matrix 1');

beautifygraph('minorticks', 'off');

ax = axes;
ax.OuterPosition = [1/2 1/2 1/2 1/2];
plotheat(Gamma1_pert, [-cl, cl]);
colormap(cmap);
colorbar;
title('Environment matrix 1, perturbed');

beautifygraph('minorticks', 'off');

ax = axes;
ax.OuterPosition = [0 0 1/2 1/2];
plotheat(Gamma2, [-cl, cl]);
colormap(cmap);
colorbar;
title('Environment matrix 2');

beautifygraph('minorticks', 'off');

ax = axes;
ax.OuterPosition = [1/2 0 1/2 1/2];
plotheat(Gamma2_pert, [-cl, cl]);
colormap(cmap);
colorbar;
title('Environment matrix 2, perturbed');

beautifygraph('minorticks', 'off');
preparegraph;

%% Calculate optimal receptor distributions

K1 = calculate_optimal_dist(S, Gamma1, Ktot, optim_args{:});
K2 = calculate_optimal_dist(S, Gamma2, Ktot, optim_args{:});

K1_pert = calculate_optimal_dist(S, Gamma1_pert, Ktot, optim_args{:});
K2_pert = calculate_optimal_dist(S, Gamma2_pert, Ktot, optim_args{:});
