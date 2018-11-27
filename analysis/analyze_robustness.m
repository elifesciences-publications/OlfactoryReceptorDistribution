% Analyze how robust the optimal receptor distribution is to removing
% receptor or odorants.
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
%   n_samples:
%       Number of samples to generate for leave-one-out simulations.
%   n_samples_frac:
%       Number of samples to generate for each choice of the fraction of
%       removed receptors.
%   receptor_fractions:
%       Receptor fractions for which to analyze robustness.
%   corr_beta:
%       Parameter to use for generating correlated environment matrices.

%% Setup

setdefault('S_choice', 'fly');
setdefault('S_size', [24, 110]);
setdefault('artif_tuning', [0.2 0.8]);
setdefault('artif_snr', 200);
setdefault('Ktot', 100);
setdefault('n_samples', 100);
setdefault('n_samples_frac', 10);
setdefault('receptor_fractions', []);
setdefault('corr_beta', 0.3);

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

% generate colormap for covariance matrix plots
% cmap_covmat = divergent([0.21 0.17 0.53], [0.98 0.40 0.17], 256);
cmap_covmat = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

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

if isempty(receptor_fractions)
    if size(S, 1) < 25
        receptor_fractions = (0:size(S, 1)-4)/size(S, 1);
    else
        receptor_fractions = linspace(0, 1-4/size(S, 1), 20);
    end
end

S_size = size(S);

%% Find how robust the optimal distribution is to missing one receptor

% Suppose a receptor array has M receptors, but we only have sensing data
% for M-1, so that we can only predict the abundances using a reduced
% sensing matrix. Since there are correlations between all receptor
% responses, the abundances of the M-1 receptors that we do have data for
% will in general depend also on the response profile of the Mth receptor.
% However, we can ask how big the difference is. In other words, if we
% can't measure all receptors (which in fact is currently the case), can we
% still make progress measuring a subset?

% make this reproducible
rng(9102);

Gamma = cell(1, n_samples);
K0 = cell(1, n_samples);
Q0 = cell(1, n_samples);

idx_dropped = zeros(1, n_samples);

Kred = cell(1, n_samples);
Qred = cell(1, n_samples);

K0red = cell(1, n_samples);
Q0red = cell(1, n_samples);

progress = TextProgress('generating leave-one-out results');
n_tries = 6;
for i = 1:n_samples
    % choose a random environment
    Gamma{i} = 0.01*generate_environment('rnd_corr', size(S, 2), ...
        'corr_beta', corr_beta);
    
    % find full receptor distribution
    [K0{i}, ~, Q0{i}] = calculate_optimal_dist(S, Gamma{i}, Ktot, optim_args{:});
    
    % sometimes the optimization does not converge
    % try to repeat the simulation if it does not
    for k = 1:n_tries
        % drop one receptor and find new optimum
        idx_dropped(i) = randi(size(S, 1));
        Sred = S([1:idx_dropped(i)-1 idx_dropped(i)+1:end], :);
        Ktot_red = Ktot*(size(S, 1) - 1) / size(S, 1);
        try
            [Kred{i}, ~, Qred{i}] = calculate_optimal_dist(Sred, Gamma{i}, Ktot_red, optim_args{:});
            worked = true;
        catch
            worked = false;
        end
        
        if worked
            break;
        end
    end
    if ~worked
        error(['Failed convergence after ' int2str(n_tries) ' tries.']);
    end
    
    % calculate versions of K0 and Q0 with the one receptor removed
    K0red{i} = K0{i}([1:idx_dropped(i)-1 idx_dropped(i)+1:end]);
    Q0red{i} = Q0{i}([1:idx_dropped(i)-1 idx_dropped(i)+1:end], ...
                     [1:idx_dropped(i)-1 idx_dropped(i)+1:end]);
                 
    progress.update(100*i/n_samples);
end
progress.done('done!');

%% Plot effect of leaving one receptor out

fig = figure;
fig.Color = [1 1 1];
fig.Units = 'inches';
fig.Position(3:4) = [3 3];

% pool together the receptor abundances so that we can compare them
K0red_pooled = cell2mat(K0red');
Kred_pooled = cell2mat(Kred');

hold on;
smartscatter(K0red_pooled, Kred_pooled, 'color', [0.177 0.459 0.733], ...
    'alpha', 0.5, 'size', 10);
maxK = max([max(K0red_pooled) max(Kred_pooled)]);
plot([0 maxK], [0 maxK], 'color', [0.3 0.3 0.3], 'linewidth', 0.5);

axis equal;
xlim([0 inf]);
ylim([0 inf]);

xlabel('Full array optimization');
ylabel('Leave-one-out optimization');

beautifygraph('fontscale', 0.667, 'ticksize', 12, 'linewidth', 0.5);
preparegraph('edge', 0);

% safeprint(fullfile('figs', ['leave_one_out_effect' postfix]));

%% Find robustness to missing a fraction of receptors

% keep this reproducible
rng(9103);

corr_values = zeros(length(receptor_fractions), n_samples_frac);

progress = TextProgress('generating receptor fraction results', ...
    'suffix', ['frac=0 (0/' int2str(length(receptor_fractions)) ' sample 0/' ...
    int2str(n_samples_frac)]);
progress_total = n_samples_frac * length(receptor_fractions);
n_tries = 6;
for j = 1:length(receptor_fractions)
    f = receptor_fractions(j);
    n_drop = round(f*size(S, 1));
    
    for i = 1:n_samples_frac
        progress_count = (j-1)*n_samples_frac + (i-1);
        progress.update(100*progress_count/progress_total, ...
            'suffix', ['frac=' num2str(f, '%.2f') ...
            ' (' int2str(j) '/' int2str(length(receptor_fractions)) ') sample ' ...
            int2str(i) '/' int2str(n_samples_frac)]);
        
        % sometimes the optimization does not converge
        % try to repeat the simulation if it does not
        for k = 1:n_tries
            % choose a random environment
            crtGamma = 0.01*generate_environment('rnd_corr', size(S, 2), ...
                'corr_beta', corr_beta);
            
            % find full receptor distribution
            try
                crtK0 = calculate_optimal_dist(S, crtGamma, Ktot, optim_args{:});
                worked = true;
            catch
                worked = false;
            end
        
            if worked
                % drop a random set of receptors and find new optimum
                idxs_dropped = randperm(size(S, 1), n_drop);
                idxs_kept = setdiff(1:size(S, 1), idxs_dropped);
                Sred = S(idxs_kept, :);
                Ktot_red = Ktot*length(idxs_kept) / size(S, 1);
                try
                    crtKred = calculate_optimal_dist(Sred, crtGamma, Ktot_red, ...
                        optim_args{:}, 'lagrate', 2e-3);
                    worked = true;
                catch
                    worked = false;
                end
                if worked
                    break;
                end
            end
        end
        if ~worked
            error(['Failed convergence after ' int2str(n_tries) ' tries.']);
        end
        
        % calculate version of K0 with the receptors removed
        crtK0red = crtK0(idxs_kept);
        
        % calculate correlation
        corr_values(j, i) = corr(crtK0red, crtKred);
    end
end
progress.done('done!');

%% Plot results of removing a fraction of receptors

% find summary statistics for the correlations
corr_low = quantile(corr_values, 0.2, 2);
corr_high = quantile(corr_values, 0.8, 2);
corr_mean = nanmean(corr_values, 2);

% set some colors
color_uncertainty = [0.9 0.9 0.9];
color_mean = [0.737 0.180 0.172];

% make the figure
fig = figure;
fig.Units = 'inches';
fig.Color = [1 1 1];
fig.Position(3:4) = [3 2];

% draw the uncertainty area
actual_fractions = round(receptor_fractions*size(S, 1))/size(S, 1);
fill([flipud(actual_fractions(:)) ; actual_fractions(:)], ...
     [flipud(corr_low(:)) ; corr_high(:)], ...
     color_uncertainty, 'linestyle', 'none');
hold on;
% draw the mean
plot(actual_fractions, corr_mean, 'color', color_mean, 'linewidth', 1);

% adjust the axes limits
% xlim(tuning_min_max);
% ylim(fixed_y_lim);

% switch to log space on x axis and adjust the ticks
% set(ax, 'xtick', tuning_min_max, ...
%     'xticklabel', {[num2str(tuning_min_max(1)) ' (narrow)'], ...
%                    [num2str(tuning_min_max(2)) ' (wide)']}, ...
%     'xscale', 'log');

% label the axes
% xlabel('receptor tuning');
% ylabel('corr(log Q_{aa}, K_a)');

% plot the 0 correlation line as a visual guide
plot(xlim, [0 0], 'k:');
ylim([0 inf]);

% beautify, making sure fonts aren't too big, and axes don't waste ink
beautifygraph('fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5);

% make the ticks a bit bigger
% ax.TickLength = 2*ax.TickLength;

%% Save the results

save(fullfile('save', ['robustness' postfix]), ...
    'Gamma', 'K0', 'K0red', 'Kred', 'Ktot', 'Q0', 'Q0red', 'S', 'S_choice', ...
    'actual_fractions', 'artif_snr', 'artif_tuning', 'corr_beta', ...
    'corr_values', 'idx_dropped', 'n_samples', 'n_samples_frac', 'n_tries', ...
    'optim_args', 'receptor_fractions');
