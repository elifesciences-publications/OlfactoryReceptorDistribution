% Sweep over sensing matrices with receptors of various tuning widths and
% calculate optimal receptor distributions in two different environments
% for each choice of sensing matrix.

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   n_odorants:
%       Number of odorants in random environments.
%   n_receptors:
%       Number of receptor types.
%   sensing_snr:
%       Signal-to-noise ratio in artificial sensing matrices.
%   corr_beta:
%       Parameter to pass to the 'corr_beta' argument for generating random
%       environments.
%   Ktot:
%       Total number of neurons.
%   tuning_values:
%       Vector of receptor tuning values to use.
%   n_samples:
%       Number of samples to draw for each tuning width.
%   postifx:
%       String to add to file name in which the results are saved.

setdefault('n_odorants', 50);
setdefault('n_receptors', 24);
setdefault('corr_beta', 50);
setdefault('Ktot', 4000);
setdefault('tuning_values', logspace(log10(0.01), log10(0.2), 48));
setdefault('n_samples', 24);
setdefault('sensing_snr', 200);
setdefault('postfix', '');

%% Generate two random environments

% keep things reproducible
rng(135753);

Gamma1 = generate_environment('rnd_corr', n_odorants, 'corr_beta', corr_beta);
Gamma2 = generate_environment('rnd_corr', n_odorants, 'corr_beta', corr_beta);

%% Run the sweep

results = make_sweep_comparisons(...
    length(tuning_values), n_samples, ...
    @(i, k) generate_random_sensing(n_receptors, n_odorants, tuning_values(i), sensing_snr), ...
    @(i, k) Gamma1, @(i, k) Gamma2, @(i, k) Ktot, ...
    'optim_opts', {'optimopts', {'MaxFunctionEvaluations', 50000}, 'method', 'lagsearch'});

%% Save the results

save(fullfile('save', ['tuning_width_sweep' postfix '.mat']), 'n_odorants', 'n_receptors', ...
    'corr_beta', 'Ktot', 'tuning_values', 'n_samples', 'sensing_snr', ...
    'Gamma1', 'Gamma2', 'results');

%% Visualize some of the results

fig = figure;
fig.Position(3) = 2*fig.Position(3);

subplot(1, 2, 1);
semilogx(tuning_values, results.corrs.log_diagQ.cp(:, :, 1)', 'color', [1.0, 0.7, 0.7]);

xlim([min(tuning_values), max(tuning_values)]);
ylim([-0.1, 1]);

xlabel('Tuning width');
ylabel('corr(log diag(Q), K)');

subplot(1, 2, 2);
semilogx(tuning_values, -results.corrs.diag_invQ.cp(:, :, 1)', 'color', [1.0, 0.7, 0.7]);

xlim([min(tuning_values), max(tuning_values)]);
ylim([-0.1, 1]);

xlabel('Tuning width');
ylabel('corr(-diag(Q^{-1}), K)');
