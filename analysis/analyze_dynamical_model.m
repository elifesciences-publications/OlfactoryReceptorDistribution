% Briefly analyze the dynamical model for optimizing receptor abundances.
% NOTE: this uses results from the 'generate_ibarra_soria_like.m' script.
%       These are loaded from 'save/ibarra_soria_like.mat'.

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   ibarra_soria_fname:
%       File name from which to load the Ibarra-Soria-like simulation
%       results.
%   n_steps
%       Number of steps to run the dynamical model simulations for.

setdefault('ibarra_soria_fname', fullfile('save', 'ibarra_soria_like.mat'));
setdefault('n_steps', 5000);

%% Load the Ibarra-Soria-like results

load(ibarra_soria_fname);

%% Calculate initial and final receptor distributions using infomax

rng(3488247);

% choose sensing and environment matrices
S_infomax = S1{1};
Gamma1_infomax = Gamma1{1};
Gamma2_infomax = Gamma2{1};

% calculate optimal receptor distributions, overlap matrices, and mutual
% information functions for one of the samples from the error bars
% simulations
% (K and Q are recalculated; we're mostly just using this to get the mutual
% information function)
[K1_infomax, ~, Q1_infomax, info1] = calculate_optimal_dist(...
    S_infomax, Gamma1_infomax, Ktot0, optim_args{:}, ...
    'lagstart', size(S_infomax, 1)/Ktot0);
[K2_infomax, ~, Q2_infomax, info2] = calculate_optimal_dist(...
    S_infomax, Gamma2_infomax, Ktot0, optim_args{:}, ...
    'lagstart', size(S_infomax, 1)/Ktot0);

%% Run the dynamical model starting from random initial conditions

% random initial receptor abundances
Kini_random = rand(size(K1_infomax));
Kini_random = Kini_random*Ktot0/sum(Kini_random);
% run the population dynamics model in the second environment, starting
% from random distribution
% learning rate chosen so that exponential growth has doubling time 2
% guesslbd = 0.01759;
% guesslbd = 0.015756;
guesslbd = 0.015742;
[~, dyn_history_from_random] = run_dyn_model(Q2_infomax, Ktot0, ...
    'guessK', Kini_random, 'tolinfo', 0, 'tolK', 0, ...
    'maxsteps', n_steps, 'guesslbd', guesslbd, 'ratelbd', 1e-8, ...
    'rate', log(2)/2); % guesslbd used to be 0.01679

%% Run the dynamical model starting from optimum for environment 1

% initial abundances set to optimum for environment 1, *plus* random jitter
jitter_amt = 0.05*Ktot0/length(K1_infomax);
Kini_env1 = abs(K1_infomax + jitter_amt*(2*rand(size(K1_infomax)) - 1));
% make sure we're still normalized to the correct total number of neurons
Kini_env1 = Kini_env1*Ktot0/sum(Kini_env1);
[~, dyn_history_from_env1] = run_dyn_model(Q2_infomax, Ktot0, ...
    'guessK', Kini_env1, 'tolinfo', 0, 'tolK', 0, ...
    'maxsteps', n_steps, 'guesslbd', guesslbd, 'ratelbd', 1e-8, ...
    'rate', log(2)/2);

%% Make some plots

% plot the dynamical evolution of receptor abundances
yl = max([max(dyn_history_from_random.K(:)) max(dyn_history_from_env1.K(:))]);
    
col1 = [0.21 0.17 0.53];
col2 = [0.98 0.80 0.17];

fig = figure;
fig.Color = [1 1 1];
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
plot(n_steps*ones(size(K2_infomax)), K2_infomax, ...
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

beautifygraph('fontscale', 0.667);

% sanity check: how far are we from maximum info
disp('Random starting point');
disp(['Maximum info found with calculate_optimal_dist: ' num2str(info2(K2_infomax), '%.4f') ...
    ' run_dyn_model: ' num2str(info2(dyn_history_from_random.K(:, end)), '%.4f')]);
disp(['(Info at initial point: ' num2str(info2(Kini_random), '%.4f') ')']);
disp(' ');

% make plot starting with distribution optimized for environment 1
ax = axes;
ax.OuterPosition = [0.5 0 0.5 1];
hold on;
% plot optimal distribution in environment 1, as predicted from infomax
% (this is the same as initial conditions, up to the perturbation we added)
plot(ones(size(K1_infomax)), K1_infomax, ...
    'd', 'markerfacecolor', col1, 'markeredgecolor', 'none', ...
    'markersize', 3);
% plot final initial conditions as predicted from infomax
plot(n_steps*ones(size(K2_infomax)), K2_infomax, ...
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

beautifygraph('fontscale', 0.667);

% sanity check: how far are we from maximum info
disp('Starting from optimum for different environment');
disp(['Maximum info found with calculate_optimal_dist: ' num2str(info2(K2_infomax), '%.4f') ...
    ' run_dyn_model: ' num2str(info2(dyn_history_from_env1.K(:, end)), '%.4f')]);
disp(['(Info at initial point: ' num2str(info2(Kini_env1), '%.4f') ')']);

preparegraph('edge', 0);

% make sure we're saving vector graphics, not rasterized
fig.Renderer = 'painters';

% safeprint(fullfile('figs', 'dynamics_example.pdf'));

%% Look at convergence times

% trying out different definition for "convergence time"

% going 90% of the way from initial to final value
t_convergence_90p = zeros(size(K1_infomax));
% convergence to within 5% of final value
t_convergence_5p_error = zeros(size(K1_infomax));
% convergence to within 5% of Ktot
t_convergence_01p_Ktot = zeros(size(K1_infomax));

for i = 1:length(K1_infomax)
    crt_trajectory = dyn_history_from_random.K(i, :);
    
    deltaK = crt_trajectory(end) - crt_trajectory(1);
    [~, t_convergence_90p(i)] = min(abs(crt_trajectory - crt_trajectory(1) - 0.9*deltaK));
    
    [~, t_convergence_5p_error(i)] = min(abs(abs(crt_trajectory - crt_trajectory(end)) ...
        - 0.05*crt_trajectory(end)));
    
    [~, t_convergence_01p_Ktot(i)] = min(abs(abs(crt_trajectory - crt_trajectory(end)) ...
        - 0.001*Ktot0));
end

%% Save results

save(fullfile('save', 'dynmodel_results.mat'), 'Gamma1_infomax', ...
    'Gamma2_infomax', 'K1_infomax', 'K2_infomax', 'Kini_env1', ...
    'Kini_random', 'Ktot0', 'Q1_infomax', 'Q2_infomax', ...
    'S_infomax', 'dyn_history_from_env1', 'dyn_history_from_random', ...
    'guesslbd', 'info1', 'info2', 'jitter_amt', 'n_steps', 'optim_args', ...
    'ibarra_soria_fname');
