% Generate optimal receptor distributions for two environments in an
% Ibarra-Soria-like in-silico experiment.

%% Setup

% This script's behavior can be modified by defining some variables before
% running. When these variables are not defined, default values are used.
%
% Options:
%   Ktot0:
%       Average value for the total number of neurons.
%   n_perturbations:
%       Number of perturbations to consider.
%   pertsize_Gamma:
%   pertsize_S:
%   pertsize_Ktot:
%       Perturbation sizes for the parameters.

setdefault('Ktot0', 2000);
setdefault('n_perturbations', 24);
setdefault('pertsize_Gamma', 0.02);
setdefault('pertsize_S', 0);
setdefault('pertsize_Ktot', 0);

%% Load mouse sensing matrix

% load response curves
mouse_responses = open(fullfile('data', 'mouse_receptor_curves.mat'));
S_mouse = mouse_responses.sensing_matrix(:, :, 11);
% set missing data to 0
S_mouse(isnan(S_mouse)) = 0;
% remove receptors that don't respond to anything
mask = ~all(S_mouse == 0, 2);
S_mouse = S_mouse(mask, :);

%% Putting error bars on abundance change, mouse sensing matrix

% random environments, but make them reproducible
rng(335632);

optim_args = {'optimopts', ...
        {'MaxFunctionEvaluations', 50000, 'Display', 'notify-detailed'}, ...
    'method', 'lagsearch', 'sumtol', 0.01};

% setting the unperturbed values of the parameters

% using the mouse sensing matrix
S0 = S_mouse;
% ... with a little bit of noise to remove the degeneracy from the many
% elements of S that are exactly zero
% S0(S0 == 0) = 0.1*randn(1, sum(S0(:) == 0));

% generate random environments
Gamma1_0 = generate_environment('rnd_corr', ...
    size(S0, 2));
Gamma2_0 = generate_environment('rnd_corr', ...
    size(S0, 2));

% initialize storage for the perturbed parameters
cell_empty = cell(1, n_perturbations);
Gamma1 = cell_empty;
Gamma2 = cell_empty;

zero_empty = zeros(n_perturbations, 1);
Ktot1 = zero_empty;
Ktot2 = zero_empty;

S1 = cell_empty;
S2 = cell_empty;

K1 = zeros(size(S0, 1), ...
    n_perturbations);
K2 = zeros(size(S0, 1), ...
    n_perturbations);

Q1 = cell_empty;
Q2 = cell_empty;

% start calculating the perturbations and how they affect the optimal
% receptor distribution
progress = TextProgress('generating perturbations');
for i = 1:n_perturbations
    % perturb environments
    if pertsize_Gamma > 0
        crt_Gamma1 = generate_environment(...
            'delta_rnd_unif', Gamma1_0, ...
            'factor_size', pertsize_Gamma);
        crt_Gamma2 = generate_environment(...
            'delta_rnd_unif', Gamma2_0, ...
            'factor_size', pertsize_Gamma);
    else
        crt_Gamma1 = Gamma1_0;
        crt_Gamma2 = Gamma2_0;
    end
    
    Gamma1{i} = crt_Gamma1;
    Gamma2{i} = crt_Gamma2;
    
    % sanity check: symmetry is preserved
    if max(abs(flatten(crt_Gamma1 - crt_Gamma1'))) > 1e-6
        error('Perturbed Gamma1 not symmetric!');
    end
    if max(abs(flatten(crt_Gamma2 - crt_Gamma2'))) > 1e-6
        error('Perturbed Gamma2 not symmetric!');
    end
    
    % sanity check: positive-definiteness is preserved
    if any(eig(crt_Gamma1) < -1e-6)
        error('Perturbed Gamma1 not positive definite!');
    end
    if any(eig(crt_Gamma2) < -1e-6)
        error('Perturbed Gamma2 not positive definite!');
    end
        
    % perturb Ktot
    crt_Ktot1 = Ktot0*(1 + pertsize_Ktot*randn);
    crt_Ktot2 = Ktot0*(1 + pertsize_Ktot*randn);
    
    Ktot1(i) = crt_Ktot1;
    Ktot2(i) = crt_Ktot2;
    
    % sanity check: Ktots are still positive
    if crt_Ktot1 <= 0
        error('Perturbed Ktot1 <= 0!');
    end
    if crt_Ktot2 <= 0
        error('Perturbed Ktot2 <= 0!');
    end
    
    % perturb S
    crt_S1 = S0 + ...
        pertsize_S *...
        randn(size(S0));
    crt_S2 = S0 + ...
        pertsize_S *...
        randn(size(S0));
    
    S1{i} = crt_S1;
    S2{i} = crt_S2;
    
    % calculate optimal repertoires
    [crt_K1, ~, crt_Q1] = calculate_optimal_dist(crt_S1, crt_Gamma1, crt_Ktot1, ...
        optim_args{:}, 'lagstart', size(crt_S1, 1)/crt_Ktot1);
    [crt_K2, ~, crt_Q2] = calculate_optimal_dist(crt_S2, crt_Gamma2, crt_Ktot2, ...
        optim_args{:}, 'lagstart', size(crt_S2, 1)/crt_Ktot2);
    
    K1(:, i) = crt_K1;
    K2(:, i) = crt_K2;
    
    Q1{i} = crt_Q1;
    Q2{i} = crt_Q2;
    
    progress.update(100*i/n_perturbations);
end
progress.done('done!');

%% ... make Ibarra-Soria-like plot

K1_mean = mean(K1, 2);
K1_std = std(K1, [], 2);

K2_mean = mean(K2, 2);
K2_std = std(K2, [], 2);

logK_ratio = log2(K2_mean ./ K1_mean);
logK_std = (K1_std ./ K1_mean) + (K2_std ./ K2_mean);

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 2.9 2.3];

errorbar(K1_mean, logK_ratio, logK_std, logK_std, K1_std, K1_std, ...
    'marker', 'none', 'linestyle', 'none', 'color', [0.7 0.7 0.7], ...
    'capsize', 0);
hold on;
plot(K1_mean, logK_ratio, '.', 'markersize', 10);

% ylim([-3.5, 3.5]);
ylim([-2.5, 2.5]);
set(gca, 'xscale', 'log');
% xlim([5, 350]);
% xlim([1 max(K1p_mean)*1.2]);
xlim([1 100]);
h = plot(xlim, [0 0], 'r');
set(h, 'color', [1 0 0 0.3], 'linewidth', 1);

xlabel('K^{env1}');
ylabel('log_2 K^{env2}/K^{env1}');

beautifygraph('fontscale', 0.667, 'ticksize', 18);
% beautifygraph('fontscale', 0.667);

set(gca, 'box', 'on', 'xminortick', 'off', 'yminortick', 'off', 'TickLength', [0 0]);

preparegraph('edge', 0);

% safe_print(fullfile('figs', 'example_ibarra_soria_style.pdf'));

%% Save the data

save(fullfile('save', 'ibarra_soria_like.mat'), 'Ktot0', 'n_perturbations', ...
    'pertsize_Gamma', 'pertsize_S', 'pertsize_Ktot', 'Gamma1_0', 'Gamma2_0', ...
    'S0', 'Gamma1', 'Gamma2', 'K1', 'K2', 'Ktot1', 'Ktot2', 'Q1', 'Q2', ...
    'S1', 'S2', 'S_mouse', 'mouse_responses', 'optim_args');
