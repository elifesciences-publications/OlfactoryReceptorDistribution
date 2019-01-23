% Check our results in a small non-linear toy model.

%% Set some parameters

n_receptors = 3;
n_odorants = 15;

% this parameter allows to trade receptor numbers for per-receptor noise
noise_scaling = 10;
% noise_scaling = 1;

%% Choose a sensing matrix with few receptors

sensing_fly = open('data/flyResponsesWithNames.mat');
S_fly = sensing_fly.relRates';

% normalize by standard deviations of background rates
S_fly_normalized = bsxfun(@rdivide, S_fly, sensing_fly.bkgStd');

% keep things reproducible
rng(912);

% choose a subset of S_fly for our sensing matrix
receptor_idxs = randperm(size(S_fly_normalized, 1), n_receptors);
odorant_idxs = randperm(size(S_fly_normalized, 2), n_odorants);
% normalize to keep responses in a reasonable range below
% (note that we could instead of scaled Gamma up by a factor of 10,000 and
% left S unchanged)
% S = (1/200) * S_fly_normalized(receptor_idxs, odorant_idxs);
S = (1/100) * S_fly_normalized(receptor_idxs, odorant_idxs);

%% Generate two (small) environments

% keep things reproducible
rng(1351);

Gamma = generate_environment('rnd_corr', size(S, 2), ...
        'corr_beta', 2.0);
    
% generate two almost non-overlapping environments
factor = sqrtm(Gamma);
overlap = 1/5;
factor1 = factor; factor1(:, ceil(end/2):end) = factor1(:, ceil(end/2):end)*overlap;
factor2 = factor; factor2(:, 1:floor(end/2)) = factor2(:, 1:floor(end/2))*overlap;

Gamma1 = factor1'*factor1;
Gamma2 = factor2'*factor2;

% choose mean concentration
c0 = 0.5*ones(1, size(S, 2));

% build generators for concentration vectors in the two environments
c_generator1 = @(k) mvnrnd(c0, Gamma1, k)';
c_generator2 = @(k) mvnrnd(c0, Gamma2, k)';
    
%% Visualize sensing matrix and environments

cmap = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

fig = figure;
fig.Position(3) = 3*fig.Position(3);
fig.Color = [1 1 1];

ax = axes;
ax.OuterPosition = [0 0 1/3 1];
cl = quantile(abs(S(:)), 0.95);
plotheat(S, [-cl, cl]);
colormap(cmap);
colorbar('southoutside');
xlabel('Odorants');
ylabel('Receptors');
title('Sensing matrix');

beautifygraph('minorticks', 'off');

cl = quantile(abs([Gamma1(:) ; Gamma2(:)]), 0.95);

ax = axes;
ax.OuterPosition = [1/3 0 1/3 1];
plotheat(Gamma1, [-cl, cl]);
colormap(cmap);
colorbar;
title('Environment matrix 1');

beautifygraph('minorticks', 'off');

ax = axes;
ax.OuterPosition = [2/3 0 1/3 1];
plotheat(Gamma2, [-cl, cl]);
colormap(cmap);
colorbar;
title('Environment matrix 2');

beautifygraph('minorticks', 'off');
preparegraph;

%% Calculate analytic response covariance matrices for single receptor

% choose noise amounts for receptors
sigma = 0.1*noise_scaling*ones(1, size(S, 1));

Q1 = S*Gamma1*S';
Q2 = S*Gamma2*S';

cov_exact1 = Q1 + diag(sigma.^2);
cov_exact2 = Q2 + diag(sigma.^2);

%% Define the nonlinear response function

response_fct_nonlinear = @(c) S*c ./ (1 + S*c);

% The chosen response function ensures that the noiseless responses are
% between -1 and 1, or between 0 and 1 if c >= 0.
marginalize_opts = {'n_samples', 1e4, 'r_range', [-0.75 1.5], 'n_r_bins', 20};

%% Calculate info as function of K for both linear (exact) and nonlinear responses

% keep things reproducible
rng(183);

% the OSN numbers K1, K2, K3 are constrained by K1 + K2 + K3 = Ktot
% we thus have a 2d space of OSN numbers to sweep
Ktot = 2*noise_scaling^2;
% bin the responses for each receptor type
Kbins = linspace(0.05*noise_scaling^2, Ktot, 11);

% store info results in matrices indexed by K1 and K2
% not all values are valid (since sometimes K1 + K2 > Ktot); the invalid
% values will stay equal to NaN
info_map_nonlinear1 = nan(length(Kbins));
info_map_nonlinear2 = nan(length(Kbins));

info_map_exact1 = nan(length(Kbins));
info_map_exact2 = nan(length(Kbins));

% this is just a very naive way of counting n_total, the total number of
% bins that obey the constrain sum(K) = Ktot; this is only used to display
% progress information
n_total = 0;
for i = 1:length(Kbins)
    K1 = Kbins(i);
    for j = 1:length(Kbins)
        K2 = Kbins(j);
        K3 = Ktot - (K1 + K2);
        if K3 < 0
            break;
        end
        n_total = n_total + 1;
    end
end

ident1 = eye(size(Q1));
ident2 = eye(size(Q2));
crt = 1;
for i = 1:length(Kbins)
    K1 = Kbins(i);
    for j = 1:length(Kbins)
        K2 = Kbins(j);
        K3 = Ktot - (K1 + K2);
        if K3 < 0
            break;
        end
        
        disp(['Working on ' int2str(crt) ' / ' int2str(n_total) '...']);
        
        % effective noise amounts are decreased by factors of sqrt(K)
        % because of averaging
        sigma_eff = sigma(:) ./ sqrt([K1 ; K2 ; K3]);
        % estimate information values using the nonlinear response in both
        % environments
        info_map_nonlinear1(i, j) = calculate_info_nonlinear(response_fct_nonlinear, ...
            sigma_eff, c_generator1, 'marginalize_opts', marginalize_opts);
        info_map_nonlinear2(i, j) = calculate_info_nonlinear(response_fct_nonlinear, ...
            sigma_eff, c_generator2, 'marginalize_opts', marginalize_opts);
       
        % calculate information values exactly for the linear response
        info_map_exact1(i, j) = 0.5*log(det(ident1 + diag(sigma_eff.^(-2))*Q1));
        info_map_exact2(i, j) = 0.5*log(det(ident2 + diag(sigma_eff.^(-2))*Q2));
        
        crt = crt + 1;
    end
end

%% Compare the linear and nonlinear info landscapes

fig = figure;
fig.Position(3) = 2*fig.Position(3);
fig.Position(4) = 2*fig.Position(4);
fig.Color = [1 1 1];

cl_ex = quantile([info_map_exact1(:) ; info_map_exact2(:)], 0.95);
cl_nl = quantile([info_map_nonlinear1(:) ; info_map_nonlinear2(:)], 0.95);

ax = axes;
ax.OuterPosition = [0 0 0.5 0.5];

plotheat(info_map_nonlinear1, [0 cl_nl]);
title('Nonlinear numeric, env 1');

set(ax, 'xtick', 1:length(Kbins), 'xticklabel', Kbins);
set(ax, 'ytick', 1:length(Kbins), 'yticklabel', Kbins);
set(ax, 'ydir', 'normal');

beautifygraph('minorticks', 'off');

ax = axes;
ax.OuterPosition = [0.5 0 0.5 0.5];

plotheat(info_map_exact1, [0 cl_ex]);
title('Linear exact, env 1');

set(ax, 'xtick', 1:length(Kbins), 'xticklabel', Kbins);
set(ax, 'ytick', 1:length(Kbins), 'yticklabel', Kbins);
set(ax, 'ydir', 'normal');

beautifygraph('minorticks', 'off');

ax = axes;
ax.OuterPosition = [0 0.5 0.5 0.5];

plotheat(info_map_nonlinear2, [0 cl_nl]);
title('Nonlinear numeric, env 1');

set(ax, 'xtick', 1:length(Kbins), 'xticklabel', Kbins);
set(ax, 'ytick', 1:length(Kbins), 'yticklabel', Kbins);
set(ax, 'ydir', 'normal');

beautifygraph('minorticks', 'off');

ax = axes;
ax.OuterPosition = [0.5 0.5 0.5 0.5];

plotheat(info_map_exact2, [0 cl_ex]);
title('Linear exact, env 1');

set(ax, 'xtick', 1:length(Kbins), 'xticklabel', Kbins);
set(ax, 'ytick', 1:length(Kbins), 'yticklabel', Kbins);
set(ax, 'ydir', 'normal');

beautifygraph('minorticks', 'off');
preparegraph;

%% Save the results

save(fullfile('save', 'nonlinear_example.mat'), 'Gamma', 'Gamma1', 'Gamma2', ...
    'Ktot', 'Q1', 'Q2', 'S', 'sigma', 'cov_exact1', 'cov_exact2', ...
    'c_generator1', 'c_generator2', 'Kbins', 'c0', ...
    'info_map_nonlinear1', 'info_map_nonlinear2', ...
    'info_map_exact1', 'info_map_exact2', 'marginalize_opts', ...
    'response_fct_nonlinear', 'n_receptors', 'n_odorants');
