% Test dynamical model

%% Setup, and solve optimization problem directly

% load fly receptor affinity data
datafile = open('flyResponsesWithNames.mat');
data = datafile.relRates;

S = data';
[M, N] = size(S);

% generate a mock environment matrix
sigma = 10.0;

Gamma = 5e-6*ones(N) + 1.5e-4*eye(N);

% calculate optimal distributions at a range of Ktot
Ktot_values = [0 logspace(log10(0.1), log10(150), 50)];
[K, info_values, Q, info_fct] = calculate_optimal_dist(S/sigma, Gamma, Ktot_values);

% choose also a single Ktot where to do the calculations
Ktot0 = 30.0;
[K0, info_value0, Q0, info_fct0] = calculate_optimal_dist(S/sigma, Gamma, Ktot0);

%% Plot how the distribution depends on Ktot

figure;

% plot the information per OSN
K_scaling = 1;

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
preparegraph;

%% Plot the distribution at Ktot0

figure;

bar(K0);
xlabel('OSN receptor type');
ylabel('Receptor abundance');

beautifygraph;
preparegraph;

%% Run the dynamical model at Ktot0 and check convergence

nsteps = 5000;
Khistory = zeros(length(K0), nsteps);
lbd_history = zeros(1, nsteps);
info_history = zeros(1, nsteps);
info_normalized_history = zeros(1, nsteps);
%lbd_history = 0.1968*ones(1, nsteps);

Khistory(:, 1) = 1.0;

alpha = 1;
beta = 1e-4;
%beta = 0;

for i = 1:nsteps
    crtK = Khistory(:, i);
    crt_lbd = lbd_history(:, i);
    
    Gamma_R = Q0 + diag(1.0./crtK);
    Gamma_R_inv = inv(Gamma_R);
    
   
%   testing some other ways of calculating diag(inv(Gamma_R))

%     Gamma_R_inv_diag = zeros(1, size(Gamma_R, 1));
%     for k = 1:size(Gamma_R, 1)
%          crt_var = Gamma_R(k, k);
%          other_mask = ~ismember(1:size(Gamma_R, 1), k);
%          crt_covar = Gamma_R(other_mask, k);
%          crt_others_covar = Gamma_R(other_mask, other_mask);
%          
%         Gamma_R_inv_diag(k) = 1.0/crt_var + (crt_covar'*inv(...
%             crt_others_covar - crt_covar*crt_covar'/crt_var)*crt_covar)...
%             /crt_var^2;
%         Gamma_R_inv_diag(k) = 1.0/(crt_var - crt_covar'*(crt_others_covar\crt_covar));
%     end
%     if norm(Gamma_R_inv_diag - diag(Gamma_R_inv)') > 1e-12
%         disp('oops');
%     end
    
    derK = alpha*(crtK - crt_lbd*crtK.^2 - diag(Gamma_R_inv));
    der_lbd = beta*(sum(crtK) - Ktot0);
    
    Khistory(:, i+1) = crtK + derK;
    lbd_history(:, i+1) = crt_lbd + der_lbd;
    info_history(:, i) = info_fct0(crtK);
    
    crtK_normalized = crtK * Ktot0 / sum(crtK);
    info_normalized_history(:, i) = info_fct0(crtK_normalized);
    
    %disp([int2str(i) ' / ' int2str(nsteps) ' ::  eps_K: ' num2str(mean(abs(derK ./ crtK)), '%.2f') ...
    %    ', eps_lbd: ' num2str(abs(der_lbd / crt_lbd), '%.2f')]);
end

%% Look at two different environments

rng('default');
%odor1 = [1 0 zeros(1, N-2)];
%odor2 = [0 1 zeros(1, N-2)];
odor1 = zeros(1, N);
odor2 = zeros(1, N);
odor1(rand(N, 1) < 0.1) = 1;
odor2(rand(N, 1) < 0.1) = 1;
Gamma1 = Gamma + 0.1*diag(odor1);
Gamma2 = Gamma + 0.1*diag(odor2);

%Ktot_values = logspace(log10(0.1), log10(15000), 50);
optim_options = {'Algorithm', 'sqp', 'MaxFunctionEvaluations', 50000, ...
    'Display', 'notify-detailed', 'StepTolerance', 1e-14};
[K1, info_values_1, Q1, ~] = calculate_optimal_dist(S/sigma, Gamma1, Ktot0, ...
    optim_options{:});
[K2, info_values_2, Q2, ~] = calculate_optimal_dist(S/sigma, Gamma2, Ktot0, ...
    optim_options{:});

%% Check dynamical model predictions for the two environments

[K1dyn, history1] = run_dyn_model(Q1, Ktot0, 'tolsteps', 200);
[K2dyn, history2] = run_dyn_model(Q2, Ktot0, 'tolsteps', 200);

%% How do the results compare?

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 6 2];
%fig.Position = [fig.Position(1:2) 3 2];

subplot(1, 2, 1);

semilogx(K1 + eps, log2((K2 + eps) ./ (K1 + eps)), 'ko', ...
    'markerfacecolor', [0.5 0.5 0.5], 'markeredgecolor', 'none');
%yrange = 1.1*max(abs(ylim));
xlim([1e-3, 1e1]);
%ylim([-0.8 0.8]);
%ylim([-yrange yrange]);
%ylim([-2, 2]);

refline(0, 0);

xlabel('K^{env1}');
ylabel('log_2 K^{env2}/K^{env1}');
title('Direct solution');

grid on;
ax = gca;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

ax.FontSize = 8;

beautifygraph;

subplot(1, 2, 2);

semilogx(K1dyn + eps, log2((K2dyn + eps) ./ (K1dyn + eps)), 'ko', ...
    'markerfacecolor', [0.5 0.5 0.5], 'markeredgecolor', 'none');
%yrange = 1.1*max(abs(ylim));
xlim([1e-3, 1e1]);
%ylim([-0.8 0.8]);
%ylim([-yrange yrange]);
%ylim([-2, 2]);

refline(0, 0);

xlabel('K^{env1}');
ylabel('log_2 K^{env2}/K^{env1}');
title('Dynamical model solution');

grid on;
ax = gca;
ax.XMinorGrid = 'off';
ax.YMinorGrid = 'off';

ax.FontSize = 8;

beautifygraph;

preparegraph;