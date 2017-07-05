function [K, info_values, Q, info_fct] = calculate_optimal_dist(S, Gamma, Ktot_values, varargin)
% calculate_optimal_dist Calculate optimal receptor distribution given a
% sensing matrix and an environment covariance matrix.
%   K = CALCULATE_OPTIMAL_DIST(S, Gamma, Ktot_values) calculates the
%   optimal receptor distribution at a series of values for `Ktot`, the
%   total number of olfactory sensory neurons, given a normalized sensing
%   matrix `S` (the sensing matrix divided by the noise level) and a
%   covariance matrix for the odor environment `Gamma`. The optimal
%   abundances for each receptor type are returned as rows in the output
%   `K` matrix, with each column corresponding to a value of `Ktot`.
%
%   [K, info_values] = CALCULATE_OPTIMAL_DIST(...) also returns the maximum
%   information (in nats) for each optimal assignment.
%
%   [K, info_values, Q] = CALCULATE_OPTIMAL_DIST(...) also returns `Q`,
%   which is essentially the response covariance matrix.
%
%   [..., info_fct] = CALCULATE_OPTIMAL_DIST(...) also returns a function
%   calculating the mutual information given a certain receptor
%   distribution. That is, info_fct(Kvals) uses the value of Q calculated
%   using S and Gamma to estimate the amount of mutual information between
%   concentrations and responses when the distribution is given by Kvals.
%
%   Additional parameters passed to CALCULATE_OPTIMAL_DIST are passed as
%   options to the fmincon optimizer.
%
%   See also: fmincon, optimoptions.

[M, ~] = size(S);

Q = S * Gamma * S';

K = zeros(M, length(Ktot_values));

options = optimoptions(@fmincon, 'Display', 'none');
for k = 1:2:length(varargin)
    options.(varargin{k}) = varargin{k+1};
end
ident = eye(size(Q));
for i = 1:length(Ktot_values)
    crtK = Ktot_values(i);
    if crtK == 0
        continue;
    end
    options.TypicalX = ones(1, M)*(crtK/M);
    [Koptim, ~, exit_flag, output] = fmincon(...
        @(Khat) -0.5*log(det(ident + diag(Khat)*Q)), ...
        ones(1, M)*(crtK/M), ...
        [], [], ...
        ones(1, M), crtK, ...
        zeros(1, M), inf(1, M), ...
        [], ...
        options); %#ok<ASGLU>

    Koptim(Koptim < eps) = 0; % eliminate rounding errors
    K(:, i) = Koptim;
    
%    disp(['Step ' int2str(i)]);
%    disp(exit_flag);
%    disp(output);
end

% don't calculate the information values if they're never used
if nargin > 1
    info_values = zeros(1, length(Ktot_values));
    info_fct = @(Kvals) 0.5*log(det(ident + diag(Kvals)*Q));
    
    for i = 1:length(Ktot_values)
        if Ktot_values(i) > 0
            info_values(i) = info_fct(K(:, i));
        end
    end
end

end