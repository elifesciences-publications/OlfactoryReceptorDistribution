function [K, info, Pr, Pr_edges, details] = calculate_optimal_nonlinear(...
    response_fct, sigma, c_generator, Ktot, varargin)
% CALCULATE_OPTIMAL_NONLINEAR Calculate optimal receptor distribution given
% a nonlinear response function and an environment sampling function.
%   K = CALCULATE_OPTIMAL_NONLINEAR(response_fct, sigma, c_generator, Ktot)
%   calculates the optimal receptor distribution at the given value for the
%   total number of olfactory sensory neurons, `Ktot`. The optimization
%   assumes that `response_fct` gives the mapping between concentration
%   vectors and response vector; that the responses are subject to Gaussian
%   noise with standard deviations given by `sigma`; and that the
%   concentration vectors are drawn from a distribution that is sampled by
%   the `c_generator` function. Specifically, `c_generator(k)` returns a
%   matrix with `k` columns, each of which being a sample from the
%   environmental distribution of odorant concentrations.
%
%   [K, info] = CALCULATE_OPTIMAL_NONLINEAR(...) also returns the maximum
%   information (in nats) for the optimal assignment.
%
%   [K, info, Pr, Pr_edges] = CALCULATE_OPTIMAL_NONLINEAR(...) also returns
%   the numerical estimate for the response distribution at the optimal
%   receptor numbers, `Pr`. This is accompanied by the information
%   regarding the bin edges, `Pr_edges` (see `marginalize_concentrations`).
%
%   Options:
%    'method':
%       Optimization method. This can be
%         'fmincon'  -- use Matlab's `fmincon`
%         'gradient` -- use a simple gradient descent algorithm with
%                       projection onto the constraint at every step (by
%                       normalizing `sum(K)`)
%    'rate'
%       Learning rate for the 'gradient' method.
%    'grad_step'
%       Amount to add to `K` when calculating gradient.
%    'maxiter'
%       Maximum number of iterations for the 'gradient' method.
%    'infotol'
%       The optimization ends when the change in information drops below
%       this value.
%    'optimopts':
%       Cell array or structure containing options to be passed to fmincon
%       or fminunc for the optimization of the receptor abundances.
%    'sumtol'
%       Tolerance for value of sum(K)/Ktot in the optimal distribution. If
%       the difference between the sum of receptor numbers and the total
%       number Ktot is larger than this tolerance, an error is issued.
%    'marginalize_opts'
%       Options to pass to `marginalize_concentrations`.
%
%   See also: calculate_info_nonlinear, marginalize_concentrations.

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('method', 'fmincon', @(s) ismember(s, {'fmincon', 'gradient'}));
parser.addParameter('rate', 1e-1, @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('grad_step', 0.01, @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('maxiter', 50, @(i) isscalar(i) && isnumeric(i) && i > 0);
parser.addParameter('infotol', 1e-5, @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('optimopts', {}, @(c) (iscell(c) && isvector(c)) || isstruct(c));
parser.addParameter('sumtol', 1e-4, @(x) isscalar(x) && isnumeric(x) && x > 0);
parser.addParameter('marginalize_opts', {}, @(c) iscell(c) && isvector(c));

if strcmp(response_fct, 'defaults') && nargin == 1
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% handle options
optimopts = optimoptions('fmincon');
optimopts.Display = 'iter-detailed';
if iscell(params.optimopts)
    optimopts = optimoptions(optimopts, params.optimopts{:});
else
    optimopts = optimoptions(optimopts, params.optimopts);
end

% define the function to be optimized, making sure it saves last Pr result
marginalize_opts = params.marginalize_opts;
    function neg_info = get_neg_info(Khat)
        disp(['Evaluating at K = ' num2str(Khat(:)', '%10.3f')]);
        [info_val, Pr, Pr_edges] = calculate_info_nonlinear(...
            response_fct, sigma(:) ./ sqrt(Khat(:)), c_generator, ...
            'marginalize_opts', marginalize_opts);
        neg_info = -info_val;
    end

% perform the optimization
M = length(sigma);
if strcmp(params.method, 'fmincon')
    optimopts.TypicalX = ones(1, M)*(Ktot/M);
    [K, neg_info, exit_flag, output] = fmincon(...
        @get_neg_info, ...
        ones(1, M)*(Ktot/M), ...
        [], [], ...
        ones(1, M), Ktot, ...
        zeros(1, M), inf(1, M), ...
        [], ...
        optimopts); %#ok<ASGLU>
else
    % gradient method
    K = ones(1, M)*(Ktot/M);
    neg_info = get_neg_info(K);
    
    info_history = zeros(params.maxiter, 1);
    result = -1;
    for i = 1:params.maxiter
        disp(['Step ' int2str(i) ', info = ' num2str(-neg_info, '%.3f') '...']);
        % estimate gradient
        gradient = zeros(1, M);
        for k = 1:M
            K_pert = K;
            K_pert(k) = K_pert(k) + params.grad_step;
            neg_info_pert = get_neg_info(K_pert);
            gradient(k) = (neg_info_pert - neg_info) / params.grad_step;
        end
        
        % make gradient step
        dK = -params.rate*gradient;
        K = K + dK;
        
        % impose constraints
        K = max(K, 0);
        K = K * Ktot / sum(K);
        
        % recalculate neg info
        last_neg_info = neg_info;
        neg_info = get_neg_info(K);
        
        % collect convergence data
        info_history(i) = -neg_info;
        
        % check for convergence
        info_change = last_neg_info - neg_info;
        if abs(info_change) < params.infotol
            result = 0;
            disp('infotol reached.');
            break;
        end
    end
    if result == -1
        disp('maxiter reached.');
    end
    
    details.info_history = info_history(1:i);
end

K(K < eps) = 0; % eliminate rounding errors

info = -neg_info;

if abs(sum(K)/Ktot - 1) > params.sumtol
    error([mfilename ':badsum'], 'The total number of receptors did not converge.');
end

if nargout > 2
    % make sure the Pr we get is for the correct K
    get_neg_info(K);
end

end
