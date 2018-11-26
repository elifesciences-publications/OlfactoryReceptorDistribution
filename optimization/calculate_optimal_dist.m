function [K, info_values, Q, info_fct] = calculate_optimal_dist(S, Gamma, Ktot_values, varargin)
% CALCULATE_OPTIMAL_DIST Calculate optimal receptor distribution given a
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
%   Options:
%    'method':
%       This can be
%           'fmincon':   use Matlab's fmincon to impose the constraint that
%                        sum(K) is fixed at Ktot.
%           'lagsearch': use a Lagrange multiplier to impose the constraint
%                        on sum(K); the value of the multiplier is fixed by
%                        a search algorithm
%    'optimopts':
%       Cell array or structure containing options to be passed to fmincon
%       or fminunc for the optimization of the receptor abundances.
%    'lagsearchopts':
%       Cell array or structure containing options to be passed to
%       fminunc for finding the value of the Lagrange multiplier.
%    'lagstart'
%       Initial value to use for Laagrange multiplier optimization.
%    'sumtol'
%       Tolerance for value of sum(K)/Ktot in the optimal distribution. If
%       the difference between the sum of receptor numbers and the total
%       number Ktot is larger than this tolerance, an error is issued.
%
%   See also: fmincon, optimoptions.

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;
%parser.KeepUnmatched = true;

parser.addParameter('method', 'fmincon', @(s) isvector(s) && ischar(s) && ...
    ismember(s, {'fmincon', 'lagsearch'}));
parser.addParameter('optimopts', {}, @(c) (iscell(c) && isvector(c)) || isstruct(c));
parser.addParameter('lagsearchopts', {}, @(c) (iscell(c) && isvector(c)) || ...
    isstruct(c));
parser.addParameter('lagstart', [], @(x) isscalar(x) && isnumeric(x));
parser.addParameter('sumtol', 1e-4, @(x) isscalar(x) && isnumeric(x) && x > 0);

if strcmp(S, 'defaults') && nargin == 1
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% setup
[M, ~] = size(S);
Q = S * Gamma * S';
K = zeros(M, length(Ktot_values));

% handle options
switch params.method
    case 'fmincon'
        optim_fct = @fmincon;
    case 'lagsearch'
        optim_fct = @fminunc;
end
optimopts = optimoptions(optim_fct);
optimopts.Display = 'none';
if strcmp(params.method, 'lagsearch')
    optimopts.SpecifyObjectiveGradient = true;
    optimopts.Algorithm = 'trust-region';
end
if iscell(params.optimopts)
    optimopts = optimoptions(optimopts, params.optimopts{:});
else
    optimopts = optimoptions(optimopts, params.optimopts);
end

lagsearchopts = optimoptions('fsolve');
lagsearchopts.Display = 'none';
if iscell(params.lagsearchopts)
    lagsearchopts = optimoptions(lagsearchopts, params.lagsearchopts{:});
else
    lagsearchopts = optimoptions(lagsearchopts, params.lagsearchopts);
end

% perform the optimization
ident = eye(size(Q));
if length(Ktot_values) > 1
    progress = TextProgress('varying Ktot...');
else
    progress = [];
end
for i = 1:length(Ktot_values)
    crtK = Ktot_values(i);
    if crtK == 0
        continue;
    end
    optimopts.TypicalX = ones(1, M)*(crtK/M);
    if strcmp(params.method, 'fmincon')
        [Koptim, ~, exit_flag, output] = fmincon(...
            @(Khat) -0.5*log(det(ident + diag(Khat)*Q)), ...
            ones(1, M)*(crtK/M), ...
            [], [], ...
            ones(1, M), crtK, ...
            zeros(1, M), inf(1, M), ...
            [], ...
            optimopts); %#ok<ASGLU>
    elseif strcmp(params.method, 'lagsearch')
        if isempty(params.lagstart)
            lagstart = M / crtK;
        else
            lagstart = params.lagstart;
        end
        lbd_optim = fsolve(@(lbd) ...
            sum(getKunc(lbd, M, crtK, Q, optimopts)) - crtK, lagstart, ...
            lagsearchopts);
        Koptim = getKunc(lbd_optim, M, crtK, Q, optimopts);
    end

    Koptim(Koptim < eps) = 0; % eliminate rounding errors
    K(:, i) = Koptim;
    
    if abs(sum(Koptim)/crtK - 1) > params.sumtol
        error([mfilename ':badsum'], 'The total number of receptors did not converge.');
    end
    
    if ~isempty(progress)
        progress.update(100*i/length(Ktot_values));
    end    
end
if ~isempty(progress)
    progress.done('done!');
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

function Kunc = getKunc(lbd, M, Ktot, Q, optimopts)

ident = eye(size(Q));

    function [I, dI] = neg_info_fct(sqrtKvals)
        Kvals = sqrtKvals(:).^2;
        shrunk = ident + diag(Kvals)*Q;
        I = 0.5*lbd*sum(Kvals) - 0.5*log(det(shrunk));
%         dI = Kvals.*(lbd - diag(Q/shrunk));
        dI = 2*sqrtKvals(:).*(lbd - diag(Q/shrunk));
    end

% info_fct = @(logKvals) 0.5*log(det(ident + diag(exp(logKvals))*Q)) ...
%     - 0.5*lbd*sum(exp(logKvals));
% [logKoptim, info_value, exit_flag, output] = fminunc(@neg_info_fct, ...
%     log(ones(1, M)*(Ktot/M)), optimopts); %#ok<ASGLU>
% Kunc = exp(logKoptim);
[sqrtKoptim, info_value, exit_flag, output] = fminunc(@neg_info_fct, ...
    sqrt(ones(1, M)*(Ktot/M)), optimopts); %#ok<ASGLU>
Kunc = sqrtKoptim.^2;

% disp(['lbd=' num2str(lbd) ', sum(K)=' num2str(sum(Kunc))]);

end
