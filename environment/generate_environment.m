function [Gamma, details] = generate_environment(env_type, N, varargin)
% GENERATE_ENVIRONMENT Generate an artificial environment covariance matrix.
%   Gamma = GENERATE_ENVIRONMENT(env_type, N) generates an N x N environment
%   covariance matrix of the given type. Possible types are
%       * 'identity'    -- identity matrix
%       * 'rnd_diag'    -- random diagonal, using lognormal elements with
%                          parameters given by 'diagmu' and 'diagsize'
%                          options (see below)
%       * 'rnd_diag_const'
%                       -- random diagonal (as above) plus small, constant
%                          covariance (amount given by 'offdiagmu')
%       * 'rnd_diag_rnd'
%                       -- random diagonal (as above) plus small, random
%                          covariances (uniform distribution around
%                          'offdiagmu', range 2*'offdiagsize')
%       * 'rnd_product' -- as product M'*M for some random matrix M, where
%                          M has 'factorrows' rows, and is drawn for a
%                          normal distribution with standard deviation
%                          'factorsize'
%       * 'rnd_corr'    -- using randcorr function, with variances
%                          generated as in 'rnd_diag', and using 'corrbeta'
%                          as parameter for randcorr
%
%   Gamma = GENERATE_ENVIRONMENT(env_type, Gamma0) generates an environment
%   covariance matrix as a perturbation to a base matrix, `Gamma0`.
%   Possible types in this case are
%       * 'delta_rnd_diag'
%                       -- same as Gamma0 plus several independent odorants
%       * 'delta_rnd_prod'
%                       -- same as Gamma0 plus several co-occurring
%                          odorants
%       * 'delta_rnd_unif'
%                       -- decompose Gamma0 into M'*M for a matrix M, then
%                          add Gaussian noise of size 'factorsize' and
%                          regenerate product
%   If one of the types that don't require a base matrix is passed with
%   this form of the function, the argument `Gamma0` is silently ignored,
%   except for its size.
%
%   [Gamma, details] = GENERATE_ENVIRONMENT(...) also returns some
%   intermediate results used in generating the environment covariance
%   matrix. The kind of results that is returned depends on the options
%   used.
%
%   Options:
%    'diagmu'
%    'diagsize'
%       Mean and standard deviation for diagonal elements. The diagonal
%       elements are drawn from a lognormal distribution. Note that,
%       unlike the parameters to lognrnd, 'diagmu' and 'diagsize' are the
%       actual mean and standard deviation of the elements, not of the
%       underlying normal distribution.
%    'offdiagmu'
%    'offdiagsize'
%       Mean and size of fluctuations for off-diagonal elements. For
%       'rnd_diag_const', the size is ignored. For 'rnd_diag_rnd', the
%       off-diagonal elements are drawn from a uniform distribution, in the
%       range `offdiagmu - offdiagsize` to `offdiagmu + offdiagsize`.
%    'deltasize'
%       Size of perturbation added for the second form of the function.
%    'Ndelta'
%       Number of odorants added for the perturbation cases (the second
%       form of the function).
%    'deltapos'
%       Positions where to add odorants for the perturbation cases (the
%       second form of tee function).
%    'factorrows'
%       Number of rows used in the factor matrix `M` when generating the
%       covariance matrix as a product using the 'rnd_product' option.
%    'factorsize'
%       Standard deviation of elements in the factor matrix `M` (for the
%       'rnd_product' option).
%    'corrbeta'
%       Parameter to use with randcorr for the 'rnd_corr' option.

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('diagmu', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('diagsize', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('offdiagmu', 0.1, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('offdiagsize', 0.1, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('deltasize', 0.5, @(x) isnumeric(x) && isscalar(x));
parser.addParameter('Ndelta', 4, @(n) isnumeric(n) && isscalar(n) && n > 0);
parser.addParameter('deltapos', [], @(v) isnumeric(v) && isvector(v) && all(v > 0));
parser.addParameter('factorrows', @(N) 10*N, @(n) (isnumeric(n) && isscalar(n) && n > 0) ...
    || isa(n, 'function_handle'));
parser.addParameter('factorsize', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
parser.addParameter('corrbeta', 5, @(x) isnumeric(x) && isscalar(x) && x > 0);

% handle displaying defaults
if nargin == 1 && strcmp(env_type, 'defaults')
    disp('Available options and their defaults:');
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% identify the form of the function that's being used
if ~(isnumeric(N) && isscalar(N))
    if isnumeric(N) && ismatrix(N) && size(N, 1) == size(N, 2)
        Gamma0 = N;
        N = size(Gamma0, 1);
    else
        error([mfilename ':badN'], 'Second argument should be a positive integer or a covariance matrix.');
    end
end

% handle some special cases
if isa(params.factorrows, 'function_handle')
    params.factorrows = params.factorrows(N);
end

details = struct;

% handle each case
switch env_type
    case 'identity'
        Gamma = eye(N);
    case {'rnd_diag', 'rnd_diag_const', 'rnd_diag_rnd'}
        % start with the diagonal part
        ln_mu = log(params.diagmu^2 / sqrt(params.diagsize^2 + params.diagmu^2));
        ln_sigma = sqrt(log(params.diagsize^2 / params.diagmu^2 + 1));
        Gamma = diag(lognrnd(ln_mu, ln_sigma, N, 1));
        
        details.ln_mu = ln_mu;
        details.ln_sigma = ln_sigma;
        
        % now add off-diagonal elements, if required
        % XXX this also affects the diagonal!
        if strcmp(env_type, 'rnd_diag_const')
            Gamma = Gamma + params.offdiagmu;
        elseif strcmp(env_type, 'rnd_diag_rnd')
            off_diag = params.offdiagmu + 2*(params.offdiagsize - 0.5)*params.offdiagsize;
            off_diag = (off_diag + off_diag')/2;
            Gamma = Gamma + off_diag;
        end
    case 'rnd_product'
        M = params.factorsize*randn(params.factorrows, N);
        Gamma = M'*M;
        
        details.factor = M;
    case 'rnd_corr'
        % start with the diagonal part
        ln_mu = log(params.diagmu^2 / sqrt(params.diagsize^2 + params.diagmu^2));
        ln_sigma = sqrt(log(params.diagsize^2 / params.diagmu^2 + 1));
        variances = lognrnd(ln_mu, ln_sigma, N, 1);
        stdevs = sqrt(variances);
        
        Gamma = randcorr(N, params.corrbeta);
        Gamma = diag(stdevs)*Gamma*diag(stdevs);
        
        details.ln_mu = ln_mu;
        details.ln_sigma = ln_sigma;
        details.variances = variances;
    case 'delta_rnd_diag'
        % select some odorants
        if isempty(params.deltapos)
            idxs = randperm(N, params.Ndelta);
        else
            idxs = params.deltapos;
        end
        
        % add variance
        delta = zeros(size(Gamma0, 1), 1);
        delta(idxs) = params.deltasize;
        Gamma = Gamma0 + diag(delta);
        
        details.idxs = idxs;
        details.delta = delta;
    case 'delta_rnd_prod'
        % select some odorants
        if isempty(params.deltapos)
            idxs = randperm(N, params.Ndelta);
        else
            idxs = params.deltapos;
        end
        
        % add variance and covariance
        Gamma = Gamma0;
        Gamma(idxs, idxs) = Gamma(idxs, idxs) + params.deltasize;
        
        details.idxs = idxs;
    case 'delta_rnd_unif'
        % decompose Gamma
        M0 = sqrtm(Gamma0);
        M0 = M0 + params.factorsize*randn(size(M0));
        Gamma = M0'*M0;
        
        details.factor = M0;
    otherwise
        error([mfilename ':badtype'], 'Unknown environment type.');
end

% check that the matrix is symmetric
if max(abs(flatten(Gamma - Gamma'))) > 0
    if max(abs(flatten(Gamma - Gamma'))) > 1e-12
        error([mfilename ':notsym'], 'The matrix generated is not symmetric! This shouldn''t happen.');
    else
        % departure from symmetry is not serious, but make sure matrix is
        % exactly symmetric
        Gamma = 0.5*(Gamma + Gamma');
    end
end

% check that the matrix is positive-definite
evals = eig(Gamma);
if any(evals < 0)
    if any(evals < -1e-12)
        warning([mfilename ':notpos'], 'The matrix generated is not positive definite!');
    else
        % departure from positive-definiteness is not serious, but make
        % sure matrix is exactly positive-definite
        [V, D] = eig(Gamma);
        Gamma = V*max(D, 0)*V';
    end
end

end