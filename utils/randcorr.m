function S = randcorr(d, param, varargin)
% RANDCORR Generate random correlation matrix.
%   S = RANDCORR(d, param) generates a `d x d` random correlation matrix
%   (i.e., covariance matrix with diagonally identically equal to 1) by
%   drawing partial correlations from a beta distribution with parameters
%   `alpha = beta = param`. Large values of `param` lead to almost diagonal
%   matrices, while small values lead to strong off-diagonal correlations.
%
%   S = RANDCORR(d, 'uniform') samples uniformly from the space of `d x d`
%   random correlation matrices. Uniform in this sense means uniform in the
%   subspace of `R^(d*(d-1))` where positive-definite matrices live.
%
%   S = RANDCORR(d, 'uniform', eta) biases the sampling distribution by
%   a power of the determinant, `det(S)^eta`. XXX or eta-1
%
%   Code adapted from
%   https://stats.stackexchange.com/questions/2746/how-to-efficiently-generate-random-positive-semidefinite-correlation-matrices

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addOptional('eta', 0, @(x) isnumeric(x) && isscalar(x));

parser.addParameter('callback', [], @(f) isa(f, 'function_handle'));
parser.addParameter('figs', false, @(b) islogical(b) && isscalar(b));
parser.addParameter('shuffle', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

if strcmp(param, 'uniform')
    beta = params.eta + (d - 1)/2;
else
    beta = param;
end

P = zeros(d);           % storing partial correlations
S = eye(d);

for k = 1:d-1
    if strcmp(param, 'uniform')
        beta = beta - 1/2;
    end
    for i = k+1:d
        P(k,i) = betarnd(beta, beta); % sampling from beta
        P(k,i) = 2*(P(k,i) - 0.5);    % linearly shifting to [-1, 1]
        p = P(k,i);
        for l = (k-1):-1:1 %// converting partial correlation to raw correlation
            p = p * sqrt((1-P(l,i)^2)*(1-P(l,k)^2)) + P(l,i)*P(l,k);
        end
        S(k,i) = p;
        S(i,k) = p;
    end
    if params.figs
        clf;
        imagesc(S, [-1 1]);
        drawnow;
    end
    if ~isempty(params.callback) 
        params.callback(S);
    end
end

if ~strcmp(param, 'uniform') && params.shuffle
    % permuting the variables to make the distribution permutation-invariant
    permutation = randperm(d);
    S = S(permutation, permutation);
end

end