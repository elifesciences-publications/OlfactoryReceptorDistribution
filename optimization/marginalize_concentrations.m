function [Pr, bin_edges] = marginalize_concentrations(...
    response_fct, sigma, c_generator, varargin)
% MARGINALIZE_CONCENTRATIONS Calculate marginal distribution over receptor
% responses in a nonlinear model.
%   Pr = MARGINALIZE_CONCENTRATIONS(response_fct, sigma, c_generator)
%   calculates the distribution over receptor responses by marginalizing
%   over the concentration distribution. `response_fct` is the response
%   function, taking in a concentration vector and returning the expected
%   value of the response vector. The responses are assumed to be normally-
%   distributed once conditioned on the concentration vector, with
%   independent errors having standard deviations given by the `sigma`
%   vector. The `c_generator` argument should be a callable that, when
%   called with an argument `k`, generates `k` sample concentration vectors
%   drawn from the concentration distribution (each as one column of the
%   output).
%
%   The return value, `Pr`, is a multi-dimensional array representing a
%   histogram of values for the probability distribution function of the
%   responses. The number of bins in each dimension and the range spanned
%   by the bins can be controlled -- see options.
%
%   The function works by generating a large number (see options) of
%   concentration vectors for each of which it accumulates the
%   corresponding contributions into the bins of `Pr`. As such, this can be
%   very time consuming.
%
%   [Pr, bin_edges] = MARGINALIZE_CONCENTRATIONS(...) also returns the
%   (lower) edges for the response bins in *each* direction.
%
%   Options:
%    'n_samples'
%       Number of concentration vectors to draw.
%    'n_r_bins'
%       Number of bins per *each* response direction.
%    'r_range'
%       Response range per *each* response direction.

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('n_samples', 1e6, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('n_r_bins', 100, @(x) isscalar(x) && isnumeric(x));
parser.addParameter('r_range', [-0.2 1.2], @(x) isnumeric(x) && ...
    isvector(x) && length(x) == 2);

if strcmp(response_fct, 'defaults') && nargin == 1
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% find the dimensionality of the response vector
n = length(sigma);

% create output array
size_opts = num2cell(repmat(params.n_r_bins, 1, n));
Pr = zeros(size_opts{:});

% find bin edges
bin_edges = linspace(params.r_range(1), params.r_range(2), params.n_r_bins+1);

% generate concentration samples
c_all = c_generator(params.n_samples);

% generate concentration vectors and accumulate Pr
progress = TextProgress('generating concentration samples');
for i = 1:params.n_samples
    % choose a concentration vector and find expected response
    r0 = response_fct(c_all(:, i));
    
    % build the response by multiplying the result for each direction
    subPr = ones(size_opts{:});
    for k = 1:n
        % this is the integral of the exponential part of the Gaussian
        % distribution, calculated over each bin (normalization below)
        binned_per_dir = -sigma(k)*sqrt(pi/2)*diff(...
            erf((r0(k) - bin_edges) / (sqrt(2)*sigma(k))));
        binned_per_dir_reshaped = reshape(binned_per_dir, [ones(1, k-1) ...
            params.n_r_bins ones(1, n-k)]);
        subPr = subPr .* binned_per_dir_reshaped;
    end
    
    % add to the distribution
    Pr = Pr + subPr;
    progress.update(100*i/params.n_samples);
end

% multiply by normalization factor
factor = prod(1 ./ (sqrt(2*pi)*sigma)) / params.n_samples;
Pr = Pr * factor;

progress.done('done!');

end
