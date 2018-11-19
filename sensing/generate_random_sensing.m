function [S, sigmas] = generate_random_sensing(M, N, tuning, snr, varargin)
% GENERATE_RANDOM_SENSING Generate a random sensing matrix.
%   S = GENERATE_RANDOM_SENSING(M, N, tuning, snr) generates an M x N random
%   sensing matrix with narrow (small `tuning` parameter) or wide (large
%   `tuning` parameter) tuning. Random noise is added as controlled by the
%   signal-to-noise parameter `snr`. The sensing matrix uses Gaussians
%   centered at random positions in odorant space (mapping to the columns of
%   the matrix), with periodic boundary conditions (i.e., treating the
%   Guassians as if they lived on a circle). The columns of the `S` matrix
%   are scrambled at the end. The peak of the Gaussians is 1, and the
%   random noise that is added is i.i.d. normal with a standard deviation
%   `1.0/snr`.
%
%   S = GENERATE_RANDOM_SENSING(M, N, [min_tuning, max_tuning], snr)
%   generates receptive fields with standard deviations chosen uniformly at
%   random from the range `[min_tuning, max_tuning]`.
%
%   [S, sigmas] = GENERATE_RANDOM_SENSING(...) also returns the standard
%   deviations `sigmas` for each receptor.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('shuffle', true, @(b) islogical(b) && isscalar(b));
parser.addParameter('noise', true, @(b) islogical(b) && isscalar(b));

% parse
parser.parse(varargin{:});
params = parser.Results;

S = zeros(M, N);
% start with Gaussian receptive fields with sizes in the given range,
% located randomly within odorant space
if isscalar(tuning)
    tuning = [tuning tuning];
end
sigmas = zeros(M, 1);
for i = 1:size(S, 1)
    crt_sigma = tuning(1) + diff(tuning)*rand;
    crt_pos = rand;
    % calculate distances as if the odorants were located on a circle (so
    % we're using periodic boundary conditions)
    crt_dist2 = 4*sin(pi*(linspace(0, 1, size(S, 2)) - crt_pos)) .^ 2;
    
    S(i, :) = exp(-0.5*crt_dist2/crt_sigma^2);
    sigmas(i) = crt_sigma;
end
% now scrable the columns to make this look more realistic
if params.shuffle
    S = S(:, randperm(size(S, 2)));
end
% and add noise
if params.noise
    S = S + randn(size(S))/snr;
end

end
