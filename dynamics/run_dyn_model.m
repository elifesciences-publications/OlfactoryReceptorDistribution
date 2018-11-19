function [Kstar, history] = run_dyn_model(Q, Ktot, varargin)
% RUN_DYN_MODEL Run dynamical model for finding optimal receptor
% distribution.
%   Kstar = RUN_DYN_MODEL(Q, Ktot) runs a dynamical model based on a
%   modified Verhulst dynamics to find the optimal receptor distribution
%   given the overlap matrix Q and the total number of receptors Ktot. It's
%   assumed that all responses have been normalized so that the noise for a
%   single receptor of any type has a standard deviation of 1.
%
%   [Kstar, history] = RUN_DYN_MODEL(Q, Ktot) also returns a history of the
%   convergence.
%
%   Options:
%    'maxsteps' <n>
%       Maximum number of steps to run in the simulation.
%       (default: 50000)
%    'tolinfo' <x>
%       Simulation ends when the change in information is lower than this
%       for `tolsteps` steps. The information is calculated using a
%       normalized version of the current Kstar such that the total number
%       of receptors is fixed at Ktot.
%       (default: 1e-6)
%    'tolsteps' <n>
%       Number of steps for which a tolerance criterion has to hold before
%       exiting.
%       (default: 50)
%    'tolK' <x>
%       Simulation ends when the maximum absolute value change in any
%       receptor abundance stays lower than this for at least `tolsteps`
%       steps.
%       (default: 1e-6)
%    'rate' <x>
%       Learning rate for receptor abundances. Larger values lead to faster
%       learning, but only as long as they are not too large. Instead, if
%       the learning rate is too large, the learning procedure fails.
%       (default: 1)
%    'ratelbd' <x>
%       Learning rate for the Lagrange multiplier that ensures that the
%       total number of receptor converges to Ktot. Large values for this
%       parameter lead to failed convergence.
%       (default: 0.1/Ktot^2)
%    'guessK' <v>
%       Starting point for receptor abundances.
%       (default: uniform -- ones(M, 1)/M, where M = size(Q, 1))
%    'guesslbd' <x>
%       Starting value for the Lagrange multiplier controlling the total
%       number of receptors.
%       (default: M/Ktot, where M = size(Q, 1))
%    'progressinterval' <x>
%       Output progress information every `x` seconds. Set to inf to
%       disable progress output.
%       (default: 10)

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

M = size(Q, 1);

parser.addParameter('maxsteps', 50000, @(x) isscalar(x) && isnumeric(x) && ...
    isreal(x) && x>0);
parser.addParameter('tolinfo', 1e-6, @(x) isscalar(x) && isnumeric(x) && ...
    isreal(x) && x>=0);
parser.addParameter('tolsteps', 50, @(x) isscalar(x) && isnumeric(x) && ...
    isreal(x) && x>0);
parser.addParameter('tolK', 1e-6, @(x) isscalar(x) && isnumeric(x) && ...
    isreal(x) && x>=0);
parser.addParameter('rate', 1, @(x) isscalar(x) && isnumeric(x) && ...
    isreal(x) && x>0);
parser.addParameter('ratelbd', 0.1/Ktot^2, @(x) isscalar(x) && isnumeric(x) && ...
    isreal(x) && x>=0);
parser.addParameter('guessK', ones(M, 1)/M, @(v) isvector(v) && isnumeric(v) && ...
    isreal(v) && all(v>=0));
parser.addParameter('guesslbd', M/Ktot, @(x) isscalar(x) && isnumeric(x) && ...
    isreal(x) && x>=0);
parser.addParameter('progressinterval', 10, @(x) isscalar(x) && isnumeric(x) && ...
    isreal(x) && x>=0);

% handle displaying defaults
if nargin == 1 && strcmp(Q, 'defaults')
    disp('Available options and their defaults:');
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% prepare the information function
info_fct = @(Kvals) 0.5*log(det(eye(M) + diag(Kvals)*Q));

% prepare the outputs
history.exitcode = -1;
history.exitreason = 'maxsteps';

history.K = zeros(M, params.maxsteps+1);
history.lbd = zeros(1, params.maxsteps+1);
history.info_naive = zeros(1, params.maxsteps+1);
history.info_normalized = zeros(1, params.maxsteps+1);

% initialize simulation
crtK = params.guessK;
crtlbd = params.guesslbd;

history.K(:, 1) = crtK;
history.lbd(1) = crtlbd;
history.info_naive(1) = info_fct(crtK);
history.info_normalized(1) = info_fct(crtK * Ktot / sum(crtK));

count_tolinfo = 0;
count_tolK = 0;

% run the simulation
t0 = tic;
displayed = false;
tlast = tic;
ending = false;
for i = 1:params.maxsteps
    % need the environment precision matrix (or at least its diagonal)
    Gamma_R = Q + diag(1.0./crtK);
    Gamma_R_inv = inv(Gamma_R);

    % calculate the derivatives using (modified) gradient descent
    derK = params.rate*(crtK - crtlbd*crtK.^2 - diag(Gamma_R_inv));
    derlbd = params.ratelbd*(sum(crtK) - Ktot);
    
    % update the state
    crtK = crtK + derK;
    crtlbd = crtlbd + derlbd;
    
    % update the history information
    history.K(:, i+1) = crtK;
    history.lbd(i+1) = crtlbd;
    history.info_naive(i+1) = info_fct(crtK);
    history.info_normalized(i+1) = info_fct(crtK * Ktot / sum(crtK));
    
    % check convergence
    if abs(history.info_normalized(i+1) - history.info_normalized(i)) < params.tolinfo
        count_tolinfo = count_tolinfo + 1;
        if count_tolinfo >= params.tolsteps
            history.exitcode = 0;
            history.exitreason = 'tolinfo';
            ending = true;
        end
    else
        count_tolinfo = 0;
    end
    
    if max(abs(derK)) < params.tolK
        count_tolK = count_tolK + 1;
        if count_tolK >= params.tolsteps
            history.exitcode = -2;
            history.exitreason = 'tolK';
            ending = true;
        end
    else
        count_tolK = 0;
    end
    
    if i == params.maxsteps
        ending = true;
    end
    
    % output progress information
    if toc(tlast) >= params.progressinterval || (displayed && ending)
        if ~displayed
            disp('   Step      info       lbd      derK     derlbd ');
        end
        disp([sprintf('%8d', i+1) '  ' sprintf('%8.3g', history.info_normalized(i+1)) ...
            '  ' sprintf('%8.3g', crtlbd) '  ' sprintf('%8.3g', max(abs(derK))) ...
            '  ' sprintf('%8.3g', abs(derlbd))]);
        displayed = true;
        tlast = tic;
    end
    
    if ending
        break;
    end
end

if i < params.maxsteps
    history.K = history.K(:, 1:i+1);
    history.lbd = history.lbd(1:i+1);
    history.info_naive = history.info_naive(1:i+1);
    history.info_normalized = history.info_normalized(1:i+1);
end

if displayed
    disp(' ');
    disp(['Simulation took ' num2str(toc(t0), '%.2f') ' seconds.']);
    disp(' ');
end

Kstar = crtK;

end