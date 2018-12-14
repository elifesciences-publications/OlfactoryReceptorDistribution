% This script calculates the optimal receptor distribution predicted by our
% model given receptor affinity data and measurements for the covariance
% structure of the olfactory environment.
%
% The script's behavior can be controlled by defining the variables
% described in the "Options" section below before running the script.
%
% The script's output is described in the "Outputs" section below. The
% script can be run one cell at a time to go through the results slowly, or
% it can be run it one go to output and save results.
%
% The script uses other functions defined in our framework. In order to
% access these, the script `setup_paths' should be run once in every Matlab
% session. Once the script is run, it does not need to be run again until
% Matlab is restarted.
%
% Options:
%   'sensing_file':
%       Set this variable to the name of a file containing sensing matrix
%       elements. This should be a Matlab .mat file. The 'sensing_field'
%       variable controls which variable in the .mat file is used for the
%       sensing matrix. Each row of this matrix should correspond to a
%       receptor type, and each column to an odorant. If this is not the
%       case, see the 'sensing_transpose' option.  Noise amounts for each
%       receptor type can be read from the same file, from the variable
%       indicated by the 'sensing_std_field' option. Otherwise they are
%       assumed to be equal to 1. The names of the receptors can also be
%       read from the file, from the variable indicated by the
%       'receptor_names_field' option.
%   'sensing_field':
%       Name of the variable in the 'sensing_file' file containing the
%       sensing matrix.
%   'sensing_std_field':
%       Name of the variable in the 'sensing_file' file containing the
%       standard deviations of the receptor noises for all receptor types.
%   'sensing_transpose'
%       Set this to `true` to transpose the sensing matrix before using.
%   'receptor_names_field':
%       Name of the variable in the 'sensing_file' file containing the
%       names of the receptors (corresponding to the rows of the sensing
%       matrix) as a cell array.
%   'odorant_names_field':
%       Name of the variable in the 'sensing_file' file containing the
%       names of the odorants (corresponding to the columns of the sensing
%       matrix) as a cell array.
%   'receptor_names'
%       Cell array of receptor names. This overrides the
%       'receptor_names_field' option.
%   'odorant_names'
%       Cell array of odorant names. This overrides the
%       'odorant_names_field' option.
%   'sensing_matrix':
%       Set this variable to directly choose the sensing matrix. This
%       overrides the 'sensing_file' option.
%   'sensing_std':
%       Set this variable to directly choose the standard deviations for
%       the receptor noises. This overrides the 'sensing_std_field' option.
%   'environment_file':
%       Set this to the name of a file containing an environment covariance
%       matrix. This should be a Matlab .mat file. The 'environment_field'
%       variable controls which variable in the .mat file is used for the
%       covariance matrix. The 'odorant_names_field' variable controls
%       which variable in the .mat file gives the names of the odorants.
%   'environment_matrix'
%       Set this variable to directly choose the environment matrix. This
%       overrides the 'environment_file' option.
%   'Ktot'
%       Total number of neurons in the epithelium.
%   'optim_args'
%       These are parameters passed to `calculate_optimal_dist` to control
%       the convergence of the optimization algorithm. Consider changing
%       these if there are convergence problems.
%   'out_file'
%       Set to an output file name to save the outputs to file. This also
%       saves the parameters used for the calculation (the sensing matrix
%       S, the environment matrix Gamma, noise amounts sigma, etc.)
%
% Outputs:
%   'K'
%       The optimal OSN abundances calculated given the parameters.
%
% Example:
%
%   To use the fly receptor affinities from Hallem&Carlson (2006) included
%   with this package, together with a randomly-generated environment, run:
%
%      sensing_file = fullfile('data', 'flyResponsesWithNames.mat');
%      sensing_field = 'relRates';
%      sensing_transpose = true;  % this file contains receptors along columns
%      sensing_std_field = 'bkgStd';
%      receptor_names_field = 'orNames';
%      odorant_names_field = 'odorNames';
%      environment_matrix = 1e-4*generate_environment('rnd_corr', 110);
%      get_predictions_from_data

%% Set default values

setdefault('sensing_file', '');
setdefault('sensing_field', 'sensing_matrix');
setdefault('sensing_std_field', '');
setdefault('sensing_transpose', false);
setdefault('receptor_names_field', '');
setdefault('odorant_names_field', '');
setdefault('receptor_names', {});
setdefault('odorant_names', {});
setdefault('sensing_matrix', []);
setdefault('sensing_std', []);
setdefault('environment_file', '');
setdefault('environment_matrix', []);
setdefault('Ktot', 1000);
setdefault('optim_args', {'optimopts', ...
        {'MaxFunctionEvaluations', 50000, 'Display', 'notify-detailed'}, ...
        'method', 'lagsearch'});
setdefault('out_file', '');

%% Load/set up the sensing matrix

% need either an input file or a direct sensing matrix
if isempty(sensing_file) && isempty(sensing_matrix)
    error('Either sensing_file or sensing_matrix must be defined before running this script.');
end

if ~isempty(sensing_matrix)
    % use the provided data
    S = sensing_matrix;
    sigma = sensing_std;
    ors = receptor_names;
    odorants = odorant_names;
else
    % load the data from file
    sensing_data = open(sensing_file);
    S = sensing_data.(sensing_field);
    
    if ~isempty(sensing_std_field)
        if ~isfield(sensing_data, sensing_std_field)
            error(['The selected sensing_std_field not found in ' sensing_file '.']);
        end
        sigma = sensing_data.(sensing_std_field);
    else
        sigma = [];
    end
    
    if ~isempty(receptor_names_field)
        if ~isfield(sensing_data, receptor_names_field)
            error(['The selecting receptor_names_field not found in ' sensing_file '.']);
        end
        ors = sensing_data.(receptor_names_field);
    else
        ors = [];
    end
    
    if ~isempty(odorant_names_field)
        if ~isfield(sensing_data, odorant_names_field)
            error(['The selecting odorant_names_field not found in ' sensing_file '.']);
        end
        odorants = sensing_data.(odorant_names_field);
    else
        odorants = [];
    end
end

% normalize vector dimensions
sigma = sigma(:);
ors = ors(:)';
odorants = odorants(:)';

% transpose the sensing matrix if necessary
if sensing_transpose
    S = S';
end

% get the number of receptors (M) and of odorants (N)
[M, N] = size(S);

% set default noise to 1 for all receptor types
if isempty(sigma)
    sigma = ones(M, 1);
end

% set default names for receptors and odorants
if isempty(ors)
    ors = arrayfun(@(i) ['R' int2str(i)], 1:M, 'uniform', false);
end
if isempty(odorants)
    odorants = arrayfun(@(i) ['odor' int2str(i)], 1:N, 'uniform', false);
end

%% Load/set up the environment covariance

% need either an input file or a direct environment covariance matrix
if isempty(environment_file) && isempty(environment_matrix)
    error('Either environment_file or environment_matrix must be defined before running this script.');
end

if ~isempty(environment_matrix)
    Gamma = environment_matrix;
else
    environment_data = open(environment_file);
    Gamma = environment_data.(environment_field);
end

% make sure the covariance matrix is square and symmetric
if size(Gamma, 1) ~= size(Gamma, 2) || max(abs(flatten(Gamma - Gamma'))) > 1e-6
    error('The environment covariance matrix should be square and symmetric.');
end

% make sure that the environment size matches the sensing matrix
if size(Gamma, 1) ~= N
    error('The size of the environment matrix does not match the number of columns in the sensing matrix.');
end

%% Calculate optimal receptor distribution

K = calculate_optimal_dist(S ./ sigma, Gamma, Ktot, optim_args{:});

%% Save results

if ~isempty(out_file)
    save(out_file, 'K', 'S', 'Gamma', 'sigma', 'Ktot', 'ors', 'odorants');
    disp(['Outputs saved to ' out_file '.']);
end

%% Make a plot of the optimal distribution

fig = figure;
fig.Position(3) = 2*fig.Position(3);
fig.Color = [1 1 1];

ax = axes;
bar(K, 'edgecolor', 'none');
ax.XTick = 1:M;
ax.XTickLabel = ors;

% xlabel('Receptor');
ylabel('Optimal OSN abundance');

beautifygraph('linewidth', 0.5, 'minorticks', 'off');
preparegraph;

%% Display the results as a table

abundance_table = table;
abundance_table.ORName = ors(:);

% clip very small abundance values to 0
Kclip = K;
Kclip(Kclip < Ktot*1e-6) = 0;

abundance_table.OptimalAbundance = Kclip;

disp(abundance_table);
