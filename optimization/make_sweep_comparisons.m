function results = make_sweep_comparisons(n, n_samples, S_fct, ...
    Gamma1_fct, Gamma2_fct, Ktot_fct, varargin)
% MAKE_SWEEP_COMPARISONS Calculate optimal receptor distribution at a
% sequence of parameters, together with some summary statistics.
%   results = MAKE_SWEEP_COMPARISONS(...
%       n, n_samples, S_fct, Gamma1_fct, Gamma2_fct, Ktot_fct)
%   calculates the optimal receptor distribution for `n` different sets of
%   parameters, running the optimization `n_samples` times for each set
%   (this is useful when at least some of the parameters depend on the
%   sample index `k`), each time using both the environment given by
%   `Gamma1_fct` and the one from `Gamma2_fct`. Setting `Gamma2_fct` to an
%   empty matrix leads to calculations only being performed for the first
%   environment.
%
%   The parameters are obtained by calling the `S_fct` (for the sensing
%   matrix), `Gamma1_fct` and `Gamma2_fct` (for the environment matrices)
%   and `Ktot_fct` (for the total neuron number) functions with two
%   arguments: `i`, representing the index of the current parameter set;
%   and `k`, representing the current sample index. Each of these parameters
%   can alternatively be a cell array, in which case the `i`th element is
%   used for all `k` samples.
%
%   The output, `results`, is a structure containing the following fields
%    'inputs'
%       A structure containing the inputs that were used. These are arrays
%       with `n_samples` rows and `n` columns. `Gamma` has a third dimension
%       with 2 elements unless `Gamma2_fct` is empty.
%        'S'
%        'Gamma'
%        'Ktot'
%           Arrays of parameters returned from `S_fct`, `Gamma1_fct`, 
%           `Gamma2_fct`,  and `Ktot_fct`. `S` and `Gamma` are cell arrays,
%           `Ktot` is a numeric array.
%        'optim_opts'
%           Cell array of options passed to `calculate_optimal_dist`.
%    'outputs'
%       A structure containing raw outputs from `calculate_optimal_dist`.
%       These are arrays with `n_samples` rows, `n` columns, and 2-element
%       third dimensions (unless `Gamma2_fct` is empty).
%        'K'
%           Cell array of receptor distributions.
%        'info_values'
%           Numeric array of maximum information values.
%        'Q'
%           Cell array of overlap matrices.
%        'info_fct'
%           Cell array of mutual information function returned by
%           `calculate_optimal_dist`.
%    'corrs'
%       A structure containing various correlations, each presented as an
%       array with `n` columns,  `n_samples` rows, and potentially a third
%       dimension. The third dimensions is missing when `Gamma2_fct` is
%       empty. The fields starting with `diff_` always lack a third
%       dimension, and only exist if `Gamma2_fct` is not empty.
%        'diagQ.cs'
%        'diagQ.cp'
%           Spearman and Pearson correlations between optimal receptor
%           distributions and the diagonal elements of the overlap matrix Q.
%        'log_diagQ.cs'
%        'log_diagQ.cp'
%           Spearman and Pearson correlations between optimal receptor
%           distributions and the logarithms of the diagonal elements of the
%           overlap matrix Q.
%        'diag_invQ.cs'
%        'diag_invQ.cp'
%           Spearman and Pearson correlations between optimal receptor
%           distributions and the diagonal elements of the inverse overlap
%           matrix A = inv(Q).
%        'diag_log_invQ.cs'
%        'diag_log_invQ.cp'
%           Spearman and Pearson correlations between optimal receptor
%           distributions and the logarithms of the diagonal elements of the
%           inverse overlap matrix A = inv(Q).
%        'diff_diagQ.cs' (only with two sets of environment matrices)
%        'diff_diagQ.cp'
%           Spearman and Pearson correlations between the change in optimal
%           receptor distributions between environments, and the change in
%           the diagonal elements of the overlap matrix Q.
%        'diff_log_diagQ.cs' (only with two sets of environment matrices)
%        'diff_log_diagQ.cp'
%           Spearman and Pearson correlations between the change in optimal
%           receptor distributions between environments, and the change in
%           the logarithms of the diagonal elements of the overlap matrix Q.
%        'diff_diag_invQ.cs' (only with two sets of environment matrices)
%        'diff_diag_invQ.cp'
%           Spearman and Pearson correlations between the change in optimal
%           receptor distributions between environments, and the change in
%           the diagonal elements of the inverse overlap matrix A = inv(Q).
%        'diff_diag_log_invQ.cs' (only with two sets of environment matrices)
%        'diff_diag_log_invQ.cp'
%           Spearman and Pearson correlations between the change in optimal
%           receptor distributions between environments, and the change in
%           the logarithms of the diagonal elements of the inverse overlap
%           matrix A = inv(Q).
%    'summaries'
%       Summarized versions of the correlations in 'corrs' obtained by
%       reducing over the first (sample) axis of the structure above using
%       either the mean, the median, the standard deviation, or the 0.2 or
%       0.8 quantiles. These are indicated by the suffixes '_mean',
%       '_median', '_std', '_low', and '_high', respectively.
%
%   Options:
%    'optim_opts'
%       Options passed to `calculate_optimal_dist`. This can be a simple
%       cell array, a cell array of cell arrays, or a function taking one
%       parameter (the index `i` of the current optimization). The same
%       options are always used for all the samples, and also for both
%       environment matrices and two are used.
%
%   See also: calculate_optimal_dist.

parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('optim_opts', []);

if strcmp(n, 'defaults') && nargin == 1
    parser.parse;
    disp(parser.Results);
    return;
end

% parse
parser.parse(varargin{:});
params = parser.Results;

% handle function vs. cell array inputs
if iscell(S_fct)
    S_fct_cell = S_fct;
    S_fct = @(i, k) S_fct_cell{i};
end
if iscell(Ktot_fct)
    Ktot_fct_cell = Ktot_fct;
    Ktot_fct = @(i, k) Ktot_fct_cell{i};
elseif isvector(Ktot_fct)
    Ktot_fct_vec = Ktot_fct;
    Ktot_fct = @(i, k) Ktot_fct_vec(i);
end
if iscell(Gamma1_fct)
    Gamma1_fct_cell = Gamma1_fct;
    Gamma1_fct = @(i, k) Gamma1_fct_cell{i};
end
if ~isempty(Gamma2_fct)
    if iscell(Gamma2_fct)
        Gamma2_fct_cell = Gamma2_fct;
        Gamma2_fct = @(i, k) Gamma2_fct_cell{i};
    end
end
if iscell(params.optim_opts)
    % handle single cell array input
    if isempty(params.optim_opts) || ~iscell(params.optim_opts{1})
        params.optim_opts = repmat({params.optim_opts}, 1, n);
    end
    
    optim_opts_cell = params.optim_opts;
    params.optim_opts = @(i, k) optim_opts_cell{i};
end

% how many Gammas do we have?
n_gamma = 2 - isempty(Gamma2_fct);

% define some convenient empty arrays
cell_empty = cell(n_samples, n);
cell_empty3 = cell(n_samples, n, n_gamma);
num_empty = zeros(n_samples, n);
num_empty3 = zeros(n_samples, n, n_gamma);

% build results.inputs structure
results.inputs.S = cell_empty;
results.inputs.Gamma = cell_empty3;
results.inputs.Ktot = num_empty;

% build results.outputs structure
results.outputs.K = cell_empty3;
results.outputs.info_values = num_empty3;
results.outputs.Q = cell_empty3;
results.outputs.info_fct = cell_empty3;

% start sweeping!
gen_progress_text = @(i, k) ...
    sprintf('sweeping %d.%d / %d.%d...', i, k, n, n_samples);
progress = TextProgress(gen_progress_text(0, 0), 'prespace', 32);
for i = 1:n
    for k = 1:n_samples
        % get the current parameters
        crt_S = S_fct(i, k);
        crt_Gamma1 = Gamma1_fct(i, k);
        if ~isempty(Gamma2_fct)
            crt_Gamma2 = Gamma2_fct(i, k);
        else
            crt_Gamma2 = [];
        end
        crt_Ktot = Ktot_fct(i, k);
        
        % ...and the current optimization options
        crt_opts = params.optim_opts(i, k);
        
        % calculate the optimal distributions
        [results.outputs.K{k, i, 1}, results.info_values(k, i, 1), ...
         results.outputs.Q{k, i, 1}, results.info_fct{k, i, 1}] = ...
            calculate_optimal_dist(crt_S, crt_Gamma1, crt_Ktot, crt_opts{:});
        if ~isempty(crt_Gamma2)
            [results.outputs.K{k, i, 2}, results.info_values(k, i, 2), ...
             results.outputs.Q{k, i, 2}, results.info_fct{k, i, 2}] = ...
                calculate_optimal_dist(crt_S, crt_Gamma2, crt_Ktot, crt_opts{:});
        end

        % store copies of the inputs
        results.inputs.S{k, i} = crt_S;
        results.inputs.Gamma{k, i, 1} = crt_Gamma1;
        if ~isempty(crt_Gamma2)
            results.inputs.Gamma{k, i, 2} = crt_Gamma2;
        end
        results.inputs.Ktot(k, i) = crt_Ktot;
        
        progress.update(100*((i-1)*n_samples + k)/ (n*n_samples), ...
            'prefix', gen_progress_text(i, k));
    end
end
progress.done('done');

% decide what we will correlate with the predicted Ks
abs_proxies = {'diagQ', @(Q) diag(Q), 'log_diagQ', @(Q) log(diag(Q)), ...
    'diag_invQ', @(Q) diag(inv(Q)), 'log_diag_invQ', @(Q) log(diag(inv(Q)))};
% decide what we will correlate with the differences K1 - K2 between the
% two environments
diff_proxies = {'diff_diagQ', @(Q1, Q2) diag(Q2) - diag(Q1), ...
    'diff_log_diagQ', @(Q1, Q2) log(diag(Q2)) - log(diag(Q1)), ...
    'diff_diag_invQ', @(Q1, Q2) diag(inv(Q2)) - diag(inv(Q1)), ...
    'diff_log_diag_invQ', @(Q1, Q2) log(diag(inv(Q2))) - log(diag(inv(Q1)))};
% decide what kinds of correlations we will calculates
corr_fcts = {'cs', @(x, y) corr(x, y, 'type', 'spearman'), ...
    'cp', @(x, y) corr(x, y, 'type', 'pearson')};

% build results.corrs structure, correlations for each environment
results.corrs = struct;
for i = 1:length(abs_proxies)/2
    crt_proxy_name = abs_proxies{2*i-1};
    
    results.corrs.(crt_proxy_name) = struct;
    for j = 1:length(corr_fcts)/2
        crt_corr_name = corr_fcts{2*j-1};
        results.corrs.(crt_proxy_name).(crt_corr_name) = num_empty3;
    end
end
if ~isempty(Gamma2_fct)
    % build results.corrs structure, correlations between environments
    for i = 1:length(diff_proxies)/2
        crt_proxy_name = diff_proxies{2*i-1};
        
        results.corrs.(crt_proxy_name) = struct;
        for j = 1:length(corr_fcts)/2
            crt_corr_name = corr_fcts{2*j-1};
            results.corrs.(crt_proxy_name).(crt_corr_name) = num_empty;
        end
    end
end

% calculate the actual correlations
for i = 1:n
    for k = 1:n_samples
        % collect the Qs and Ks
        crt_Qs = results.outputs.Q(k, i, 1);
        crt_Ks = results.outputs.K(k, i, 1);
        if ~isempty(Gamma2_fct)
            crt_Qs{2} = results.outputs.Q{k, i, 2};
            crt_Ks{2} = results.outputs.K{k, i, 2};
        end

        % calculate the per-environment correlations
        for j = 1:n_gamma
            crt_K = crt_Ks{j};
            for p = 1:length(abs_proxies)/2
                crt_proxy_name = abs_proxies{2*p-1};
                crt_proxy_fct = abs_proxies{2*p};
                crt_proxy = crt_proxy_fct(crt_Qs{j});
                for q = 1:length(corr_fcts)/2
                    crt_corr_name = corr_fcts{2*q-1};
                    crt_corr_fct = corr_fcts{2*q};
                    crt_corr = crt_corr_fct(crt_proxy, crt_K);
                    results.corrs.(crt_proxy_name).(crt_corr_name)(k, i, j) = crt_corr;
                end
            end
        end
        
        if ~isempty(Gamma2_fct)
            % calculate the between-environments correlations
            for p = 1:length(diff_proxies)/2
                crt_proxy_name = diff_proxies{2*p-1};
                crt_proxy_fct = diff_proxies{2*p};
                crt_proxy = crt_proxy_fct(crt_Qs{1}, crt_Qs{2});
                crt_diffK = crt_Ks{2} - crt_Ks{1};
                for q = 1:length(corr_fcts)/2
                    crt_corr_name = corr_fcts{2*q-1};
                    crt_corr_fct = corr_fcts{2*q};
                    crt_corr = crt_corr_fct(crt_proxy, crt_diffK);
                    results.corrs.(crt_proxy_name).(crt_corr_name)(k, i) = crt_corr;
                end
            end
        end
    end
end

% summarize correlations (calculate means, standard deviations, etc.)
results.summaries = struct;
summary_fcts = {'mean', @(c) mean(c, 1), 'median', @(c) median(c, 1), ...
    'std', @(c) std(c, [], 1), 'low', @(c) quantile(c, 0.2, 1), ...
    'high', @(c) quantile(c, 0.8, 1)};
proxy_fields = fieldnames(results.corrs);
for i = 1:length(proxy_fields)
    crt_corrs_struct = results.corrs.(proxy_fields{i});
    corr_fields = fieldnames(crt_corrs_struct);
    for j = 1:length(corr_fields)
        crt_corrs = crt_corrs_struct.(corr_fields{j});
        for k = 1:length(summary_fcts)/2
            summary_fct_name = summary_fcts{2*k-1};
            summary_fct = summary_fcts{2*k};
            summary = summary_fct(crt_corrs);
            
            summary_name = [corr_fields{j} '_' summary_fct_name];
            results.summaries.(proxy_fields{i}).(summary_name) = summary;
        end
    end
end

end
