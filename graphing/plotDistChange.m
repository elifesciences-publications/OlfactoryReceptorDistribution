function plotDistChange(K1, K2, varargin)
% plotDistChange Make a plot showing a change in receptor distribution.
%   plotDistChange(K1, K2) makes a plot showing how the receptor
%   distribution changes from K1 to K2.
%
%   Options:
%    'beautifyopts':
%       Cell array of options to be passed to beautifygraph.
%    'method':
%       How to show the change. This can be
%        'plot': a plot where thin lines connect the distribution K1 drawn
%                in dots on the left, to the distribution K2 drawn in dots
%                on the right
%        'bar':  a barplot showing the change K2 - K1 with different colors
%                for increases vs. decreases in abundance
%    'receptornames':
%       Cell array of receptor names. For now this is only used with the
%       'bar' option.

% parse optional arguments
parser = inputParser;
parser.CaseSensitive = true;
parser.FunctionName = mfilename;

parser.addParameter('beautifyopts', {}, @(c) isempty(c) || (iscell(c) && isvector(c)));
parser.addParameter('method', 'plot', @(s) ischar(s) && isvector(s) && ...
    ismember(s, {'plot', 'bar'}));
parser.addParameter('receptornames', {}, @(c) isempty(c) || (iscell(c) && isvector(c)));

% handle the 'defaults' option
if nargin == 1 && strcmp(K1, 'defaults')
    parser.parse;
    disp(parser.Results);
    return;
end

% parse special arguments
parser.parse(varargin{:});
params = parser.Results;

if ~isvector(K1) || ~isvector(K2)
    error([mfilename ':badinp'], 'The inputs are not vectors.');
end

col1 = [0.21 0.17 0.53];
col2 = [0.98 0.80 0.17];

M = length(K1);
if length(K2) ~= M
    error([mfilename ':badsz'], 'The lengths of the two distributions do not match.');
end

% normalize
K1 = K1(:)' / sum(K1);
K2 = K2(:)' / sum(K2);

switch params.method
    case 'plot'
        pos1 = zeros(1, M);
        pos2 = ones(1, M);
        plot([pos1 ; pos2], [K1 ; K2], 'k-');
        hold on;
        plot(pos1, K1, 'd', 'markerfacecolor', col1, 'markeredgecolor', 'none');
        plot(pos2, K2, 'd', 'markerfacecolor', col2, 'markeredgecolor', 'none');
        
        ylabel('Fraction of OSN');
        
        ylim([0 2/M]);
        
        beautifygraph(params.beautifyopts{:});
        
        ax = gca;
        ax.Box = 'off';
        ax.XMinorTick = 'off';
        ax.YMinorTick = 'off';
        ax.TickLength = [0 0];
        ax.XTick = [0 1];
        ax.XTickLabel = {'env1', 'env2'};
    case 'bar'
        diffK = K1 - K2;
        bar(diffK.*(diffK > 0), 0.8, 'facecolor', col1, 'edgecolor', 'none');
        hold on;
        bar(diffK.*(diffK < 0), 0.8, 'facecolor', col2, 'edgecolor', 'none');
        
        crt_yl = max(abs(diffK));
        xlim([0 M+1]);
        ylim([-crt_yl crt_yl]);
        
        beautifygraph(params.beautifyopts{:});
        
        ax = gca;
        ax.YTick = [];
        ax.XTick = 1:M;
        if ~isempty(params.receptornames)
            ax.XTickLabel = params.receptornames;
        end
        ax.XMinorTick = 'off';
        ax.FontSize = 6;
        
        xtickangle(60);        
end

end