function setup_paths
% SETUP_PATHS Temporarily add paths to all the olfactory receptor
% distribution code to Matlab's search path.
%   SETUP_PATHS temporarily adds all the folders related to the olfactory
%   receptor distribution code to Matlab's search path. This gets reset
%   when Matlab quits.

PWD = pwd;
folders = {'analysis', 'data', 'dynamics', 'graphing', 'optimization', ...
    'sandbox', 'utils'};

for i = 1:length(folders)
    addpath(fullfile(PWD, folders{i}));
end

end
