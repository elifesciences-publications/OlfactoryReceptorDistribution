% Make plots meant to show generic properties of the optimization.

%% Load Hallem&Carlson sensing data

sensing_fly = open('data/flyResponsesWithNames.mat');
S_fly = sensing_fly.relRates';

%% Plot a sample of the Hallem&Carlson sensing data

odor_choices = {'1-hexanol', 'E2-hexenol', 'E3-hexenol', 'Z3-hexenol', ...
    '2-heptanone', 'butyl acetate', 'pentyl acetate', 'methyl hexanoate', ...
    'ethyl hexanoate', 'methyl benzoate', 'ethyl benzoate'};
or_choices = {'Or2a', 'Or19a', 'Or35a', 'Or47b', 'Or67a', 'Or85b'};

odor_idxs = cellfun(@(s) find(strcmp(sensing_fly.odorNames, s)), odor_choices);
or_idxs = cellfun(@(s) find(strcmp(sensing_fly.orNames, s)), or_choices);

fig = figure;
ax = axes;
ax.Position = [0 0 1 1];

% cmap = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64);
cmap = divergent([0.177 0.459 0.733], [0.737 0.180 0.172], 64, [0.969 0.957 0.961]);

cl = quantile(flatten(S_fly(or_idxs, odor_idxs)), 0.9);

plotMatrix(S_fly(or_idxs, odor_idxs), [-cl cl]);
colormap(cmap);

axis off;
