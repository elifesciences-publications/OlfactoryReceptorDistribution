% Make plots showing dependence of number of expressed receptor types on
% total number of neurons.

%% Load the results

load(fullfile('save', 'Ktot_dependence.mat'));

%% Make the plots

% mouse --> mammal
% artificial --> tuning?
% skip artificial narrow & wide

to_plot = {'fly', 'fly_scrambled', 'mouse', 'mouse_scrambled', ...
    'artificial', 'gaussian', 'binary', 'signed'};
renames = containers.Map;
renames('mouse') = 'mammal';
renames('mouse_scrambled') = 'mammal scrambled';
renames('artificial') = 'tuning';

fig = figure;
fig.Units = 'inches';
fig.Color = [1 1 1];
fig.Position(3:4) = [5.5 2.3];

nx = 4;
ny = 2;
n_plots = 8;

for i = 1:n_plots
    crt_S_choice = to_plot{i};
    
    crt_x = mod(i-1, nx);
    crt_y = floor((i-1)/nx);
    
    ax = axes;
    ax.OuterPosition = [crt_x/nx 1-(crt_y+1)/ny 1/nx 1/ny];

    % make plot showing how the number of receptor expressed in non-zero
    % numbers of neurons grows with the total number of neurons
    crt_results = results(crt_S_choice);
    crt_Ktot_values = crt_results.Ktot_values;
    crt_minK_ordered = crt_results.steps_Ktot_idx;
    crt_S = crt_results.S;
    crt_K = crt_results.K;
    
%     crt_K_nz = (bsxfun(@rdivide, crt_K, crt_Ktot_values) > 0.004);
    crt_K_nz = (bsxfun(@rdivide, crt_K, crt_Ktot_values) > 0.01);
    
    rectype_count = sum(crt_K_nz, 1);
    semilogx(crt_Ktot_values, rectype_count, 'linewidth', 1);
    xlabel('Number of OSNs');
    ylabel('# receptor types');
    xlim([crt_Ktot_values(2), crt_Ktot_values(end)]);
    ylim([0 size(crt_S, 1)+1]);
    
    if renames.isKey(crt_S_choice)
        crt_name = renames(crt_S_choice);
    else
        crt_name = crt_S_choice;
        crt_name(crt_name == '_') = ' ';
    end
    title(crt_name);
    
    beautifygraph(ax, 'fontscale', 0.667, 'ticksize', 10, 'linewidth', 0.5, ...
        'titlesize', 11, 'labelsize', 10);

    % select a good crop of ticks
    ax.XTick = [10^1 10^2 10^3 10^4 10^5];

    ax.XMinorTick = 'off';
    ax.YMinorTick = 'off';
    ax.TickLength = 3*ax.TickLength;
end
preparegraph('edge', 0);

safeprint(fullfile('figs', 'ktot_dependence_si'));
