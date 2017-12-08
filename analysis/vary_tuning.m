%% Generate sample sensing matrices

S_low = generate_random_sensing(24, 50, 0.05, 200, 'shuffle', false);
S_mid = generate_random_sensing(24, 50, 0.20, 200, 'shuffle', false);
S_high = generate_random_sensing(24, 50, 1.0, 200, 'shuffle', false);

S_low_shuffled = S_low(:, randperm(50));
S_mid_shuffled = S_mid(:, randperm(50));
S_high_shuffled = S_high(:, randperm(50));

%% Image and save them

all_mats = containers.Map;

all_mats('S_low_tuning') = S_low;
all_mats('S_mid_tuning') = S_mid;
all_mats('S_high_tuning') = S_high;

all_mats('S_low_tuning_shuffled') = S_low_shuffled;
all_mats('S_mid_tuning_shuffled') = S_mid_shuffled;
all_mats('S_high_tuning_shuffled') = S_high_shuffled;

map_keys = all_mats.keys;
for i = 1:length(map_keys)
    crt_name = map_keys{i};
    crt_S = all_mats(crt_name);

    fig = figure;
    fig.Units = 'inches';
    
    fig_x = 2.5;
    fig_y = 1.4;
    
    fig.Position = [fig.Position(1:2) fig_x fig_y];
    
    ax = axes;
    ax.Units = 'inches';
    
    edge_x = 0.25;
    edge_y = 0.25;
    
    margin_x = 0.1;
    
    ax_x = (fig_x - edge_x - margin_x);
    ax_y = ax_x*24/50;
    
    ax.Position = [edge_x edge_y ax_x ax_y];
    
    imagesc(crt_S, [0 1]);
    
    beautifygraph;
    axis equal;
    box on;
    set(gca, 'xminortick', 'off', 'yminortick', 'off', 'xtick', [], 'ytick', []);
    xlabel('Odorants');
    ylabel('Receptors');
    
    preparegraph;
    
    safe_print(fullfile('figs', 'tuning_figs', [crt_name '.pdf']));
end