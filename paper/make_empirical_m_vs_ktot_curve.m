% Look at animal data of number of receptor types vs. size (of epithelium,
% of animal, of brain, etc.)

%% Data from Niimura et al. (2014)

animal_data = table;
animal_data.Order = {'Proboscidea' ; 'Cetartiodactyla' ; 'Carnivora' ; ...
    'Perissodactyla' ; 'Lagomorpha' ; 'Rodentia' ; 'Rodentia' ; 'Rodentia' ; ...
    'Primates' ; 'Primates' ; 'Primates' ; 'Primates' ; 'Primates'};
animal_data.Name = {'african elephant' ; 'cow' ; 'dog' ; 'horse' ; 'rabbit' ; ...
    'guinea pig' ; 'rat' ; 'mouse' ; 'marmoset' ; 'macaque' ; 'orangutan' ; ...
    'chimpanzee' ; 'human'};
animal_data.M = [1948 ; 1186 ; 811 ; 1066 ; 768 ; 796 ; 1207 ; 1130 ; ...
    366 ; 309 ; 296 ; 380 ; 396];
animal_data.M_trunc = [89 ; 41 ; 11 ; 23 ; 22 ; 26 ; 52 ; 0 ; 27 ; ...
    17 ; 37 ; 19 ; 0];
animal_data.M_pseudo = [2230 ; 1057 ; 278 ; 1569 ; 256 ; 1340 ; 508 ; ...
    236 ; 231 ; 280 ; 488 ; 414 ; 425];

% XXX weights from Google
animal_data.MassLo = [2700 ; 450 ;  2 ;  380 ; 2 ; 0.7 ; 0.3 ; 0.005 ; 0.10 ; 5 ;  30 ; 25 ;  50];
animal_data.MassHi = [5400 ; 820 ; 90 ; 1000 ; 4 ; 1.2 ; 0.5 ; 0.020 ; 0.25 ; 8 ; 170 ; 60 ; 100];
animal_data.MassMid = exp(0.5*(log(animal_data.MassLo) + log(animal_data.MassHi)));

% XXX number of neurons (whole body or full CNS) from Wikipedia
% XXX using African elephant number
% XXX cow estimate from a different website
% XXX horse full brain estimate based on number for cerebral cortex
% XXX cortical estimate for cow and rabbit are from full brain estimates
animal_data.Neurons =           [257e9 ;   3e9 ; 2.253e9 ; 6.4e9 ; 494e6 ; 240e6 ; 200e6 ; 71e6 ; 636e6 ; 6.38e9 ; 32.6e9 ;  28e9 ; 86e9];
animal_data.CorticalNeurons =   [5.6e9 ; 0.6e9 ;   530e6 ; 1.2e9 ; 114e6 ;  44e6 ;  31e6 ; 14e6 ; 245e6 ; 1.71e9 ;  8.9e9 ; 6.2e9 ; 16e9];

% from Moulton et al. (1967)
animal_data.OlfactorySurface = [NaN ; NaN ; 9.76 ; NaN ; NaN ; 1.12 ; 1.12 ; NaN ; NaN ; NaN ; NaN ; NaN ; 3.08];

animal_data.Properties.VariableUnits = {'', '', '', '', '', 'kg', 'kg', 'kg', '', '', 'cm^2'};

%% Ratio of olfactory surface to animal surface works well

fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [2.4 1.4];
fig.Color = [1 1 1];

opts = {'fitopts', {'legend', false, 'style', {'--', 'color', [0.7 0.7 0.7], 'linewidth', 1}}, ...
    'scatteropts', {'color', hex2color('BC2E2C'), 'alpha', 1}};

mask = isfinite(animal_data.OlfactorySurface);
x_values = animal_data.OlfactorySurface(mask) ./ animal_data.MassMid(mask).^(2/3);
y_values = animal_data.M(mask);
names = animal_data.Name(mask);

scatterfit(x_values, y_values, opts{:});
xlabel('area_{epithelium} / mass^{2/3}');
ylabel('# intact OR genes');

ylim([300, 1300]);
xlim([0 2.5]);

for i = 1:length(names)
    shift_x = 0.08;
    shift_y = -20;
    if strcmp(names(i), 'guinea pig')
        shift_x = -0.85;
        shift_y = 60;
    end
    text(x_values(i) + shift_x, y_values(i) + shift_y, names{i}, 'fontsize', 8);
end

beautifygraph('linewidth', 0.5, 'minorticks', 'off', 'fontscale', 0.667, 'ticksize', 12);
preparegraph;

safeprint(fullfile('figs', 'empirical_m_vs_ktot'));

%% Make some more plots

fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [6 4];
fig.Color = [1 1 1];

opts = {'fitopts', {'legend', 'cp', 'style', {'--', 'color', [0.7 0.7 0.7], 'linewidth', 1}}, ...
    'scatteropts', {'color', hex2color('BC2E2C'), 'alpha', 1}};

subplot(2, 2, 1);
scatterfit(log10(animal_data.MassMid), animal_data.M, opts{:});
xlabel('log_{10} mass/kg')
ylabel('# intact OR genes');

subplot(2, 2, 2);
scatterfit(log10(animal_data.OlfactorySurface), animal_data.M, opts{:});
xlabel('log_{10} area_{epithelium}/cm^2')
ylabel('# intact OR genes');

subplot(2, 2, 3);
scatterfit(log10(animal_data.Neurons), animal_data.M, opts{:});
xlabel('log_{10} n_{neurons}')
ylabel('# intact OR genes');

subplot(2, 2, 4);
scatterfit(log10(animal_data.CorticalNeurons), animal_data.M, opts{:});
xlabel('log_{10} n_{cortical neurons}')
ylabel('# intact OR genes');
