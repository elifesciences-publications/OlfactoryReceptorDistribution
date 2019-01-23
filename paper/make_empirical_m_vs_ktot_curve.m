% Look at animal data of number of receptor types vs. size (of epithelium,
% of animal, of brain, etc.)

%% Generate a data table

animal_data = table;

% choose animals
animal_data.Name = {'sheep' ; 'guinea pig' ; 'cat' ; 'dog' ; 'monkey' ; 'human' ; ...
    'rat' ; 'mouse' ; 'marmoset'};

% set up fields for olfactory surface area, animal mass, and number of OR genes
animal_data.OlfactorySurface = nan(height(animal_data), 1);
animal_data.Mass = nan(height(animal_data), 1);
animal_data.M = nan(height(animal_data), 1);

% set units
animal_data.Properties.VariableUnits = {'', 'cm^2', 'kg', ''};

% set the olfactory surface areas (in cm^2)
% sheep, guinea pig, monkey from Moulton et al. (1967)
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'sheep')) = 8.84;
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'guinea pig')) = 1.12;
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'monkey')) = 4.16;
% cat, dog (German shepherd), human from Pihlstroem et al. (2005)
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'cat')) = 27.91;
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'dog')) = 139.00;
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'human')) = 11.25;
% rat, mouse from Gross et al. (1981)
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'rat')) = 6.75;
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'mouse')) = 1.37;
% marmoset from Smith et al. (2014)
animal_data.OlfactorySurface(strcmp(animal_data.Name, 'marmoset')) = 0.4869;

% set the number of intact OR genes from Niimura et al. (2014)
animal_data.M(strcmp(animal_data.Name, 'guinea pig')) = 796;
animal_data.M(strcmp(animal_data.Name, 'dog')) = 811;
animal_data.M(strcmp(animal_data.Name, 'human')) = 396;
animal_data.M(strcmp(animal_data.Name, 'rat')) = 1207;
animal_data.M(strcmp(animal_data.Name, 'mouse')) = 1130;
animal_data.M(strcmp(animal_data.Name, 'marmoset')) = 366;
% more numbers of OR genes from Hughes et al. (2018)
animal_data.M(strcmp(animal_data.Name, 'cat')) = 856;
animal_data.M(strcmp(animal_data.Name, 'sheep')) = 962;

% set the masses (in kg)
% guinea pig, human, sheep, and cat from Rousseeuw, P.J. & Leroy, A.M. (1987) Robust Regression and Outlier Detection. Wiley, p. 57.
animal_data.Mass(strcmp(animal_data.Name, 'guinea pig')) = 1.04;
animal_data.Mass(strcmp(animal_data.Name, 'human')) = 62;
animal_data.Mass(strcmp(animal_data.Name, 'sheep')) = 55.5;
animal_data.Mass(strcmp(animal_data.Name, 'cat')) = 3.3;
% German shepherd from FEDERATION CYNOLOGIQUE INTERNATIONALE (AISBL)
% (http://www.fci.be/Nomenclature/Standards/166g01-en.pdf)
animal_data.Mass(strcmp(animal_data.Name, 'dog')) = 30.0; % German shepherd
% rat, mouse from Gross et al. (1981)
animal_data.Mass(strcmp(animal_data.Name, 'rat')) = 0.288;
animal_data.Mass(strcmp(animal_data.Name, 'mouse')) = 0.033;
% marmoset (Callithrix jacchus) from Smith et al. (2014)
animal_data.Mass(strcmp(animal_data.Name, 'marmoset')) = 0.250;

%% Ratio of olfactory surface to animal surface works well

fig = figure;
fig.Units = 'inches';
fig.Position(3:4) = [2.4 1.4];
fig.Color = [1 1 1];

opts = {'fitopts', {'legend', false, 'style', {'--', 'color', [0.7 0.7 0.7], 'linewidth', 1}}, ...
    'scatteropts', {'color', hex2color('BC2E2C'), 'alpha', 1}};

mask = isfinite(animal_data.OlfactorySurface);
% x_values = animal_data.OlfactorySurface(mask) ./ animal_data.Mass(mask).^(2/3);
% x_values = animal_data.OlfactorySurface(mask) ./ animal_data.Mass(mask);
% x_values = log10(animal_data.OlfactorySurface(mask) ./ animal_data.Mass(mask).^0.3);
x_values = log10(animal_data.OlfactorySurface(mask) ./ animal_data.Mass(mask).^0.2);
y_values = animal_data.M(mask);
names = animal_data.Name(mask);

scatterfit(x_values, y_values, opts{:});
xlabel('area_{epithelium} / mass^{0.2}');
ylabel('# intact OR genes');

ylim([250, 1300]);
xlim([-0.3 2.0]);

for i = 1:length(names)
    shift_x = 0.06;
    shift_y = 10;
    if ~strcmp(names{i}, 'dog')
        switch names{i}
            case 'marmoset'
                shift_x = 0.04;
                shift_y = -60;
            case 'dog'
                shift_x = 0;
                shift_y = -90;
        end
        text(x_values(i) + shift_x, y_values(i) + shift_y, names{i}, 'fontsize', 8);
    else
        text(x_values(i), y_values(i) - 30, 'German', ...
            'fontsize', 8, 'horizontalalignment', 'center', ...
            'verticalalignment', 'top');
        text(x_values(i), y_values(i) - 150, 'shepherd', ...
            'fontsize', 8, 'horizontalalignment', 'center', ...
            'verticalalignment', 'top');
    end
end

set(gca, 'xtick', [0 1 2], 'xticklabel', {'1', '10', '100'});

beautifygraph('linewidth', 0.5, 'minorticks', 'off', 'fontscale', 0.667, 'ticksize', 12);
preparegraph;

safeprint(fullfile('figs', 'empirical_m_vs_ktot'));
