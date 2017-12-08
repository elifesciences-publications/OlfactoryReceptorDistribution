% check the effects of a change in Ktot

%% Load fly sensing matrix

% fly sensing matrix
datafile = open('flyResponsesWithNames.mat');
data = datafile.relRates;

S_fly = data';
sigma_fly = 1000;

%% Make environment

% random environment, but make it reproducible
rng(2334);

Gamma = generate_environment('rnd_diag_const', size(S_fly, 2));

%% Get distribution at various values for Ktot

Ktot_values = [0 logspace(log10(0.1), log10(15000), 100)];
[K, info_values, Q, info_fct] = calculate_optimal_dist(S_fly/sigma_fly, Gamma, Ktot_values);

[M, N] = size(S_fly);

% find where each receptor kicks in
minK_by_type = ones(M, 1);
K_nz = bsxfun(@rdivide, K, Ktot_values) > 0.004;
for i = 1:M
    idx = find(~K_nz(i, :), 1, 'last');
    if ~isempty(idx)
        minK_by_type(i) = idx+1;
    end
end
[minK_ordered, rec_order] = sort(minK_by_type);

%% Plot "natural" receptor usage

fig = figure;
fig.Units = 'inches';
fig.Position = [fig.Position(1:2) 5.8 1.9];

% plot the receptor distribution
% subplot(1, 3, 2);
ax1 = axes;
ax1.Units = 'normalized';
ax1.OuterPosition = [1/3 0 1/3 1];
%idx_max = 62;
idx_max = 32;
% scale noise to bring K values to reasonable range
K_scaling = 50;
h_area = area(K_scaling*Ktot_values(1:idx_max), K_scaling*K(:, 1:idx_max)', 'linewidth', 1);
hold on;
for i = 1:M-1
    if minK_ordered(i+1) <= length(Ktot_values) && minK_ordered(i+1) <= idx_max
        line(repmat(K_scaling*Ktot_values(minK_ordered(i+1)), [2 1]), ylim, ...
            'linestyle', ':', 'color', [0.7 0.7 0.7], 'linewidth', 1);
    end
end
% axis equal;
xlim([0 K_scaling*Ktot_values(idx_max)]);
ylim([0 K_scaling*Ktot_values(idx_max)]);

xlh = xlabel('OSN number');
xlh.Units = 'characters';
xlh.Position(2) = xlh.Position(2) - 0.2;
ylabel('OSNs by receptor type');
% find out which receptors we've displayed, and in what order they appeared
rec_dispd = rec_order(minK_ordered <= idx_max);
h_leg = legend(h_area(rec_dispd), datafile.orNames(rec_dispd), 'location', 'northwest', ...
    'fontsize', 6);
h_leg.Position = [h_leg.Position(1) - 0.025 h_leg.Position(2) + 0.08 h_leg.Position(3:4)];
legend('boxoff');
beautifygraph('fontscale', 0.833);

ax = gca;
ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;

% plot the information per OSN
% subplot(1, 3, 1);
ax2 = axes;
ax2.Units = 'normalized';
ax2.OuterPosition = [0 0 1/3 1];
h_area = area(K_scaling*Ktot_values, K_scaling*K', 'linewidth', 0.5);
hold on;
% for i = 1:M-1
%     if minK_ordered(i+1) <= length(Ktot_values)
%         line(repmat(K_scaling*Ktot_values(minK_ordered(i+1)), [2 1]), ylim, ...
%             'linestyle', ':', 'color', [0.7 0.7 0.7], 'linewidth', 1);
%     end
% end
% axis equal;
xlim([0 K_scaling*Ktot_values(end)]);
ylim([0 K_scaling*Ktot_values(end)]);

xlabel('OSN number');
ylabel('OSNs by receptor type');
% find out which receptors we've displayed, and in what order they appeared
% rec_dispd = rec_order(minK_ordered <= idx_max);
% legend(h_area(rec_dispd), datafile.orNames(rec_dispd), 'location', 'northwest', ...
%     'fontsize', 8);
% legend('boxoff');
beautifygraph('fontscale', 0.833);

ax = gca;
min_nz_value = min(ax.XTick(ax.XTick > 0));
scale_exp = floor(log10(min_nz_value));
iif = @(varargin) varargin{2 * find([varargin{1:2:end}], 1, 'first')}();
ticks_to_str = @(x) iif(abs(x) > eps, ...
    [num2str(x/10^scale_exp) '\times10^{' int2str(scale_exp) '}'], ...
    true, '0');

ax = gca;
ax.XTick = [0 4*10^5];
ax.XTickLabel = arrayfun(ticks_to_str, ax.XTick, 'uniform', false);
%ax.YTickLabel = arrayfun(ticks_to_str, ax_high_snr.YTick, 'uniform', false);

ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;

% subplot(1, 3, 3);
ax3 = axes;
ax3.Units = 'normalized';
ax3.OuterPosition = [2/3 0 1/3 1];
rectype_count = zeros(size(Ktot_values));
for i = 1:length(minK_ordered)
    rectype_count(minK_ordered(i):end) = rectype_count(minK_ordered(i):end) + 1;
end
semilogx(K_scaling*Ktot_values, rectype_count, 'linewidth', 1);
xlabel('OSN number');
ylabel('Active receptor types');
xlim([K_scaling*Ktot_values(2), K_scaling*Ktot_values(end)]);
ylim([0 M+1]);
beautifygraph('fontscale', 0.833);

ax = gca;
ax.XTick = [10^1 10^2 10^3 10^4 10^5];

ax.XMinorTick = 'off';
ax.YMinorTick = 'off';
ax.TickLength = 3*ax.TickLength;

preparegraph('edge', 0);

ymax = max([ax1.Position(2) ax2.Position(2) ax3.Position(2)]);
hmin = min([ax1.Position(4) ax2.Position(4) ax3.Position(4)]);
ax1.Position(2) = ymax;
ax2.Position(2) = ymax;
ax3.Position(2) = ymax;
ax1.Position(4) = hmin;
ax2.Position(4) = hmin;
ax3.Position(4) = hmin;
% disp(ax1.Position);
% disp(ax2.Position);
% disp(ax3.Position);

safe_print(fullfile('figs', 'natural_receptor_dist.pdf'));