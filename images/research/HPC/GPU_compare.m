clc; clear all; close all;

% GPU model names
models = {'K20', 'K40', 'P100', 'V100', 'A100', 'H100'};
years  = [ 2012   2013   2016    2017    2020    2022];

% Metric data
metrics = {
    % Name            Data                    Format    YLabel               Title                          Filename
    'CUDA Cores',     [2496 2880 3584 5120 6912 16896],  '%d',  'CUDA Cores',       'Compute Core Comparison',       'GPU_Cores.png';  % 网页6][6,7,9,12](@ref)[^14
    'FP64 TFLOPS',    [1.17 1.4 4.7 7.0 9.7 34.0],      '%.1f', 'FP64 (TFLOPS)',    'FP64 Performance Comparison',   'GPU_FP64.png';    % 网页6][7,9,12,14](@ref)[^16
    'Memory',         [5 12 16 32 40 80],                '%d',  'Memory (GB)',      'Memory Capacity Comparison',   'GPU_Memory.png'; % 网页6][8,9,12](@ref)[^14
    'Bandwidth',      [208 288 732 900 1555 3350],       '%d',  'Bandwidth (GB/s)', 'Memory Bandwidth Comparison',   'GPU_Bandwidth.png'; % 网页8][9,12,14](@ref)[^16
    };

% Color scheme (one per GPU model)
colors = [
    0.55 0.10 0.20    % K20
    0.15 0.45 0.70    % K40
    0.30 0.65 0.75    % P100
    0.55 0.77 0.85    % V100
    0.95 0.55 0.25    % A100
    0.40 0.70 0.40    % H100
    ];

% Generate plots for each metric
for m = 1:size(metrics, 1)
    % Create figure with consistent sizing
    figure('Color', 'w', 'Position', [200+(m/2-1)*300, 200+(m/2-1)*225, 600, 450]); hold on;
    hold on;

    % Create bar chart
    data = metrics{m, 2};
    hBar = bar(data);

    % Apply styling
    hBar.FaceColor = 'flat';
    hBar.CData = colors;

    % Configure axes
    ax = gca;
    ax.XTick = 1:6;
    ax.XTickLabel = models;
    ax.FontSize = 18;
    ax.FontName = 'Arial';
    ax.Position = [0.20 0.20 0.75 0.75];
    
    yMax = max(data);

    % Set dynamic Y-axis limits
    if(m == 1) 
        ylim([0 20000])
    end
    if(m == 2)
        ylim([0 40])
    end
    if(m == 3)
        ylim([0 100])
    end
    if(m == 4)
        ylim([0 4000])
    end
    xlim([0.3 6.7])

    % Add labels
    xlabel('GPU Architecture', 'FontSize', 24);
    ylabel(metrics{m, 4}, 'FontSize', 24);

    % Add data labels
    labelOffset = 0.01 * yMax; % Dynamic positioning
    for i = 1:6
        text(i, data(i) + labelOffset,...
            num2str(data(i), metrics{m, 3}),...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'bottom',...
            'FontSize', 14,...
            'FontName', 'Arial');
    end

    for i = 1:size(data, 2)
        y_pos = cumsum(data(:,i));
        for j = 1:size(data, 1)
            text(i, y_pos(j) - data(j,i)/2, ...
                sprintf('%.0f', years(i)), ...
                'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', ...
                'FontSize', 14, ...
                'Color', 'w', ...
                'FontWeight', 'bold');
        end
    end

    xlims = xlim;
    ylims = ylim;
    h_border = plot([xlims(1), xlims(2), xlims(2), xlims(1), xlims(1)], ...
        [ylims(1), ylims(1), ylims(2), ylims(2), ylims(1)], ...
        'k-', 'LineWidth', 2);
    set(get(get(h_border, 'Annotation'), 'LegendInformation'), ...
        'IconDisplayStyle', 'off');
    set(gca, 'Position', [0.18, 0.18, 0.75, 0.75]);  % Margins adjustment
    set(gca, 'LineWidth', 2, 'TickDir', 'in', 'TickLength', [0.02 0.02])

    % Save image
    print('-dpng', '-r600', metrics{m, 6});
end