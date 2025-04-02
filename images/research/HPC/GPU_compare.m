clc; clear all; close all;

% GPU model names
models = {'K20', 'K40', 'P100', 'V100', 'A100', 'H100'};

% Metric data (hypothetical values - verify with actual specs)
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
    figure('Color', [1 1 1],...
           'Position', [100 100 360 280],...
           'PaperUnits', 'centimeters',...
           'PaperPosition', [0 0 18 14],...
           'PaperSize', [18 14],...
           'PaperPositionMode', 'manual');
    
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
    ax.YGrid = 'on';
    ax.GridLineStyle = '--';
    ax.FontSize = 8;
    ax.FontName = 'Arial';
    ax.Position = [0.15 0.18 0.75 0.75];
    
    % Set dynamic Y-axis limits
    yMax = 1.1 * max(data);
    ylim([0 yMax]);
    
    % Add labels
    xlabel('GPU Architecture', 'FontSize', 9, 'FontWeight', 'bold');
    ylabel(metrics{m, 4}, 'FontSize', 9, 'FontWeight', 'bold');
    title(metrics{m, 5}, 'FontSize', 10, 'FontWeight', 'bold');
    
    % Add data labels
    labelOffset = 0.05 * yMax; % Dynamic positioning
    for i = 1:6
        text(i, data(i) + labelOffset,...
            num2str(data(i), metrics{m, 3}),...
            'HorizontalAlignment', 'center',...
            'VerticalAlignment', 'bottom',...
            'FontSize', 8,...
            'FontName', 'Arial');
    end
    
    % Save image
    print('-dpng', '-r600', metrics{m, 6});
end