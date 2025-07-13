clc; clear; close all;

% 数据设置：设备名称, 延迟(ns), 带宽(GB/s)
devices = {'Intel Xeon E5-2680v3', 'AMD EPYC 7763', ...
           'NVIDIA V100', 'NVIDIA A100', 'NVIDIA H100'};

latency_ns = [80, 90, 400, 500, 550];
bandwidth_gbps = [70, 204.8, 900, 1555, 3350];
% data for sure: 900, 155, 3350

types = {'CPU', 'CPU', 'GPU', 'GPU', 'GPU'};

% 配色与形状（CPU 蓝色方块，GPU 橙色圆点）
colors = {[0.2, 0.4, 0.7], [0.2, 0.4, 0.7], ...
          [0.9, 0.4, 0.1], [0.9, 0.4, 0.1], [0.9, 0.4, 0.1]};
markers = {'s', 's', 'o', 'o', 'o'};

% 绘图初始化
figure('Color','w', 'Position',[100,100,600,450]); hold on;

% 绘制每个设备的数据点
for i = 1:length(devices)
    loglog(latency_ns(i), bandwidth_gbps(i), markers{i}, ...
        'MarkerSize', 12, ...
        'MarkerFaceColor', colors{i}, ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 1.5);
    
    % 添加文本标注（带白色背景框，避免重叠）
    text(latency_ns(i)*1.2, bandwidth_gbps(i), devices{i}, ...
        'FontSize', 18, ...
        'Margin', 2);
end

xticks = [100 1000];
xlabels = {'$10^{2}$', '$10^{3}$'};

yticks = [100 1000];
ylabels = {'$10^{2}$', '$10^{3}$'};

% 设置坐标轴和标题
set(gca, 'XScale', 'log', 'YScale', 'log', ...
    'XTick', xticks, 'XTickLabel', xlabels, ...
    'YTick', yticks, 'YTickLabel', ylabels, ...
    'FontSize', 28,  ...
    'TickLength', [0.025 0.025], ...
    'LineWidth', 2, ...
    'TickLabelInterpreter', 'latex');
xlabel('访问延迟 (ns)', 'FontSize', 28);
ylabel('带宽 (GB/s)', 'FontSize', 28);

% 坐标范围和网格
xlim([10 5000]);
ylim([10 5000]);

xlims = xlim;
ylims = ylim;
h_border = plot([xlims(1), xlims(2), xlims(2), xlims(1), xlims(1)], ...
                [ylims(1), ylims(1), ylims(2), ylims(2), ylims(1)], ...
                'k-', 'LineWidth', 2);
set(get(get(h_border, 'Annotation'), 'LegendInformation'), ...
    'IconDisplayStyle', 'off');
set(gca, 'Position', [0.20, 0.25, 0.75, 0.70]);  % Margins adjustment

% 保存为出版质量图像
print('compareLatencyBandwidth.png', '-dpng', '-r600');
