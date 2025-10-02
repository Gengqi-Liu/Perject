% MATLAB Code: Visualization of Hybrid Response Across Angles 10, 20, 30, 40, 50 (仅时域)
% -------------------------------------------------------------------------
% 目的: 
% 1. 加载所有角度的 FIR/IIR 混合响应结果。
% 2. 绘制每个通道的指定五个角度 (10, 20, 30, 40, 50) 的时域对比图。
% -------------------------------------------------------------------------

% 0. 初始化与数据加载
% -------------------------------------------------------------------------
data_filename = 'Cross_Angle_Hybrid_N7500_M7500_Locked_A53_Final.mat';
try
    load(data_filename);
    
    % 检查关键结构体
    if ~exist('All_Hybrid_Results', 'var') || ~exist('fs', 'var') || ~exist('t_cut_ms', 'var')
        error('加载的文件中缺少 All_Hybrid_Results, fs, 或 t_cut_ms 等关键变量。');
    end
    
    disp(['数据文件 ', data_filename, ' 加载成功。']);
    
    % 提取共享参数
    N_fir_samples = length(All_Hybrid_Results.A1.Ch11.h_hybrid);
    channel_names = fieldnames(All_Hybrid_Results.A1);
    
    t_samples = (0:N_fir_samples-1) / fs;
    t_cut_s = t_cut_ms / 1000; % 80ms 拼接点
    
catch ME
    disp(['错误: 无法加载文件 ', data_filename, '。请检查文件路径和内容。']);
    return;
end

% 1. 定义要对比的角度和绘图样式
% -------------------------------------------------------------------------
angles_to_plot = [10, 20, 30, 40, 50];
num_angles_plot = length(angles_to_plot);

% 为五个角度定义颜色和线型
line_specs = {
    'b-',   % 角度 10: 蓝色实线
    'g-',   % 角度 20: 绿色实线
    'r-',   % 角度 30: 红色实线
    'c--',  % 角度 40: 青色虚线
    'm--'   % 角度 50: 洋红色虚线
};

disp(['正在生成角度 ', num2str(angles_to_plot), ' 的多通道时域对比图...']);

% 2. 时域对比图
% -------------------------------------------------------------------------
figure('Name', '多角度混合响应时域对比', 'Position', [50, 50, 1500, 700]);
sgtitle(['FIR/IIR 混合响应时域对比 (角度: ', num2str(angles_to_plot), ')'], 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:length(channel_names)
    ch_name = channel_names{i};
    ch_label = All_Hybrid_Results.A1.(ch_name).label; % 任何角度的数据都可以用于获取标签

    subplot(2, 2, i);
    hold on;
    h_time_legend = gobjects(num_angles_plot, 1); % 存储图柄用于 legend

    for j = 1:num_angles_plot
        angle_idx = angles_to_plot(j);
        angle_tag = ['A', num2str(angle_idx)];
        
        % 检查角度数据是否存在
        if isfield(All_Hybrid_Results, angle_tag) && isfield(All_Hybrid_Results.(angle_tag), ch_name)
            hybrid_data = All_Hybrid_Results.(angle_tag).(ch_name);
            
            % 绘制时域响应
            h_time_legend(j) = plot(t_samples, hybrid_data.h_hybrid, line_specs{j}, 'LineWidth', 1, ...
                'DisplayName', ['角度 ', num2str(angle_idx)]);
        else
            warning(['通道 ', ch_name, ' 的角度 ', num2str(angle_idx), ' 数据缺失。跳过绘图。']);
        end
    end
    
    % 绘制 80ms 垂直线，标记拼接点
    plot([t_cut_s, t_cut_s], ylim, 'k:', 'LineWidth', 1, 'DisplayName', [num2str(t_cut_ms), 'ms 拼接点']);
    
    title(['通道: ', ch_label], 'FontSize', 12);
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('幅度', 'FontSize', 10);
    xlim([0, max(t_samples)]);
    legend(h_time_legend(isvalid(h_time_legend)), 'Location', 'northeast', 'FontSize', 9);
    grid on;
    hold off; % 结束当前通道的子图
end

disp('所有多角度时域对比图已生成。');
% -------------------------------------------------------------------------