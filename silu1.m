% III. 可视化 II：时域冲激响应对比 (直接对比，红蓝主题)
% -------------------------------------------------------------------------
% 目的：直接对比原始完整 FIR 响应和 IIR 逼近滤波器自身的冲激响应（均从 t=0 开始），
%       使用统一的红蓝配色。
% -------------------------------------------------------------------------

% 确保所有变量（如 fs, N_fir_samples, IIR_Results, Results, t_cut_ms）
% 在加载文件后都可用。
load('WLS_IIR_Coefficients_N7500_M7500_A53_Full_Process.mat');

% 重新定义 channel_names 和 num_channels，以防文件加载后未定义
channel_names = fieldnames(IIR_Results);
num_channels = length(channel_names);

disp('开始生成时域对比图：IIR 滤波器自身响应 vs 原始 FIR 响应...');

figure('Name', 'IIR vs. 原始 FIR 时域响应直接对比 (红蓝)', 'Position', [50, 50, 1500, 700]);
sgtitle('IIR 滤波器冲激响应 vs. 原始 FIR 响应（均从 t=0 开始）', 'FontSize', 16, 'FontWeight', 'bold');

t_samples = (0:N_fir_samples-1) / fs; % 时间轴

% 定义统一的颜色和线型
FIR_Color = 'b';
IIR_LineSpec = 'r--'; 

for i = 1:num_channels
    ch_name = channel_names{i};
    fir_data = Results.(ch_name);
    iir_data = IIR_Results.(ch_name);
    
    h_fir_full = fir_data.h_fir_full; % 原始完整 FIR 响应
    h_iir_raw = iir_data.h_iir_raw;   % IIR 从 t=0 开始的冲激响应

    % A. 绘制对比图
    subplot(2, 2, i);
    
    % 绘制原始完整 FIR 响应 (蓝色实线)
    plot(t_samples, h_fir_full, FIR_Color, 'LineWidth', 1, 'DisplayName', ['原始 FIR (', fir_data.label, ')']);
    hold on;
    
    % 绘制 IIR 逼近响应（红色虚线，从 t=0 开始，无任何偏移）
    plot(t_samples, h_iir_raw, IIR_LineSpec, 'LineWidth', 1.5, 'DisplayName', ['WLS IIR 响应']);
    
    % 绘制 80ms 垂直线，标记原始 FIR 的线性部分和拖尾部分的分界点
    t_cut_s = t_cut_ms / 1000;
    plot([t_cut_s, t_cut_s], ylim, 'k:', 'LineWidth', 1, 'DisplayName', [num2str(t_cut_ms), 'ms 标记']);
    
    hold off;
    title(['通道: ', fir_data.label], 'FontSize', 12);
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('幅度', 'FontSize', 10);
    xlim([0, max(t_samples)]); % 显示完整时长
    legend('Location', 'northeast', 'FontSize', 9);
    grid on;
end

disp('所有时域直接对比图已生成。');
% -------------------------------------------------------------------------