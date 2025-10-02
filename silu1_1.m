% MATLAB 完整代码（精简版）：FIR 截断、WLS IIR 拟合与时域/频域对比 (四通道)
% -------------------------------------------------------------------------
% 思路: 1. 截取原始 FIR 80ms 后的拖尾。
%       2. 对拖尾部分进行 7500 阶 WLS IIR 滤波器拟合。
%       3. 绘制频域响应对比图。
%       4. 计算 IIR 冲激响应并与原始 FIR 拼接对比（时域）。
% -------------------------------------------------------------------------

% 0. 初始化与数据加载
% -------------------------------------------------------------------------
% 1. 加载原始 FIR 数据文件
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat');

% 原始数据参数
fs = 44100;           % 采样频率 (Hz)
N_fir_samples = size(mIRInt, 1); % 脉冲响应的样本总数
angle_idx = 53;       % 固定角度索引

% 定义截断时间 (80ms)
t_cut_ms = 80;
t_cut_s = t_cut_ms / 1000;
idx_cut_start = floor(t_cut_s * fs) + 1; % 80ms 后第一个样本的索引

% 定义 IIR 滤波器阶数
N_order = 5000; % 分子 (B) 系数阶数
M_order = 5000; % 分母 (A) 系数阶数

disp(['原始 FIR 长度: ', num2str(N_fir_samples), ' 个样本。']);
disp(['IIR 拟合的截断起始点: ', num2str(t_cut_ms), ' ms (索引 ', num2str(idx_cut_start), ')']);
disp(['IIR 目标阶数: N=', num2str(N_order), ', M=', num2str(M_order)]);
disp('-------------------------------------------------------------------------');

% 2. 通道定义与初始化
% -------------------------------------------------------------------------
channels_to_process = {
    {1, 1, '左耳, 扬声器 1'},
    {1, 2, '左耳, 扬声器 2'},
    {2, 1, '右耳, 扬声器 1'},
    {2, 2, '右耳, 扬声器 2'}
};
num_channels = length(channels_to_process);

% 初始化结果存储结构体
Results = struct();
IIR_Results = struct();
comparison_data = cell(num_channels, 1);

% 3. 共享频率轴准备 (只需计算一次)
% -------------------------------------------------------------------------
f_nyquist = fs / 2;
num_points = floor(N_fir_samples/2) + 1;
f_half = linspace(0, f_nyquist, num_points);
omega = f_half / f_nyquist * pi; % 归一化频率点 (0 到 pi)
W = ones(size(f_half)); % WLS 权重 (标准最小二乘法)

% -------------------------------------------------------------------------
% I. FIR 截断/FFT，WLS IIR 拟合与频域分析
% -------------------------------------------------------------------------

disp('开始进行 FIR 截断、FFT 和 IIR 拟合...');

for i = 1:num_channels
    % 提取通道信息
    ear_idx = channels_to_process{i}{1};
    spk_idx = channels_to_process{i}{2};
    label = channels_to_process{i}{3};
    ch_name = ['Ch', num2str(ear_idx), num2str(spk_idx)];
    
    % A. FIR 预处理 (截断 80ms 之前，零填充)
    h_fir_full = squeeze(mIRInt(:, ear_idx, spk_idx, angle_idx));
    h_fir_truncated = h_fir_full(idx_cut_start:N_fir_samples);
    h_fir_processed = [zeros(idx_cut_start - 1, 1); h_fir_truncated];
    
    % B. FFT & 频域准备 (FIR 目标响应)
    H_half_fir = fft(h_fir_processed);
    H_half_fir = H_half_fir(1:num_points);
    Magnitude_dB_fir = 20 * log10(abs(H_half_fir));
    Phase_unwrapped_deg_fir = rad2deg(unwrap(angle(H_half_fir)));
    
    % C. WLS IIR 滤波器设计
    [B_iir, A_iir] = invfreqz(H_half_fir, omega, N_order, M_order, W);
    
    % D. IIR 频域响应计算
    [H_iir_complex, ~] = freqz(B_iir, A_iir, omega);
    Magnitude_iir_dB = 20 * log10(abs(H_iir_complex));
    Phase_iir_unwrapped_deg = rad2deg(unwrap(angle(H_iir_complex)));
    
    % E. IIR 时域响应计算 (用于时域对比)
    h_iir_raw = impz(B_iir, A_iir, N_fir_samples, fs);
    
    % F. 存储 IIR 结果
    IIR_Results.(ch_name).label = label;
    IIR_Results.(ch_name).B_iir = B_iir;
    IIR_Results.(ch_name).A_iir = A_iir;
    IIR_Results.(ch_name).Magnitude_iir_dB = Magnitude_iir_dB;
    IIR_Results.(ch_name).Phase_iir_unwrapped_deg = Phase_iir_unwrapped_deg;
    IIR_Results.(ch_name).h_iir_raw = h_iir_raw; % 存储原始IIR时域响应

    % G. 存储 FIR 结果
    Results.(ch_name).label = label;
    Results.(ch_name).Magnitude_dB = Magnitude_dB_fir;
    Results.(ch_name).Phase_unwrapped_deg = Phase_unwrapped_deg_fir;
    Results.(ch_name).h_fir_full = h_fir_full; % 存储原始完整FIR
    
end

disp('IIR 滤波器设计完成，开始生成频域对比图...');

% -------------------------------------------------------------------------
% II. 可视化 I：频域响应对比
% -------------------------------------------------------------------------
figure('Position', [50, 50, 1600, 800]);
sgtitle(['FIR 目标 (80ms后截断) vs. WLS IIR 逼近 (N=', num2str(N_order), ', M=', num2str(M_order), ')'], 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:num_channels
    ch_name = channel_names{i};
    fir_data = Results.(ch_name);
    iir_data = IIR_Results.(ch_name);
    
    % --- 子图 1: 幅值响应比较 (上排) ---
    subplot(2, 4, i);
    semilogx(f_half, fir_data.Magnitude_dB, 'b', 'LineWidth', 1.5, 'DisplayName', '目标 FIR 响应');
    hold on;
    semilogx(f_half, iir_data.Magnitude_iir_dB, 'r--', 'LineWidth', 1.5, 'DisplayName', 'WLS IIR 逼近');
    hold off;
    title(['幅度: ', fir_data.label], 'FontSize', 12);
    xlabel('频率 (Hz)', 'FontSize', 10);
    ylabel('幅值 (dB)', 'FontSize', 10);
    xlim([20, f_nyquist]);
    legend('Location', 'southwest', 'FontSize', 8);
    grid on;

    % --- 子图 2: 相位响应比较 (下排) ---
    subplot(2, 4, i + 4);
    semilogx(f_half, fir_data.Phase_unwrapped_deg, 'b', 'LineWidth', 1.5, 'DisplayName', '目标 FIR 相位');
    hold on;
    semilogx(f_half, iir_data.Phase_iir_unwrapped_deg, 'r--', 'LineWidth', 1.5, 'DisplayName', 'WLS IIR 相位');
    hold off;
    title(['相位: ', fir_data.label], 'FontSize', 12);
    xlabel('频率 (Hz)', 'FontSize', 10);
    ylabel('相位 (度)', 'FontSize', 10);
    xlim([20, f_nyquist]);
    legend('Location', 'southwest', 'FontSize', 8);
    grid on;
end

% -------------------------------------------------------------------------
% III. 可视化 II：时域冲激响应对比 (拼接)
% -------------------------------------------------------------------------

disp('开始生成时域对比图...');

figure('Name', 'IIR vs. FIR 时域冲激响应对比 (四通道)', 'Position', [50, 50, 1500, 700]);
sgtitle(['IIR vs. 原始 FIR 时域对比 (IIR 对齐到 ', num2str(t_cut_ms), 'ms)'], 'FontSize', 16, 'FontWeight', 'bold');

t_samples = (0:N_fir_samples-1) / fs; % 时间轴

line_colors = {'b', 'r', 'g', 'm'}; % 为四个通道指定颜色

for i = 1:num_channels
    ch_name = channel_names{i};
    fir_data = Results.(ch_name);
    iir_data = IIR_Results.(ch_name);
    
    % A. 时域对齐 (IIR 响应从 idx_cut_start 开始)
    h_iir_raw = iir_data.h_iir_raw;
   % h_iir_aligned = zeros(N_fir_samples, 1);
    
    % IIR 响应应该被放置的部分长度
    len_iir_part = N_fir_samples - idx_cut_start + 1;
    h_iir_part = h_iir_raw(1:len_iir_part);

    % 将其放置到 h_iir_aligned 的正确位置
    h_iir_aligned(idx_cut_start:N_fir_samples) = h_iir_part;

    % B. 绘制对比图
    subplot(2, 2, i);
    % 绘制原始完整 FIR 响应
    plot(t_samples, fir_data.h_fir_full, line_colors{i}, 'LineWidth', 1, 'DisplayName', ['原始完整 FIR (', fir_data.label, ')']);
    hold on;
    % 绘制对齐后的 IIR 逼近响应
    plot(t_samples, h_iir_aligned, [line_colors{i}, '--'], 'LineWidth', 1.5, 'DisplayName', ['WLS IIR 逼近 (对齐)']);
    
    % 绘制 80ms 垂直线
    plot([t_cut_s, t_cut_s], ylim, 'k:', 'LineWidth', 1, 'DisplayName', [num2str(t_cut_ms), 'ms 拼接点']);
    
    hold off;
    title(['通道: ', fir_data.label], 'FontSize', 12);
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('幅度', 'FontSize', 10);
    xlim([0, max(t_samples)]); % 显示完整时长
    legend('Location', 'northeast', 'FontSize', 9);
    grid on;
end

% -------------------------------------------------------------------------
% IV. 数据保存
% -------------------------------------------------------------------------
iir_filename = ['WLS_IIR_Coefficients_N', num2str(N_order), '_M', num2str(M_order), '_A', num2str(angle_idx), '_Full_Process.mat'];

% 保存 IIR 结果、FIR 目标信息、设计参数和频率轴
save(iir_filename, 'IIR_Results', 'Results', 'N_order', 'M_order', 'fs', 'f_half', 't_cut_ms', 'idx_cut_start'); 

disp('-------------------------------------------------------------------------');
disp(['所有结果已统一保存到文件: ', iir_filename]);
disp('所有时域和频域对比图已生成。');
% -------------------------------------------------------------------------