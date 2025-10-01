% MATLAB Code: IIR vs. FIR Time Domain Impulse Response Comparison
% -------------------------------------------------------------------------
% 0. 初始化参数与数据加载
% -------------------------------------------------------------------------

% 设定两个滤波器的目标采样频率
fs_target = 44100; % Hz

% 1. 加载 IIR 滤波器系数 (B_iir, A_iir)
iir_filename = 'WLS_IIR_Filter_Coefficients.mat'; 
try
    load(iir_filename);
    if ~exist('B_iir', 'var') || ~exist('A_iir', 'var')
        error('加载的文件中缺少 B_iir 或 A_iir 系数。');
    end
catch ME
    disp(['错误: 无法加载文件 ', iir_filename, '，请确保该文件存在并包含 B_iir 和 A_iir。']);
    return; % 终止程序
end

% 2. 加载原始 FIR 数据 (作为比较的目标)
fir_data_file = 'filters/Room_Home_221025_5_1_HP_HD800_221025.mat';
try
    load(fir_data_file); 
    % 提取特定的冲激响应序列: 左耳(1), 扬声器 1(1), 角度索引 53
    h_fir_full = squeeze(mIRInt(:, 1, 1, 53));  
catch ME
    disp(['错误: 无法加载或解析原始 FIR 数据文件 ', fir_data_file, '。']);
    return; % 终止程序
end

% 确定原始 FIR 的长度和持续时间 (以进行精确对比)
N_fir_samples = length(h_fir_full); 
duration_sec = N_fir_samples / fs_target;

disp(['原始 FIR 长度: ', num2str(N_fir_samples), ' 个样本。']);
disp(['对比持续时间: ', num2str(duration_sec), ' 秒 (基于 FIR 目标)。']);

% -------------------------------------------------------------------------
% 1. 计算 IIR 滤波器冲激响应
% -------------------------------------------------------------------------

% 使用 impz(B, A, N, fs) 计算 IIR 滤波器在相同长度 (N_fir_samples) 下的冲激响应
% N 必须匹配原始 FIR 的长度 (N_fir_samples)
disp('正在计算 IIR 逼近的冲激响应...');
[h_iir_time, t_samples] = impz(B_iir, A_iir, N_fir_samples, fs_target);
disp('IIR 冲激响应计算完成。');

% -------------------------------------------------------------------------
% 2. 可视化对比 (时域)
% -------------------------------------------------------------------------

figure('Name', 'IIR vs. FIR 时域冲激响应对比', 'Position', [100, 100, 800, 600]); 

% 绘制对比图
plot(t_samples, h_fir_full, 'b', 'LineWidth', 1.5, 'DisplayName', '目标 FIR 响应 (线性相位)');
hold on;
plot(t_samples, h_iir_time, 'r--', 'LineWidth', 1, 'DisplayName', 'WLS IIR 逼近 (非线性相位)');
hold off;

title('IIR vs. FIR 时域冲激响应对比 (完整持续时间)', 'FontSize', 14);
xlabel('时间 (秒)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
legend('Location', 'northeast');
grid on;

% -------------------------------------------------------------------------
% 3. (可选) 细节放大图 (观察波形差异)
% -------------------------------------------------------------------------

figure('Name', 'IIR vs. FIR 时域响应细节对比', 'Position', [950, 100, 800, 600]); 

plot(t_samples, h_fir_full, 'b', 'LineWidth', 1.5, 'DisplayName', '目标 FIR 响应');
hold on;
plot(t_samples, h_iir_time, 'r--', 'LineWidth', 1, 'DisplayName', 'WLS IIR 逼近');
hold off;

% 限制 X 轴到一个较短的时间窗口 (例如，前 50 毫秒) 来查看瞬态细节
% 这对于观察 FIR 和 IIR 在波形上的时间对齐差异非常有用
detail_limit_sec = 0.05; 
xlim([0, detail_limit_sec]); 
ylim_range = max(abs([h_fir_full(t_samples<detail_limit_sec); h_iir_time(t_samples<detail_limit_sec)]));
ylim([-ylim_range*1.1, ylim_range*1.1]); % 自动调整 Y 轴范围

title(['IIR vs. FIR 冲激响应细节对比 (前 ', num2str(detail_limit_sec*1000), ' 毫秒)'], 'FontSize', 14);
xlabel('时间 (秒)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
legend('Location', 'northwest');
grid on
