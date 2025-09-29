% MATLAB Code for FIR Frequency-Domain to Time-Domain Transformation (Comparison)
% -------------------------------------------------------------------------
% 0. Setup and Data Loading
% -------------------------------------------------------------------------

% 1. 加载上一步保存的频域数据 (H_half, fs, 等)
load('FIR_Frequency_Response_Data.mat'); 

% 2. 重新加载原始的 FIR 数据 (用于对比)
% ！！！注意：这里需要确保您的 filters/ 路径和原始文件与您第一次运行的代码一致
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat');
full_fir = squeeze(mIRInt(:, 1, 1, 53));  

% 检查原始 IR 的长度 N_original
N_original = length(full_fir);

% -------------------------------------------------------------------------
% 1. 频域对称性重构 (Reconstruct Full Spectrum)
% -------------------------------------------------------------------------

% 单边谱 H_half 包含 DC (H(1)) 到 Nyquist (H(end))
N_half = length(H_half);
N_full = 2 * (N_half - 1); % IFFT需要的点数 N_full (本例中 N=63000)

% --- 重构 H_reconstructed (与原 IFFT 脚本一致) ---
if mod(N_full, 2) == 0 
    H_conjugate = conj(flipud(H_half(2:end-1)));
else 
    H_conjugate = conj(flipud(H_half(2:end)));
end
H_reconstructed = [H_half; H_conjugate];

% -------------------------------------------------------------------------
% 2. 逆快速傅里叶变换 (IFFT)
% -------------------------------------------------------------------------

h_ifft = ifft(H_reconstructed);
h_reconstructed_fir = real(h_ifft); 

% -------------------------------------------------------------------------
% 3. 可视化对比 (Plotting Comparison)
% -------------------------------------------------------------------------

% 创建时间轴 (秒)
t = (0:N_full-1) / fs; 

figure('Position', [100, 100, 1200, 500]); 

% --- 子图 1: 完整的 IR 对比 ---
subplot(1, 2, 1);
% 绘制原始 IR (蓝色实线)
plot(t, full_fir, 'b', 'LineWidth', 1.5, 'DisplayName', '原始 FIR (full\_fir)');
hold on;
% 绘制重构 IR (红色虚线)
plot(t, h_reconstructed_fir, 'r--', 'LineWidth', 1, 'DisplayName', '重构 FIR (IFFT 结果)');
hold off;

title('完整脉冲响应对比：原始 vs. 重构', 'FontSize', 14);
xlabel('时间 (秒)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);
legend('Location', 'northeast');
grid on;

% --- 子图 2: 放大 IR 的关键起始部分 (用于观察精度) ---
subplot(1, 2, 2);
plot(t, full_fir, 'b', 'LineWidth', 1.5, 'DisplayName', '原始 FIR');
hold on;
plot(t, h_reconstructed_fir, 'r--', 'LineWidth', 1, 'DisplayName', '重构 FIR');
hold off;

title('脉冲响应起始部分放大对比', 'FontSize', 14);
xlabel('时间 (秒)', 'FontSize', 12);
ylabel('幅度', 'FontSize', 12);

% 将 X 轴限制在前 50 毫秒 (0.05 秒) 以清晰对比起始波形
xlim([0, 0.05]); 
legend('Location', 'northwest');
grid on;

disp('已完成原始和重构脉冲响应的可视化对比。');