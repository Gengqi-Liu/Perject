% MATLAB Code for WLS IIR Filter Approximation and Comparison
% -------------------------------------------------------------------------
% 0. Setup and Data Loading
% -------------------------------------------------------------------------

% 1. 加载保存的频域数据 (H_half, f_half, fs, Magnitude_dB, Phase_unwrapped_deg)
load('FIR_Frequency_Response_Data.mat'); 

% -------------------------------------------------------------------------
% 1. WLS IIR 滤波器设计 (使用 invfreqz)
% -------------------------------------------------------------------------

% 确定 IIR 滤波器的阶数 (N: 分子阶数, M: 分母阶数)
% 选择合适的阶数对逼近质量至关重要。阶数越高，逼近越精确，但滤波器越复杂。
N_order = 10000; % 分子 (B) 系数阶数 (通常用于逼近FIR的长脉冲响应)
M_order = 10000; % 分母 (A) 系数阶数 (用于创建极点，实现更高效的逼近)

% invfreqz 输入数据准备
% H：目标复数频响 (H_half 的实部和虚部)
% W：频率权重 (设置为1，即标准最小二乘，或加权)
% F：频率点 (f_half，但需要归一化到 [0, pi] 弧度)

% 归一化频率点 omega (从 0 到 pi)
omega = f_half / (fs / 2) * pi;

% 权重矩阵 W：我们采用标准最小二乘法 (权重全为 1) 作为起点。
% W 的长度必须与 H_half 和 omega 相同。
W = ones(size(H_half)); 

% !!! 修正后的 invfreqz 调用 !!!
% invfreqz(H, W_rad, N, M, Wt)
disp('正在使用 invfreqz 进行加权最小二乘 IIR 滤波器设计...');
[B_iir, A_iir] = invfreqz(H_half, omega, N_order, M_order, W); 
disp('IIR 滤波器设计完成。');

% -------------------------------------------------------------------------
% 2. IIR 频响计算与对比数据准备
% -------------------------------------------------------------------------

% 计算 IIR 滤波器在相同频率点上的频率响应
% freqz(B, A, W) 使用 W (弧度) 作为频率点
H_iir_complex = freqz(B_iir, A_iir, omega);

% 计算 IIR 逼近的幅频响应 (dB)
Magnitude_iir_dB = 20 * log10(abs(H_iir_complex));

% -------------------------------------------------------------------------
% 3. 可视化对比 (Plotting Comparison)
% -------------------------------------------------------------------------

figure('Position', [100, 100, 1200, 500]); 

% --- 子图 1: 幅频响应对比 (Magnitude Response) ---
subplot(1, 2, 1);
% 绘制原始 FIR 的幅频响应 (蓝色)
semilogx(f_half, Magnitude_dB, 'b', 'LineWidth', 1.5, 'DisplayName', '目标 FIR 响应');
hold on;
% 绘制 IIR 逼近的幅频响应 (红色)
semilogx(f_half, Magnitude_iir_dB, 'r--', 'LineWidth', 1.5, 'DisplayName', 'WLS IIR 逼近');
hold off;

title('幅频响应对比：FIR 目标 vs. WLS IIR 逼近', 'FontSize', 14);
xlabel('频率 (Hz)', 'FontSize', 12);
ylabel('幅度 (dB)', 'FontSize', 12);
xlim([20, fs/2]); % 从 20 Hz 到 奈奎斯特频率
legend('Location', 'southwest');
grid on;

% -------------------------------------------------------------------------
% 4. 结果保存
% -------------------------------------------------------------------------

% 定义保存的文件名
%iir_filename = 'WLS_IIR_Filter_Coefficients.mat'; 

% 保存 IIR 滤波器的系数 B 和 A
%save(iir_filename, 'B_iir', 'A_iir', 'N_order', 'M_order', 'fs'); 

%disp(['IIR 滤波器系数已成功保存到文件: ', iir_filename]);
%disp(['IIR 分子阶数 N: ', num2str(N_order), ', 分母阶数 M: ', num2str(M_order)]);