% MATLAB 代码：WLS IIR 滤波器逼近与比较 (新增相位对比)
% -------------------------------------------------------------------------
% 0. 设置与数据加载
% -------------------------------------------------------------------------
% 1. 加载已保存的频域数据 (H_half, f_half, fs, Magnitude_dB, Phase_unwrapped_deg)
load('FIR_Frequency_Response_Data.mat'); 
% -------------------------------------------------------------------------
% 1. WLS IIR 滤波器设计 (使用 invfreqz 函数)
% -------------------------------------------------------------------------
% 确定 IIR 滤波器阶数 (N: 分子阶数, M: 分母阶数)
N_order = 7500; % 分子 (B) 系数阶数
M_order = 7500; % 分母 (A) 系数阶数

% invfreqz 输入数据准备
% 归一化频率点 omega (从 0 到 pi)
omega = f_half / (fs / 2) * pi;
% 权重矩阵 W: 标准最小二乘法 (所有权重为 1)
W = ones(size(H_half)); 

disp('正在使用 invfreqz 进行加权最小二乘 IIR 滤波器设计...');
[B_iir, A_iir] = invfreqz(H_half, omega, N_order, M_order, W); 
disp('IIR 滤波器设计完成。');
% -------------------------------------------------------------------------
% 2. IIR 频率响应计算和数据准备
% -------------------------------------------------------------------------
% 在相同的频率点上计算 IIR 滤波器的频率响应
% H_iir_complex 是复数响应
[H_iir_complex, ~] = freqz(B_iir, A_iir, omega); 

% 计算 IIR 逼近的幅值响应 (dB)
Magnitude_iir_dB = 20 * log10(abs(H_iir_complex));

% ***************************************************************
% 新增：计算 IIR 逼近的解卷绕相位响应 (度)
% ***************************************************************
% MATLAB 的 angle 函数返回的是卷绕相位，我们需要进行解卷绕以进行清晰对比
Phase_iir_rad = angle(H_iir_complex);
% unwrap 函数执行相位解卷绕
Phase_iir_unwrapped_deg = unwrap(Phase_iir_rad) * (180/pi); 
% -------------------------------------------------------------------------
% 3. 可视化比较 (绘图比较)
% -------------------------------------------------------------------------
figure('Position', [100, 100, 1000, 800]); 

% --- 子图 1: 幅值响应比较 (顶部) ---
subplot(2, 1, 1);
semilogx(f_half, Magnitude_dB, 'b', 'LineWidth', 1.5, 'DisplayName', '目标 FIR 响应');
hold on;
semilogx(f_half, Magnitude_iir_dB, 'r--', 'LineWidth', 1.5, 'DisplayName', 'WLS IIR 逼近');
hold off;
title('幅值响应比较：FIR 目标 vs. WLS IIR 逼近', 'FontSize', 14);
xlabel('频率 (Hz)', 'FontSize', 12);
ylabel('幅值 (dB)', 'FontSize', 12);
xlim([20, fs/2]); 
legend('Location', 'southwest');
grid on;

% --- 子图 2: 相位响应比较 (底部) ---
subplot(2, 1, 2);
% 绘制原始 FIR 解卷绕相位响应 (蓝色)
semilogx(f_half, Phase_unwrapped_deg, 'b', 'LineWidth', 1.5, 'DisplayName', '目标 FIR 线性相位');
hold on;
% 绘制 IIR 逼近解卷绕相位响应 (红色)
semilogx(f_half, Phase_iir_unwrapped_deg, 'r--', 'LineWidth', 1.5, 'DisplayName', 'WLS IIR 相位');
hold off;
title('相位响应比较：FIR 目标 (线性) vs. WLS IIR (非线性)', 'FontSize', 14);
xlabel('频率 (Hz)', 'FontSize', 12);
ylabel('相位 (度)', 'FontSize', 12);
xlim([20, fs/2]); 
legend('Location', 'southwest');
grid on;

% -------------------------------------------------------------------------
% 4. 保存结果
% -------------------------------------------------------------------------
% 定义保存的文件名
iir_filename = 'WLS_IIR_Filter_Coefficients.mat'; 
% 保存 IIR 滤波器系数 B 和 A
save(iir_filename, 'B_iir', 'A_iir', 'N_order', 'M_order', 'fs'); 
disp(['IIR 滤波器系数已成功保存到文件: ', iir_filename]);
disp(['IIR 分子阶数 N: ', num2str(N_order), ', 分母阶数 M: ', num2str(M_order)]);
