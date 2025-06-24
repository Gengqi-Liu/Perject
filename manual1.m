% ========================================
% 高精度 IIR 拟合脚本：将 FIR 响应拟合为“几乎一致”的 IIR
% 使用频域 Chebyshev + IRLS + 极点半径控制
% ========================================

clear; clc;

%% === 1. 加载 FIR 响应数据 ===
load('fir_response.mat');  % 确保 y_data 存在
y = y_data(:);              % 保证是列向量
N = length(y);
disp(['加载 FIR 长度：' num2str(N)]);

%% === 2. 设置频域目标 ===
Nfft = 4096;
H = fft(y, Nfft);
D = H(1:Nfft/2);
om = linspace(0, pi, Nfft/2).';
W = ones(size(om));             % 统一权重
W(om < pi/4) = 5;               % 强化低频拟合
W(om > 3*pi/4) = 0.5;           % 弱化高频影响

%% === 3. 设置 IIR 拟合参数 ===
nb = 100;       % 分子阶数
na = 100;       % 分母阶数
rmax = 0.99;    % 极点最大半径

%% === 4. 初始化 ===
a = ones(na+1, 1);  % 初始分母系数
E = exp(-1j * om * (1:na));

%% === 5. IRLS + 极点控制拟合循环 ===
for iter = 1:25
    % 构造系数矩阵 A：D(k)*E + exp(jωk)
    A = zeros(length(om), na + nb + 1);
    for k = 1:length(om)
        A(k, 1:na) = D(k) * E(k, :);
        A(k, na+1:end) = exp(-1j * om(k) * (0:nb));
    end

    % IRLS 加权解
    h = (A .* W) \ (D .* W);  % 最小加权误差

    % 更新 a, b
    a = [1; h(1:na)];
    b = h(na+1:end);

    % 极点稳定化
    rts = roots(a);
    unstable = abs(rts) > rmax;
    rts(unstable) = rmax * rts(unstable) ./ abs(rts(unstable));
    a = real(poly(rts)).';
end

%% === 6. 响应对比绘图 ===
[Hw, w] = freqz(b, a, Nfft);

figure;
plot(w(1:Nfft/2)/pi, 20*log10(abs(H(1:Nfft/2)) + 1e-6), 'k', 'LineWidth', 1.2); hold on;
plot(w/pi, 20*log10(abs(Hw) + 1e-6), 'r--', 'LineWidth', 1.2);
xlabel('归一化频率 (\times\pi rad/sample)');
ylabel('幅度 (dB)');
legend('原始 FIR', '拟合 IIR');
title('FIR 与 IIR 频响拟合对比（Chebyshev IRLS 拟合）');
grid on;

%% === 7. 可选：保存结果
save('iir_fitted_result.mat', 'b', 'a');
