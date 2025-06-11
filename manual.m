load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat');  % 自动加载
dims = size(mIRInt);   % 获取维度大小

fprintf('Loaded mIRInt with size: [%d samples, %d ears, %d speakers, %d directions]\n', dims);

% 你可以自定义你要选取的耳/扬声器/角度
ear_id       = 1;       % 1 = 左耳, 2 = 右耳
speaker_id   = 1;       % 1~dims(3)
direction_id = 53;      % 1~dims(4)

% 自动提取该组合的 FIR 响应
h_fir_full = squeeze(mIRInt(:, ear_id, speaker_id, direction_id));

% 限定长度（如前 1000 点）
L = 1000;
h_fir = h_fir_full(1:min(L, length(h_fir_full)));

%% === Step 2: 设置 IIR 阶数 ===
M = 100;  % 分子 b 阶数
N = 100;  % 分母 a 阶数（不含 a0）

%% === Step 3: 构造最小解误差方程组 A x = b ===
row_count = L - N;
col_count = N + M + 1;
A = zeros(row_count, col_count);
b_vec = zeros(row_count, 1);

for i = N+1:L
    % 构造一行
    for j = 1:N
        A(i - N, j) = -h_fir(i - j);
    end
    if (i - N) <= M + 1
        A(i - N, N + (i - N)) = 1;
    end
    b_vec(i - N) = h_fir(i);
end

% 手写最小二乘解（伪逆）
AtA = A' * A;
Atb = A' * b_vec;
x = AtA \ Atb;

% 提取 IIR 系数
a_coeffs = [1; x(1:N)];
b_coeffs = x(N+1:end);

%% === Step 4: 手写计算 IIR 冲激响应 ===
x_impulse = zeros(L, 1);
x_impulse(1) = 1;
h_iir = zeros(L, 1);

for n = 1:L
    % 输入项
    for k = 0:M
        if n - k >= 1
            h_iir(n) = h_iir(n) + b_coeffs(k+1) * x_impulse(n - k);
        end
    end
    % 反馈项
    for k = 1:N
        if n - k >= 1
            h_iir(n) = h_iir(n) - a_coeffs(k+1) * h_iir(n - k);
        end
    end
end

%% === Step 5: 手写计算 DFT 频率响应 ===
K = 1024;
H_fir = zeros(K, 1);
H_iir = zeros(K, 1);

for k = 0:K-1
    omega = 2 * pi * k / K;
    for n = 1:L
        H_fir(k+1) = H_fir(k+1) + h_fir(n) * exp(-1j * omega * (n-1));
        H_iir(k+1) = H_iir(k+1) + h_iir(n) * exp(-1j * omega * (n-1));
    end
end
f = (0:K-1) / K;

%% === Step 6: 绘图比较 ===
figure;
subplot(2,1,1);
plot(h_fir, 'b', 'DisplayName', 'FIR');
hold on;
plot(h_iir, 'r--', 'DisplayName', 'Fitted IIR');
xlabel('Sample Index');
ylabel('Amplitude');
title('Impulse Response Comparison');
legend;
grid on;

subplot(2,1,2);
plot(f, abs(H_fir), 'b', 'DisplayName', 'FIR');
hold on;
plot(f, abs(H_iir), 'r--', 'DisplayName', 'Fitted IIR');
xlabel('Normalized Frequency (\times\pi rad/sample)');
ylabel('Magnitude');
title('Frequency Response Comparison');
legend;
grid on;