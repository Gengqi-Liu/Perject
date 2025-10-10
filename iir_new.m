%% ------------------------------------------------------------------------
% 10.10最新版本获得5.1声道iirFIR to IIR Conversion (Tail Fit via Least-Squares / invfreqz)
% ------------------------------------------------------------------------
clear; clc;

%% 1. Load FIR Data
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat'); % mIRInt: [samples, ears, speakers, angles]

fs = 44100;                % Sampling frequency
cut_time = 0.08;           % 80 ms
cut_idx = round(cut_time * fs); % 对应样本点索引

% FIR dimensions
[num_samples, num_ears, num_speakers, num_angles] = size(mIRInt);

%% 2. IIR Design Parameters
B_order = 7500; % 分子阶数
A_order = 7500; % 分母阶数

%% 3. Preallocate Storage
mIIR_B = cell(num_ears, num_speakers); 
mIIR_A = cell(num_ears, num_speakers);

%% 4. Loop Over Each Ear and Speaker
for ear = 1:num_ears
    for spk = 1:num_speakers
        % 选择角度索引，这里假设用角度 53
        angle_idx = 53;
        
        % 获取完整 FIR
        h_fir_full = squeeze(mIRInt(:, ear, spk, angle_idx));  
        
        % ---------------------------
        % 切断尾部
        h_fir_head = h_fir_full(1:cut_idx);    % 前80ms保留原样
        h_fir_tail = h_fir_full(cut_idx+1:end);% 尾部进行 FFT
        
        % ---------------------------
        % FFT 准备
        N_fft = length(h_fir_tail);
        H_tail = fft(h_fir_tail);             % 复数频率响应
        f_half = (0:floor(N_fft/2))*(fs/N_fft);
        H_half = H_tail(1:length(f_half));
        
        omega = 2*pi*f_half/fs;               % 归一化角频率 [0, pi]
        W = ones(size(H_half));               % 权重
        
        % ---------------------------
        % invfreqz 最小二乘法拟合 IIR
        disp(['Designing IIR for Ear ', num2str(ear), ' Speaker ', num2str(spk), ' ...']);
        [B_iir, A_iir] = invfreqz(H_half, omega, B_order, A_order, W);
        
        % ---------------------------
        % 构造完整 IIR 响应，从0开始
        h_iir_full = filter(B_iir, A_iir, [h_fir_head; zeros(length(h_fir_tail),1)]);
        
        % 可选：验证拟合效果
        % freqz(B_iir, A_iir, 1024, fs);
        % figure; plot((0:length(h_fir_full)-1)/fs, 20*log10(abs(h_fir_full)), 'b'); hold on;
        % plot((0:length(h_iir_full)-1)/fs, 20*log10(abs(h_iir_full)), 'r--'); grid on;
        % legend('FIR', 'IIR');
        
        % ---------------------------
        % 保存 IIR 系数
        mIIR_B{ear, spk} = B_iir;
        mIIR_A{ear, spk} = A_iir;
    end
end

%% 5. Save IIR Filters
save('IIR_filters.mat', 'mIIR_B', 'mIIR_A');

disp('All IIR filters designed and saved successfully.');