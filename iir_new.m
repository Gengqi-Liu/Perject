%% ------------------------------------------------------------------------
% Generate IIR_filters.mat compatible with hybrid FIR+IIR real-time script
% ------------------------------------------------------------------------
clear; clc;

%% 1. Load FIR Data
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat'); % mIRInt: [samples, ears, speakers, angles]

fs = 44100;                % Sampling frequency
cut_time = 0.08;           % 80 ms cutoff for tail
cut_idx = round(cut_time * fs);

[num_samples, num_ears, num_speakers, num_angles] = size(mIRInt);

%% 2. IIR Design Parameters
B_order = 7500; % 分子阶数
A_order = 7500; % 分母阶数

%% 3. Preallocate mIIR struct array
mIIR = cell(num_ears, num_speakers);

%% 4. Loop over ears and speakers
for ear = 1:num_ears
    for spk = 1:num_speakers
        % 使用角度索引 53
        angle_idx = 53;
        
        % 完整 FIR
        h_fir_full = squeeze(mIRInt(:, ear, spk, angle_idx));
        
        % ---------------------------
        % 切断尾部
        h_fir_head = h_fir_full(1:cut_idx);       % 保留前80ms
        h_fir_tail = h_fir_full(cut_idx+1:end);   % 尾部拟合
        
        % ---------------------------
        % FFT
        N_fft = length(h_fir_tail);
        H_tail = fft(h_fir_tail);
        f_half = (0:floor(N_fft/2))*(fs/N_fft);
        H_half = H_tail(1:length(f_half));
        
        omega = 2*pi*f_half/fs;
        W = ones(size(H_half));
        
        % ---------------------------
        % invfreqz 拟合 IIR
        disp(['Designing IIR for Ear ', num2str(ear), ' Speaker ', num2str(spk)]);
        [B_iir, A_iir] = invfreqz(H_half, omega, B_order, A_order, W);
        
        % ---------------------------
        % 保存到 struct
        mIIR{ear, spk} = struct('b', B_iir, 'a', A_iir);
    end
end

%% 5. Save
save('IIR_filters.mat', 'mIIR');

disp('IIR_filters.mat generated successfully!');
