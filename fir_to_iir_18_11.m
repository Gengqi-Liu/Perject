%% ------------------------------------------------------------------------
% FIR to IIR Conversion (Hybrid 5.1, 2 ears × 6 Tx)
% 输出: IIR_filters.mat  (mIIR_B, mIIR_A, cut_idx, tIIRDelay)
% ------------------------------------------------------------------------
clear; clc;

%% 1. Load FIR Data
% mIRInt: [samples, ears, speakers, angles]
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat'); 

fs       = 44100;      % Sampling frequency
tIIRDelay = 0.08;      % 80 ms
cut_idx  = round(tIIRDelay * fs); % 对应样本点索引

num_ears = 2;
num_tx   = 6;          % 5.1声道: 6 Tx
angle_idx = 53;        % 固定中间角度 (103 个角里的中点)

%% 2. IIR Design Parameters
B_order = 7500;
A_order = 7500;

%% 3. Preallocate storage as 3D arrays
mIIR_B = zeros(num_ears, num_tx, B_order+1);
mIIR_A = zeros(num_ears, num_tx, A_order+1);

%% 4. Loop over ears and Tx
for iCRx = 1:num_ears
    for iCTx = 1:num_tx
        % 对应耳朵 & Tx 的 FIR (固定 angle_idx)
        h_fir_full = squeeze(mIRInt(:, iCRx, iCTx, angle_idx));
        
        % 防止长度不足 cut_idx
        if length(h_fir_full) <= cut_idx+10
            error('IR 太短: Rx=%d, Tx=%d, length=%d', ...
                  iCRx, iCTx, length(h_fir_full));
        end
        
        % Split head and tail
        h_head = h_fir_full(1:cut_idx);          %#ok<NASGU> % 这里只是说明用法
        h_tail = h_fir_full(cut_idx+1:end);      % 80ms 之后的尾部

        % 对 tail 单独做 FFT (去掉前 80ms 延迟)
        N_fft  = length(h_tail);
        H_tail = fft(h_tail);
        f_half = (0:floor(N_fft/2))*(fs/N_fft);
        H_half = H_tail(1:length(f_half));
        
        omega  = 2*pi*f_half/fs;   % normalized frequency [0,pi]
        W      = ones(size(H_half));  % weights
        
        % invfreqz least-squares IIR
        fprintf('Designing IIR for Rx %d, Tx %d ...\n', iCRx, iCTx);
        [B_iir, A_iir] = invfreqz(H_half, omega, B_order, A_order, W);
        
        % 防止 invfreqz 失败
        if any(isnan(B_iir)) || any(isnan(A_iir))
            error('invfreqz failed: NaN in coefficients (Rx=%d, Tx=%d).', ...
                  iCRx, iCTx);
        end
        
        % Store as 3D array
        mIIR_B(iCRx, iCTx, 1:length(B_iir)) = B_iir;
        mIIR_A(iCRx, iCTx, 1:length(A_iir)) = A_iir;
    end
end

%% 5. Save as .mat for direct use
save('IIR_filters.mat','mIIR_B','mIIR_A','cut_idx','tIIRDelay');
disp('IIR_filters.mat saved (mIIR_B, mIIR_A, cut_idx, tIIRDelay).');