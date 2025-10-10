%% ------------------------------------------------------------------------
% FIR to IIR Conversion (Hybrid 5.1, 2 ears × 6 Tx)
% Output: mIIR_B, mIIR_A as 3D arrays [2,6,nB/nA]
% ------------------------------------------------------------------------
clear; clc;

%% 1. Load FIR Data
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat'); % mIRInt: [samples, ears, speakers, angles]
fs = 44100;                % Sampling frequency
cut_time = 0.08;           % 80 ms
cut_idx = round(cut_time * fs); % 对应样本点索引

num_ears = 2;
num_tx   = 6;              % 5.1声道: 6 Tx
angle_idx = 53;            % 固定角度
%% 2. IIR Design Parameters
B_order = 7500;
A_order = 7500;

%% 3. Preallocate storage as 3D arrays
mIIR_B = zeros(num_ears, num_tx, B_order+1);
mIIR_A = zeros(num_ears, num_tx, A_order+1);

%% 4. Loop over ears and Tx
for iCRx = 1:num_ears
    for iCTx = 1:num_tx
        % Get FIR
        h_fir_full = squeeze(mIRInt(:, iCRx, iCTx, angle_idx));
        
        % Split head and tail
        h_head = h_fir_full(1:cut_idx);
        h_tail = h_fir_full(cut_idx+1:end);
        
        % FFT of tail
        N_fft = length(h_tail);
        H_tail = fft(h_tail);
        f_half = (0:floor(N_fft/2))*(fs/N_fft);
        H_half = H_tail(1:length(f_half));
        omega = 2*pi*f_half/fs; % normalized frequency [0,pi]
        W = ones(size(H_half));  % weights
        
        % invfreqz least-squares IIR
        disp(['Designing IIR for Rx ', num2str(iCRx), ', Tx ', num2str(iCTx)]);
        [B_iir, A_iir] = invfreqz(H_half, omega, B_order, A_order, W);
        
        % Store as 3D array
        mIIR_B(iCRx, iCTx, 1:length(B_iir)) = B_iir;
        mIIR_A(iCRx, iCTx, 1:length(A_iir)) = A_iir;
    end
end

%% 5. Save as .mat for direct use
save('IIR_filters.mat','mIIR_B','mIIR_A');
disp('IIR_filters.mat saved (3D arrays, compatible with real-time script).');
