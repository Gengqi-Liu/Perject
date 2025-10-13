%% ------------------------------------------------------------------------
% FIR to IIR Conversion (Hybrid 5.1, 2 ears × 6 Tx × 101 Angles)
% Output: mIIR_B, mIIR_A as 4D arrays [2,6,101,nB/nA]
% ------------------------------------------------------------------------
clear; clc;

%% 1. Load FIR Data
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat'); % mIRInt: [samples, ears, speakers, angles]
fs = 44100;                % Sampling frequency
cut_time = 0;              % e.g. 80 ms → 0.08
cut_idx = round(cut_time * fs);

num_ears   = 2;
num_tx     = 6;            % 5.1 channels
num_angles = size(mIRInt, 4); % typically 101

%% 2. IIR Design Parameters
B_order = 5000;
A_order = 5000;

%% 3. Preallocate storage as 4D arrays
mIIR_B = zeros(num_ears, num_tx, num_angles, B_order+1);
mIIR_A = zeros(num_ears, num_tx, num_angles, A_order+1);

%% 4. Loop over ears, Tx, and angles
for iCRx = 1:num_ears
    for iCTx = 1:num_tx
        for iAng = 1:num_angles
            % Get FIR
            h_fir_full = squeeze(mIRInt(:, iCRx, iCTx, iAng));

            % Split head and tail
            if cut_idx > 0
                h_head = h_fir_full(1:cut_idx);
                h_tail = h_fir_full(cut_idx+1:end);
            else
                h_head = [];
                h_tail = h_fir_full;
            end

            % FFT of tail
            N_fft = length(h_tail);
            H_tail = fft(h_tail);
            f_half = (0:floor(N_fft/2))*(fs/N_fft);
            H_half = H_tail(1:length(f_half));
            omega = 2*pi*f_half/fs; % normalized frequency [0,pi]
            W = ones(size(H_half));  % equal weights

            % invfreqz least-squares IIR
            disp(['Designing IIR for Rx ', num2str(iCRx), ...
                  ', Tx ', num2str(iCTx), ', Angle ', num2str(iAng), ...
                  ' / ', num2str(num_angles)]);

            try
                [B_iir, A_iir] = invfreqz(H_half, omega, B_order, A_order, W);
            catch ME
                warning(['invfreqz failed at Rx ', num2str(iCRx), ...
                         ', Tx ', num2str(iCTx), ', Angle ', num2str(iAng), ...
                         ': ', ME.message]);
                B_iir = zeros(1,B_order+1);
                A_iir = zeros(1,A_order+1);
            end

            % Store
            mIIR_B(iCRx, iCTx, iAng, 1:length(B_iir)) = B_iir;
            mIIR_A(iCRx, iCTx, iAng, 1:length(A_iir)) = A_iir;
        end
    end
end

%% 5. Save as .mat for direct use
save('IIR_filters_allAngles.mat','mIIR_B','mIIR_A','-v7.3');
disp('IIR_filters_allAngles.mat saved (4D arrays [2×6×101×nCoeff]).');
