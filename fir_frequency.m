% MATLAB Code for FIR Time-Domain to Frequency-Domain Transformation

% -------------------------------------------------------------------------
% 0. Setup and Data Loading
% -------------------------------------------------------------------------

% 1. Load the .mat file (Ensure the file path is correct)
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat');

% 2. Select the specific Impulse Response (IR)
%    Selection: Left Ear (1), Speaker 1 (1), Angle Index 53
h_fir_full = squeeze(mIRInt(:, 1, 1, 53));  

% Original Data Parameters
fs = 48000;           % Sampling Frequency (Hz)
N = length(h_fir_full); % Number of Samples (63000)

% -------------------------------------------------------------------------
% 1. Frequency Domain Transformation (DFT/FFT)
% -------------------------------------------------------------------------

% Perform the Discrete Fourier Transform (FFT)
H_complex = fft(h_fir_full);

% -------------------------------------------------------------------------
% 2. Prepare Plotting Data: Select Positive Frequency Axis (0 to Nyquist)
% -------------------------------------------------------------------------

% Calculate the Nyquist Frequency (fs/2)
f_nyquist = fs / 2;

% Determine the number of points for the single-sided spectrum (DC to Nyquist)
% Since N=63000 is even, Nyquist is at N/2 + 1
num_points = floor(N/2) + 1; 

% Create the Frequency Axis (Hz): from 0 to f_nyquist
f_half = linspace(0, f_nyquist, num_points); 

% Keep only the complex response for the positive frequencies
H_half = H_complex(1:num_points);

% -------------------------------------------------------------------------
% 3. Magnitude Response Calculation
% -------------------------------------------------------------------------

% Convert magnitude to Decibels (dB)
Magnitude_dB = 20 * log10(abs(H_half));

% -------------------------------------------------------------------------
% 4. Phase Response Calculation
% -------------------------------------------------------------------------

% Calculate phase in radians
Phase_rad = angle(H_half);

% Unwrap the phase and convert to degrees to remove +/- 180 degree jumps
Phase_unwrapped_deg = rad2deg(unwrap(Phase_rad)); 

% -------------------------------------------------------------------------
% 5. Visualization (Plotting)
% -------------------------------------------------------------------------

% Create a figure window with a specified size
figure('Position', [100, 100, 1200, 500]); 

% --- Plot Magnitude Response (Left Subplot) ---
subplot(1, 2, 1);
% Use semilogx for a logarithmic frequency axis (standard in acoustics)
semilogx(f_half, Magnitude_dB, 'LineWidth', 1.5);
title('FIR Filter Magnitude Response', 'FontSize', 14);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Magnitude (dB)', 'FontSize', 12);

% Set x-axis limit to start from 20 Hz to ignore DC/very low-frequency components
xlim([20, f_nyquist]); 
grid on;
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off'); % Keep grid clean

% --- Plot Unwrapped Phase Response (Right Subplot) ---
subplot(1, 2, 2);
% Logarithmic frequency axis
semilogx(f_half, Phase_unwrapped_deg, 'LineWidth', 1.5);
title('FIR Filter Unwrapped Phase Response', 'FontSize', 14);
xlabel('Frequency (Hz)', 'FontSize', 12);
ylabel('Phase (Degrees)', 'FontSize', 12);

% Match frequency range with the magnitude plot
xlim([20, f_nyquist]); 
grid on;
set(gca, 'XMinorGrid', 'off', 'YMinorGrid', 'off');

% -------------------------------------------------------------------------
% 6. Data Saving (保存结果)
% -------------------------------------------------------------------------
% 定义保存的文件名
filename = 'FIR_Frequency_Response_Data.mat'; 

% 保存核心频域数据 (频率轴和复数响应)
% 同时保存幅频和相频结果，以备直接调用
save(filename, 'f_half', 'H_half', 'Magnitude_dB', 'Phase_unwrapped_deg', 'fs'); 

disp(['所有频域结果已成功保存到文件: ', filename]);
disp('您可以随时使用 load() 函数调用这些数据。');