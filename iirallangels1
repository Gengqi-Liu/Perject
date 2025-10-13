%% ------------------------------------------------------------------------
% FIR to IIR Conversion (Hybrid 5.1, 2 ears × 6 Tx × 101 Angles)
% Auto-resume + incremental save version
% ------------------------------------------------------------------------
clear; clc;

%% 1. Load FIR Data
load('filters/Room_Home_221025_5_1_HP_HD800_221025.mat'); % mIRInt: [samples, ears, speakers, angles]
fs = 44100;                  % Sampling frequency
cut_time = 0;                % 头部截断长度 (s)
cut_idx = round(cut_time * fs);

num_ears   = 2;
num_tx     = 6;
num_angles = size(mIRInt, 4); % 通常为101

%% 2. IIR Design Parameters
B_order = 5000;
A_order = 5000;

%% 3. 文件名
save_file = 'IIR_filters_allAngles.mat';
temp_file = [save_file '.tmp'];

%% 4. 如果已有文件，则续算
if isfile(save_file)
    load(save_file, 'mIIR_B', 'mIIR_A', 'doneMask');
    disp('✅ 已加载已有进度，将从上次未完成的角度继续...');
else
    % 初始化
    mIIR_B = zeros(num_ears, num_tx, num_angles, B_order+1);
    mIIR_A = zeros(num_ears, num_tx, num_angles, A_order+1);
    doneMask = false(num_ears, num_tx, num_angles); % 已完成标志
end

%% 5. 并行池设置（可选）
if isempty(gcp('nocreate'))
    parpool('threads');  % 或 parpool('local')，依机器情况选择
end

%% 6. 主循环：耳 × Tx × 角度
for iCRx = 1:num_ears
    for iCTx = 1:num_tx

        parfor iAng = 1:num_angles  % 可换为 for 以调试
            if doneMask(iCRx, iCTx, iAng)
                continue; % 已完成，跳过
            end

            try
                % --- 获取FIR ---
                h_fir_full = squeeze(mIRInt(:, iCRx, iCTx, iAng));
                if cut_idx > 0
                    h_tail = h_fir_full(cut_idx+1:end);
                else
                    h_tail = h_fir_full;
                end

                % --- 频域数据 ---
                N_fft = length(h_tail);
                H_tail = fft(h_tail);
                f_half = (0:floor(N_fft/2))*(fs/N_fft);
                H_half = H_tail(1:length(f_half));
                omega = 2*pi*f_half/fs;
                W = ones(size(H_half));

                % --- IIR 拟合 ---
                [B_iir, A_iir] = invfreqz(H_half, omega, B_order, A_order, W);

                % --- 存入临时变量 ---
                localB = zeros(1, B_order+1);
                localA = zeros(1, A_order+1);
                localB(1:length(B_iir)) = B_iir;
                localA(1:length(A_iir)) = A_iir;

                % --- 写入结果 (并行安全) ---
                parsave(temp_file, iCRx, iCTx, iAng, localB, localA);

            catch ME
                warning(['❌ invfreqz failed at Rx=',num2str(iCRx),...
                    ', Tx=',num2str(iCTx),', Ang=',num2str(iAng),...
                    ': ',ME.message]);
            end
        end
    end
end

disp('✅ 全部IIR设计完成！');
