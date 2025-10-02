% MATLAB Code: Cross-Angle FIR/IIR Hybrid Response (Fixed IIR Tail)
% -------------------------------------------------------------------------
% 目的: 
% 1. 加载角度 53 的 IIR 系数，作为四个通道的锁定 IIR 拖尾源。
% 2. 遍历所有角度，将每个 (通道, 角度) 的 FIR 头部 (< 80ms) 
%    与对应通道锁定的 IIR 尾部 (>= 80ms) 进行拼接，生成混合响应。
% -------------------------------------------------------------------------

% 0. 初始化与数据加载 (与之前保持一致)
% -------------------------------------------------------------------------
% 1. 加载 IIR 锁定系数和参数 (来自角度 53 的拟合结果)
iir_lock_filename = 'WLS_IIR_Coefficients_N7500_M7500_A53_Full_Process.mat';
try
    load(iir_lock_filename);
    
    % 关键参数来自加载的文件
    fs = fs;
    N_order = N_order;
    M_order = M_order;
    t_cut_ms = t_cut_ms;
    idx_cut_start = idx_cut_start;
    
    % 提取锁定 IIR 系数和 t=0 开始的 IIR 冲激响应
    Locked_IIR_Coefficients = IIR_Results; 
    
    channel_names = fieldnames(IIR_Results);
    num_channels = length(channel_names);
    fixed_angle_idx = 53; % 锁定 IIR 尾部的角度索引
    
    disp(['IIR 系数从文件 ', iir_lock_filename, ' 加载成功。']);
    disp(['IIR 拖尾锁定于角度 ', num2str(fixed_angle_idx), '。']);
catch ME
    disp(['错误: 无法加载文件 ', iir_lock_filename, '。请检查文件路径和内容。']);
    return;
end

% 2. 加载原始 FIR 数据 (包含所有角度)
fir_data_file = 'filters/Room_Home_221025_5_1_HP_HD800_221025.mat';
try
    load(fir_data_file);
    N_fir_samples = size(mIRInt, 1); % 脉冲响应的样本总数
    
    % **验证和定义角度总数 (107)**
    num_angles_loaded = size(mIRInt, 4);
    if num_angles_loaded ~= 107
        warning(['警告: 加载的数据文件中的角度总数是 ', num2str(num_angles_loaded), '，而不是预期的 107。']);
    end
    num_angles = num_angles_loaded; % 使用实际加载的值，但已记录警告
    
catch ME
    disp(['错误: 无法加载原始 FIR 数据文件 ', fir_data_file, '。']);
    return;
end

% 定义通道信息 (方便从 mIRInt 中提取数据)
channels_info = {
    {1, 1, '左耳, 扬声器 1', 'Ch11'},
    {1, 2, '左耳, 扬声器 2', 'Ch12'},
    {2, 1, '右耳, 扬声器 1', 'Ch21'},
    {2, 2, '右耳, 扬声器 2', 'Ch22'}
};

% 3. 共享频率轴准备
f_nyquist = fs / 2;
num_points = floor(N_fir_samples/2) + 1;
f_half = linspace(0, f_nyquist, num_points);

disp(['原始数据总角度数: ', num2str(num_angles)]);
disp('-------------------------------------------------------------------------');

% 1. 批处理循环：构建所有角度的混合响应
% -------------------------------------------------------------------------
All_Hybrid_Results = struct(); % 结构体用于存储所有角度的混合响应结果

disp(['开始构建所有 ', num2str(num_angles), ' 个角度的混合响应...']);

% 定义 IIR 拖尾的长度 (用于确保索引不超出范围)
len_iir_part = N_fir_samples - idx_cut_start + 1;

for k = 1:num_angles
    current_angle_idx = k;
    angle_tag = ['A', num2str(current_angle_idx)];
    All_Hybrid_Results.(angle_tag) = struct();

    for i = 1:num_channels
        [ear_idx, spk_idx, label, ch_name] = channels_info{i}{:};
        
        % A. 获取锁定的 IIR 拖尾
        % IIR 冲激响应 h_iir_raw 是从 t=0 开始的
        h_iir_raw = Locked_IIR_Coefficients.(ch_name).h_iir_raw;
        
        % **修正点: 提取 IIR 响应中对应于 t >= 80ms 的信号**
        % IIR 拖尾是从 idx_cut_start 开始到末尾的部分
        h_iir_tail = h_iir_raw(idx_cut_start : N_fir_samples);
        
        % B. 获取当前角度的 FIR 头部
        h_fir_full = squeeze(mIRInt(:, ear_idx, spk_idx, current_angle_idx));
        
        % 提取 FIR 的线性相位头部 (0ms 到 80ms 之前)
        h_fir_pre = h_fir_full(1 : idx_cut_start - 1);
        
        % C. 拼接：创建混合 (Hybrid) 脉冲响应
        % 确保 h_fir_pre 的长度 + h_iir_tail 的长度 = N_fir_samples
        if (length(h_fir_pre) + length(h_iir_tail)) ~= N_fir_samples
            % 这不应该发生，但作为健壮性检查
            error(['通道 ', ch_name, ' 角度 ', num2str(k), ' 的拼接长度不匹配原始长度。']);
        end
        h_hybrid = [h_fir_pre; h_iir_tail];
        
        % D. 频域分析
        H_hybrid_half = fft(h_hybrid);
        H_hybrid_half = H_hybrid_half(1:num_points);
        
        Magnitude_dB_hybrid = 20 * log10(abs(H_hybrid_half));
        Phase_unwrapped_deg_hybrid = rad2deg(unwrap(angle(H_hybrid_half)));
        
        % E. 存储结果
        All_Hybrid_Results.(angle_tag).(ch_name).label = label;
        All_Hybrid_Results.(angle_tag).(ch_name).h_hybrid = h_hybrid;
        All_Hybrid_Results.(angle_tag).(ch_name).Magnitude_dB = Magnitude_dB_hybrid;
        All_Hybrid_Results.(angle_tag).(ch_name).Phase_unwrapped_deg = Phase_unwrapped_deg_hybrid;
    end
end

disp('所有角度的混合响应构建完成。');
disp('-------------------------------------------------------------------------');

% 2. 可视化 (仅绘制用于锁定的角度 53 的结果) (保持原样)
% -------------------------------------------------------------------------
disp(['正在绘制角度 ', num2str(fixed_angle_idx), ' 的时域和频域对比图...']);

% 获取角度 53 的数据
A53_Hybrid_Data = All_Hybrid_Results.(['A', num2str(fixed_angle_idx)]);
t_samples = (0:N_fir_samples-1) / fs;
t_cut_s = t_cut_ms / 1000;
FIR_Color = 'b';
Hybrid_LineSpec = 'r--'; 

% --- 频域对比图 (所有通道) ---
figure('Position', [50, 50, 1600, 800]);
sgtitle(['FIR/IIR 混合响应频域对比 (锁定 IIR: 角度 ', num2str(fixed_angle_idx), ')'], 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:num_channels
    [ear_idx, spk_idx, ~, ch_name] = channels_info{i}{:};
    hybrid_data = A53_Hybrid_Data.(ch_name);
    
    % 获取原始 FIR 响应 (作为对比)
    h_fir_full = squeeze(mIRInt(:, ear_idx, spk_idx, fixed_angle_idx));
    H_fir_half = fft(h_fir_full);
    H_fir_half = H_fir_half(1:num_points);
    Magnitude_dB_fir = 20 * log10(abs(H_fir_half));

    % 子图 1: 幅值响应 (上排)
    subplot(2, 4, i);
    semilogx(f_half, Magnitude_dB_fir, FIR_Color, 'LineWidth', 1, 'DisplayName', '原始 FIR');
    hold on;
    semilogx(f_half, hybrid_data.Magnitude_dB, Hybrid_LineSpec, 'LineWidth', 1.5, 'DisplayName', '混合响应');
    title(['幅度: ', hybrid_data.label], 'FontSize', 12);
    xlabel('频率 (Hz)', 'FontSize', 10);
    ylabel('幅值 (dB)', 'FontSize', 10);
    xlim([20, f_nyquist]);
    legend('Location', 'southwest', 'FontSize', 8);
    grid on;

    % 子图 2: 相位响应 (下排)
    subplot(2, 4, i + 4);
    % 相位需要重新计算
    Phase_unwrapped_deg_fir = rad2deg(unwrap(angle(H_fir_half))); 
    semilogx(f_half, Phase_unwrapped_deg_fir, FIR_Color, 'LineWidth', 1, 'DisplayName', '原始 FIR');
    hold on;
    semilogx(f_half, hybrid_data.Phase_unwrapped_deg, Hybrid_LineSpec, 'LineWidth', 1.5, 'DisplayName', '混合响应');
    title(['相位: ', hybrid_data.label], 'FontSize', 12);
    xlabel('频率 (Hz)', 'FontSize', 10);
    ylabel('相位 (度)', 'FontSize', 10);
    xlim([20, f_nyquist]);
    legend('Location', 'southwest', 'FontSize', 8);
    grid on;
end

% --- 时域对比图 (所有通道) ---
figure('Name', ['FIR/IIR 混合响应 vs. 原始 FIR (角度 ', num2str(fixed_angle_idx), ')'], 'Position', [50, 50, 1500, 700]);
sgtitle(['FIR/IIR 混合响应 vs. 原始 FIR 时域对比 (角度 ', num2str(fixed_angle_idx), ', 拼接点: ', num2str(t_cut_ms), 'ms)'], 'FontSize', 16, 'FontWeight', 'bold');

for i = 1:num_channels
    [ear_idx, spk_idx, ~, ch_name] = channels_info{i}{:};
    hybrid_data = A53_Hybrid_Data.(ch_name);
    
    % 获取原始 FIR 响应
    h_fir_full = squeeze(mIRInt(:, ear_idx, spk_idx, fixed_angle_idx));
    
    subplot(2, 2, i);
    plot(t_samples, h_fir_full, FIR_Color, 'LineWidth', 1, 'DisplayName', ['原始 FIR (', hybrid_data.label, ')']);
    hold on;
    plot(t_samples, hybrid_data.h_hybrid, Hybrid_LineSpec, 'LineWidth', 1.5, 'DisplayName', ['混合 FIR/IIR 响应']);
    
    plot([t_cut_s, t_cut_s], ylim, 'k:', 'LineWidth', 1, 'DisplayName', [num2str(t_cut_ms), 'ms 拼接点']);
    
    hold off;
    title(['通道: ', hybrid_data.label], 'FontSize', 12);
    xlabel('时间 (秒)', 'FontSize', 10);
    ylabel('幅度', 'FontSize', 10);
    xlim([0, max(t_samples)]);
    legend('Location', 'northeast', 'FontSize', 9);
    grid on;
end

% 3. 数据保存 (保持原样)
% -------------------------------------------------------------------------
final_filename = ['Cross_Angle_Hybrid_N', num2str(N_order), '_M', num2str(M_order), '_Locked_A', num2str(fixed_angle_idx), '_Final.mat'];

% 保存锁定的 IIR 系数和所有角度的混合响应结果
save(final_filename, 'Locked_IIR_Coefficients', 'All_Hybrid_Results', 'fs', 'f_half', 't_cut_ms', 'idx_cut_start', 'fixed_angle_idx'); 

disp('-------------------------------------------------------------------------');
disp(['所有角度的混合响应和锁定的 IIR 系数已保存到文件: ', final_filename]);
disp(['已生成角度 ', num2str(fixed_angle_idx), ' 的时域和频域对比图。']);