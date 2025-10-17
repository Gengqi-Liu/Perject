%%10.17这需要使用sos的脚本。

addpath('utils');
SetParameters;

% --- 新增：定义关键参数 (假设 frameLength 在 SetParameters 中定义) ---
% 假设 SetParameters 或外部定义了 frameLength
if ~exist('frameLength', 'var')
    frameLength = 1024; % 默认值，如果未定义
end
Fs = 44100; % 假设采样率
iNoRx = 2; % 接收通道数量（耳）
% -------------------------------------------------------------------------


%% Initialize IIR Filters (for each ear × transmitter pair)
% -------------------------------------------------------------------------
% Example placeholder file; replace later with your real IIR coefficients
if exist('IIR_filters.mat','file')
  load('IIR_filters.mat','mIIR_B','mIIR_A');
  disp('Loaded IIR filters (mIIR_B, mIIR_A)');
else
  % Placeholder: simple lowpass-type IIR (for testing)
  [b0,a0] = butter(2, 0.4); % 2nd-order, normalized cutoff 0.4
  mIIR_B = repmat(reshape(b0,1,1,[]),[2,6,1]); % 2 Rx × 6 Tx
  mIIR_A = repmat(reshape(a0,1,1,[]),[2,6,1]);
  disp('Using placeholder IIR filters (butterworth 2nd order)');
end

% --- 修改：使用 dsp.BiquadFilter 系统对象 (推荐) ---
% IIR 滤波器系统对象 (2 Rx * 6 Tx = 12 个滤波器)
% 我们使用 dsp.BiquadFilter 以获得更好的数值稳定性
IIRFilters = cell(iNoRx, iNoTx); % iNoRx=2, iNoTx=6 (假设)

for iCRx = 1:iNoRx
    for iCTx = 1:iNoTx
        b = squeeze(mIIR_B(iCRx,iCTx,:));
        a = squeeze(mIIR_A(iCRx,iCTx,:));
        
        % 转换为二阶分段 (SOS) 形式
        [sos, g] = tf2sos(b, a);
        
        IIRFilters{iCRx, iCTx} = dsp.BiquadFilter('SOSMatrix', sos, ...
                                                  'ScaleValues', g);
    end
end
% 原始的 mIIRReg 不再需要，但为保持变量兼容性，保留其初始化
maxOrder = max(size(mIIR_A,3), size(mIIR_B,3)) - 1;
mIIRReg = zeros(maxOrder,iNoRx,6);  % [order, iCRx, iCTx]


% --- IIR 尾部连接所需的延迟缓冲区初始化 (80ms) ---
Delay_sec = 0.080;          % 80 ms
Delay_samples = round(Delay_sec * Fs); % 3528 个采样点 (假设 Fs=44100)

% IIR 输入延迟缓冲区 (用于实现 80ms 后的尾部连接)
% 维度：[总大小 x iNoTx]
% 这里 iNoTx 必须在 SetParameters 或 NUPOLS 初始化前被定义
if ~exist('iNoTx', 'var')
    % 假设 iNoTx = 6，如果 NUPOLS 初始化前还未定义
    iNoTx = 6;
end
DelayBuffer = zeros(Delay_samples + frameLength, iNoTx);
DelayBufferPointer = 1; % 环形缓冲区指针
% ----------------------------------------------------------------------


%% Initialize NUPOLS


% only for bBLEAudio on Windows -> 48kHz
bBLEAudio = false;
if ispc
  if bBLEAudio
    % 假设这里将 Fs 更改为 48000
    mIRInt48 = zeros(round(size(mIRInt,1)*48/44.1),size(mIRInt,2),size(mIRInt,3),size(mIRInt,4));
    for iCRx = 1:2
      for iCTx = 1:iNoTx
        for iCA = 1:size(mIRInt,4)
          mIRInt48(:,iCRx,iCTx,iCA) = SRC(mIRInt(:,iCRx,iCTx,iCA),44.1e3,48e3);
        end
      end
    end
    mIRInt = mIRInt48;
    Fs = 48000; % 更新 Fs
    disp('Resampled to 48 kHz, required for WASAPI');
  end
end

iNoTx = size(mIRInt,3);
iNoAngleEl = size(mIRInt,5);
vTxInd = 1:iNoTx;
[H1tp,mFDL_buf1_old,mFDL_buf1_cur,x_in_buf1_old,x_in_buf1_cur,iC1,...
  H2tp,mFDL_buf2_old,mFDL_buf2_cur,x_in_buf2_old,x_in_buf2_cur,iC2,...
  H3tp,mFDL_buf3_old,mFDL_buf3_cur,x_in_buf3_old,x_in_buf3_cur,iC3,...
  H4tp,mFDL_buf4_old,mFDL_buf4_cur,x_in_buf4_old,x_in_buf4_cur,iC4,...
  x_ring,y_ring_cur,y_ring_old,B,Bi,Pi,N_RB_x,N_RB_y,N,...
  vBlockDelay,vSchedOffset]...
  = InitializeNUPOLS(frameLength,mIRInt);
H1 = H1tp(:,:,:,:,:,(iNoAngleEl+1)/2);
H2 = H2tp(:,:,:,:,:,(iNoAngleEl+1)/2);
H3 = H3tp(:,:,:,:,:,(iNoAngleEl+1)/2);
H4 = H4tp(:,:,:,:,:,(iNoAngleEl+1)/2);
if iNoAngleEl > 1
  disp('ELevation correction active');
else
  disp('Elevation correction inactive');
end
% clear mIRInt;


%% Crossover filter
bTrinaural = false;
crossFilt = crossoverFilter('NumCrossovers',1,'CrossoverFrequencies',5000, ...
    'CrossoverSlopes',12);
[b1,a1,b2,a2] = getFilterCoefficients(crossFilt,1);


%% Initialize UDP receiver for headtracker
if exist('u','var')
  clear u
end
echoudp("off")
if ispc
  echoudp("on",5005)
else % macOS
  echoudp("on",5006)
end
u = udpport("datagram",'LocalHost','127.0.0.1','LocalPort',5005);


% --- 新增：初始化 Audio System Objects ---
% 替换 fileReader 和 deviceWriterActive/Speaker/Headphone
audioReader = audioDeviceReader('SampleRate', Fs, 'SamplesPerFrame', frameLength, 'NumChannels', iNoTx); 
deviceWriterHeadphone = audioDeviceWriter('SampleRate', Fs, 'ChannelMapping', [1 2]); % 假设立体声输出
deviceWriterSpeaker = audioDeviceWriter('SampleRate', Fs, 'ChannelMapping', 1:iNoTx); % 假设多声道输出

% 初始激活的设备
deviceWriterActive = deviceWriterHeadphone;

% --------------------------------------


%% Initialize counters etc.
iCount        = 0;
fAngleHor     = 0;
fAngleVer     = 0;
fAngleHorCal  = 0; iCalCount = 0;
vAngleHorSave = zeros(1,1e5,'single');
vAngleHorSave2 = zeros(1,1e5,'single');
vAngleHorPred = zeros(1,1e5,'single');
bHeadphone    = true; % start with headphone ON
iCountMax     = 0;
fMaxAmpl      = 0;
iRunTimeLen   = 100000;
vRunTime      = zeros(1,iRunTimeLen,'single');
vUnderrun     = false(1,iRunTimeLen);

% For speaker implementation
iDelay        = 50;
vIR           = [zeros(1,iDelay+1),1];
mReg          = zeros(length(vIR)-1,iNoTx);
% For fading implementation: Output is calculated twice (old/current) and
% mixed
vWeightsUp    = (1:frameLength).'/frameLength;
mWeightsUp    = repmat(vWeightsUp,1,iNoRx); % 确保维度正确
mWeightsDown  = 1-mWeightsUp;

%%

%%Fsos = dsp.SOSFilter(sosCoefficients(:,1:3),sosCoefficients(:,4:6));

%% Real-time processing
disp('Real-time convolving starts ... ')
while true % endless

  iCount = iCount + 1;

  %% Receiving headtracker data via UDP
  if u.NumDatagramsAvailable > 0
    data = read(u,u.NumDatagramsAvailable,"char");
    sAngleNew = data(end).Data;
    fAngleHor = str2double(sAngleNew(1:8));
    fAngleVer = str2double(sAngleNew(9:16));
  end 

  % Kalman filtering of raw data will be added here later
  vAngleHorPred(mod(iCount-1,length(vAngleHorPred))+1) = fAngleHor;

  bAutoCal = false;
  if bAutoCal
    % automatic center calibration
    fOldAngleHorWithCal = vAngleHorPred(mod(iCount-2,length(vAngleHorPred))+1);
    fNewAngleHorWithCal = fAngleHor-fAngleHorCal;
    if abs(fOldAngleHorWithCal-fNewAngleHorWithCal) < 5 && abs(fNewAngleHorWithCal) > 5
      iCalCount = iCalCount + 1;
    else
      iCalCount = 0;
    end
    if iCalCount > round(5*Fs/frameLength) % 假设 fSamplFreq = Fs
      fAngleHorCal = fAngleHor;
      iCalCount = 0;
    end
     
    vAngleHorPred(mod(iCount-1,length(vAngleHorPred))+1) = mod(fAngleHor-fAngleHorCal-180,360)-180;
  end

  %% Toggle between headphone and loudspeakers
  if mod(iCount,5)==1
    if bHeadphone
      if fAngleVer<-50
        bHeadphone = false;
        release(deviceWriterActive);
        deviceWriterActive = deviceWriterSpeaker;
        mReg   = zeros(length(vIR)-1,iNoTx);
        mOut   = zeros(frameLength,iNoTx);
      end
    else % Speakers are on
      if fAngleVer>-55
        bHeadphone = true;
        release(deviceWriterActive);
        deviceWriterActive = deviceWriterHeadphone;
      end
    end
  end

  %% Read data
  % --- 修改：使用 audioDeviceReader 读取数据 ---
  mIn       = audioReader();

  % crossover filter

  if iNoTx > 2
    if mod(iCount,100)==1
      if exist('status/Trinaural.txt','file')
        fid = fopen('status/Trinaural.txt', 'r');
        phi_degree = fscanf(fid,'%f');
        if phi_degree < 90 && phi_degree >= 0
          phi = phi_degree*pi/180;
          bTrinaural = true;
        else
          bTrinaural = false;
        end
        fclose(fid);
      else
        bTrinaural = false;
      end
    end
  end

  % trinaural synthesis
  if bTrinaural

    w = 1;
     
    a = 1/2*(sin(phi) + w);
    b = 1/2*(sin(phi) - w);
    c = 1/sqrt(2)*cos(phi);
    P = [a,b;b,a;c,c];
    mIn(:,1:3) = (P*mIn(:,1:2)')';
    bTwoFreqBands = false;
    if bTwoFreqBands
      phi2 = atan(sqrt(2)); % Gerzon: 54.74°
      a2 = 1/2*(sin(phi2) + w);
      b2 = 1/2*(sin(phi2) - w);
      c2 = 1/sqrt(2)*cos(phi2);
      P2 = [a2,b2;b2,a2;c2,c2];
      mIn2(:,1:3) = (P2*mIn2(:,1:2)')';
      % sum
      mIn = mIn + mIn2;
    end
  end


  %% Check whether headphone is on
  tic; % for performance analysis
  if ~bHeadphone
    %% Headphone is OFF
    if bShowDisplay && mod(iCount,iNoIterShowDisplay)==1
      PrintStatus(sRoomName,sHeadphoneName,'OFF');
    end
    for iCTx=1:iNoTx
      [mOut(:,iCTx),mReg(:,iCTx)] = filter(vIR,1,mIn(:,iCTx),mReg(:,iCTx));
    end
  else
    %% Headphone is ON
    if bShowDisplay && mod(iCount,iNoIterShowDisplay)==1
      PrintStatus(sRoomName,sHeadphoneName,'ON',N,iNoTx,fAngleHor-fAngleHorCal,fAngleVer,...
        vAngle,vRunTime,vUnderrun,fUpdateTime,iCount);
    end

    %% Frequency-domain real-time convolution
    iCMod                 = 1+mod(iCount-1,length(vAngleHorSave));
    vAngleHorSave(iCMod)  = fAngleHor; % save angle for debugging
    
    [~,iAngleIndOld]      = min(abs(vAngle-vAngleHorPred(1+mod(iCMod-2,length(vAngleHorSave)))));
    [~,iAngleIndCur]      = min(abs(vAngle-vAngleHorPred(iCMod)));

    %% Update ring buffer
    iCircPt_x             = mod(iCount-1,N_RB_x)+1;
    iCircPt_y             = mod(iCount-1,N_RB_y)+1;
    x_ring(:,iCircPt_x,:) = mIn;
    iCircPt_y_Del                 = mod(iCircPt_y-1-1,N_RB_y)+1;
    y_ring_old(:,iCircPt_y_Del,:) = 0;
    y_ring_cur(:,iCircPt_y_Del,:) = 0;

    %% SEGMENT 1
    if mod(iCount,Bi(1)/B==0) % when it is available, here every block
      iC1 = iC1 + 1;
      vInd1_x = iCircPt_x;
      mIn_seg = x_ring(:,vInd1_x,:);
      vIndUpdate = [iAngleIndOld,iAngleIndCur]; % will be called in any iteration
      % H1(:,:,:,:,vIndUpdate) = interpElevation(H1tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer); 
      [y_part1_old,x_in_buf1_old,mFDL_buf1_old] = UPConv(mIn_seg,x_in_buf1_old,mFDL_buf1_old,vTxInd,H1,iC1,iAngleIndOld,Bi(1),Pi(1));
      [y_part1_cur,x_in_buf1_cur,mFDL_buf1_cur] = UPConv(mIn_seg,x_in_buf1_cur,mFDL_buf1_cur,vTxInd,H1,iC1,iAngleIndCur,Bi(1),Pi(1));
      vInd1_y = iCircPt_y;
      y_ring_old(:,vInd1_y,:) = y_ring_old(:,vInd1_y,:) + reshape(y_part1_old,B,Bi(1)/B,[]);
      y_ring_cur(:,vInd1_y,:) = y_ring_cur(:,vInd1_y,:) + reshape(y_part1_cur,B,Bi(1)/B,[]);
    end
    %% SEGMENT 2
    if numel(Bi)>1 && mod(iCount,Bi(2)/B)==0 % when it is available
      iC2 = iC2 + 1;
      vInd2_x = mod(iCircPt_x+(-Bi(2)/B+1:0)-1,N_RB_x)+1;
      mIn_seg = reshape(x_ring(:,vInd2_x,:),Bi(2),[]);
      vIndUpdate = [iAngleIndOld,iAngleIndCur];
      % H2(:,:,:,:,vIndUpdate) = interpElevation(H2tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer); 
      [y_part2_old,x_in_buf2_old,mFDL_buf2_old] = UPConv(mIn_seg,x_in_buf2_old,mFDL_buf2_old,vTxInd,H2,iC2,iAngleIndOld,Bi(2),Pi(2));
      [y_part2_cur,x_in_buf2_cur,mFDL_buf2_cur] = UPConv(mIn_seg,x_in_buf2_cur,mFDL_buf2_cur,vTxInd,H2,iC2,iAngleIndCur,Bi(2),Pi(2));
      vInd2_y = mod(iCircPt_y+vBlockDelay(2)+(0:Bi(2)/B-1)-1,N_RB_y)+1;
      y_ring_old(:,vInd2_y,:) = y_ring_old(:,vInd2_y,:) + reshape(y_part2_old,B,Bi(2)/B,[]);
      y_ring_cur(:,vInd2_y,:) = y_ring_cur(:,vInd2_y,:) + reshape(y_part2_cur,B,Bi(2)/B,[]);
    end
    %% SEGMENT 3
    if numel(Bi)>2 && mod(iCount,Bi(3)/B)==0+vSchedOffset(3) % when it is available
      iC3 = iC3 + 1;
      vInd3_x = mod(iCircPt_x-vSchedOffset(3)+(-Bi(3)/B+1:0)-1,N_RB_x)+1;
      mIn_seg = reshape(x_ring(:,vInd3_x,:),Bi(3),[]);
      vIndUpdate = [iAngleIndOld,iAngleIndCur];
      % H3(:,:,:,:,vIndUpdate) = interpElevation(H3tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer); 
      [y_part3_old,x_in_buf3_old,mFDL_buf3_old]  = UPConv(mIn_seg,x_in_buf3_old,mFDL_buf3_old,vTxInd,H3,iC3,iAngleIndOld,Bi(3),Pi(3));
      [y_part3_cur,x_in_buf3_cur,mFDL_buf3_cur]  = UPConv(mIn_seg,x_in_buf3_cur,mFDL_buf3_cur,vTxInd,H3,iC3,iAngleIndCur,Bi(3),Pi(3));
      vInd3_y = mod(iCircPt_y-vSchedOffset(3)+vBlockDelay(3)+(0:Bi(3)/B-1)-1,N_RB_y)+1;
      y_ring_old(:,vInd3_y,:) = y_ring_old(:,vInd3_y,:) + reshape(y_part3_old,B,Bi(3)/B,[]);
      y_ring_cur(:,vInd3_y,:) = y_ring_cur(:,vInd3_y,:) + reshape(y_part3_cur,B,Bi(3)/B,[]);
    end
    %% SEGMENT 4
    if numel(Bi)>=4 && mod(iCount,Bi(4)/B)==0+vSchedOffset(4) % when it is available
      iC4 = iC4 + 1;
      vInd4_x = mod(iCircPt_x-vSchedOffset(4)+(-Bi(4)/B+1:0)-1,N_RB_x)+1;
      mIn_seg     = reshape(x_ring(:,vInd4_x,:),Bi(4),[]);
      % H4(:,:,:,:,vIndUpdate) = interpElevation(H4tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer);       
      [y_part4_old,x_in_buf4_old,mFDL_buf4_old]  = UPConv(mIn_seg,x_in_buf4_old,mFDL_buf4_old,vTxInd,H4,iC4,iAngleIndOld,Bi(4),Pi(4));
      [y_part4_cur,x_in_buf4_cur,mFDL_buf4_cur]  = UPConv(mIn_seg,x_in_buf4_cur,mFDL_buf4_cur,vTxInd,H4,iC4,iAngleIndCur,Bi(4),Pi(4));
      vInd4_y = mod(iCircPt_y-vSchedOffset(4)+vBlockDelay(4)+(0:Bi(4)/B-1)-1,N_RB_y)+1;
      y_ring_old(:,vInd4_y,:) = y_ring_old(:,vInd4_y,:) + reshape(y_part4_old,B,Bi(4)/B,[]);
      y_ring_cur(:,vInd4_y,:) = y_ring_cur(:,vInd4_y,:) + reshape(y_part4_cur,B,Bi(4)/B,[]);
    end

    %% Take block from ring buffer (FIR OUTPUT)
    mOut_Old  = squeeze(y_ring_old(:,iCircPt_y,:));
    mOut_Cur  = squeeze(y_ring_cur(:,iCircPt_y,:));
    % FIR 滤波结果：Combine old and current output
    mOut_FIR  = mWeightsDown.*mOut_Old + mWeightsUp.*mOut_Cur;


    % ----------------------------------------------------------------------
    %% APPLY IIR POST-FILTERING (FIR → IIR HYBRID STAGE)
    % IIR 作为 80ms 后的尾部，与 FIR 输出并行相加。
    
    % 1. 从延迟缓冲区取出 80ms 后的输入
    startIdx = DelayBufferPointer;
    endIdx = startIdx + frameLength - 1;
    mOut_Delayed = DelayBuffer(startIdx:endIdx, :); % [frameLength x iNoTx]

    % 2. IIR 滤波块 (使用延迟后的 mOut_Delayed)
    mOut_IIR = zeros(size(mOut_FIR)); % IIR 输出
    for iCRx = 1:iNoRx
      for iCTx = 1:iNoTx
        % --- 修改：使用 dsp.BiquadFilter 系统对象 ---
        % 系统对象自动管理状态
        mOut_IIR_temp = IIRFilters{iCRx, iCTx}(mOut_Delayed(:,iCTx));
        % --------------------------------------------
        
        % IIR 输出累加 (保持原有的累加逻辑)
        mOut_IIR(:,iCTx) = mOut_IIR(:,iCTx) + mOut_IIR_temp; 
      end
    end

    % 3. 更新延迟缓冲区 (将当前 FIR 的输入 mIn 存入)
    DelayBuffer(DelayBufferPointer:endIdx, :) = mIn;
    
    % 更新指针 (环形缓冲区)
    DelayBufferPointer = DelayBufferPointer + frameLength;
    if DelayBufferPointer > Delay_samples
        DelayBufferPointer = DelayBufferPointer - Delay_samples;
    end

    % 4. 最终合并：FIR (前 80ms) + IIR (后 80ms 尾部)
    mOut = mOut_FIR + mOut_IIR; 
    % ----------------------------------------------------------------------


%targetRMS = 0.1; % 可调
%currentRMS = rms(mOut(:));
%mOut = mOut * (targetRMS/currentRMS);
    %% headphone equalization
    if 0%bHPEQ
      for iCRx=1:2
        % 注意：如果这里也要改用系统对象，则需要初始化相应的 HPEQ 滤波器对象
        [mOut(:,iCRx),mHPReg(:,iCRx)] = filter(vHPIR,1,mOut(:,iCRx),mHPReg(:,iCRx));
      end
    end

    % %% Nicos equalization
    % ... (保持不变) ...


    % Signalverarbeitung (z. B. Block mit N Samples)
    for iCRx = 1:2
%         mOut(:,iCRx) = Fsos(msOut(:,iCRx));  % Filterung mit Zustand
    end

    %% Preamplifier
    mOut      = 256*mOut;
    fMaxAmpl  = max(max(abs(mOut(:))),fMaxAmpl);
    if bShowDisplay && mod(iCount,iNoIterShowDisplay)==1
      if fMaxAmpl>0.5
        iCountMax = iCountMax + 1;
        if iCountMax == 3000
          fMaxAmpl = 0;
          iCountMax = 0;
        end
        disp(['Critical maximal amplitude: ',num2str(fMaxAmpl)]);
      end
    end
    % save run time per iteration for speed analysis
    vRunTime(mod(iCount-1,length(vRunTime))+1) = toc;
  end
  %% Write data to output buffer
  % --- 修改：使用 audioDeviceWriter 写入数据 ---
  nUnderrun = deviceWriterActive(mOut(:,1:iNoRx)); % 假设 headphone/speaker 只有 iNoRx (2) 个通道
  vUnderrun(mod(iCount-1,length(vRunTime))+1) = false;
  if nUnderrun > 0
    fprintf('Audio writer queue was underrun by %d samples.\n',...
      nUnderrun);
    vUnderrun(mod(iCount-1,length(vRunTime))+1) = true;
  end
end
