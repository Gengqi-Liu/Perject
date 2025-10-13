%10.10全新版本，调整了逻辑且使用5.1声道iir
addpath('utils');
SetParameters;
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

% Initialize filter memory for all 12 filters
maxOrder = max(size(mIIR_A,3), size(mIIR_B,3)) - 1;
mIIRReg = zeros(maxOrder,2,6);  % [order, iCRx, iCTx]

%% Initialize NUPOLS


% only for bBLEAudio on Windows -> 48kHz
bBLEAudio = false;
if ispc
  if bBLEAudio
    mIRInt48 = zeros(round(size(mIRInt,1)*48/44.1),size(mIRInt,2),size(mIRInt,3),size(mIRInt,4));
    for iCRx = 1:2
      for iCTx = 1:iNoTx
        for iCA = 1:size(mIRInt,4)
          mIRInt48(:,iCRx,iCTx,iCA) = SRC(mIRInt(:,iCRx,iCTx,iCA),44.1e3,48e3);
        end
      end
    end
    mIRInt = mIRInt48;
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


%% Initialize counters etc.
iCount        = 0;
fAngleHor     = 0;
fAngleVer     = 0;
fAngleHorCal  = 0; iCalCount = 0;
vAngleHorSave = zeros(1,1e5,'single');
vAngleHorSave2 = zeros(1,1e5,'single');
vAngleHorPred = zeros(1,1e5,'single');
bHeadphone    = true; % start with headphone ON
iCountMax     = 0;
fMaxAmpl      = 0;
iRunTimeLen   = 100000;
vRunTime      = zeros(1,iRunTimeLen,'single');
vUnderrun     = false(1,iRunTimeLen);

% For speaker implementation
iDelay        = 50;
vIR           = [zeros(1,iDelay+1),1];
mReg          = zeros(length(vIR)-1,iNoTx);
% For fading implementation: Output is calculated twice (old/current) and
% mixed
vWeightsUp    = (1:frameLength).'/frameLength;
mWeightsUp    = repmat(vWeightsUp,1,2);
mWeightsDown  = 1-mWeightsUp;

%%

%%

Fsos = dsp.SOSFilter(sosCoefficients(:,1:3),sosCoefficients(:,4:6));

disp('Real-time IIR filtering starts ... ')
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

  %% Toggle between headphone and loudspeakers
  if mod(iCount,5)==1
    if bHeadphone
      if fAngleVer<-50
        bHeadphone = false;
        release(deviceWriterActive);
        deviceWriterActive = deviceWriterSpeaker;
      end
    else % Speakers are on
      if fAngleVer>-55
        bHeadphone = true;
        release(deviceWriterActive);
        deviceWriterActive = deviceWriterHeadphone;
      end
    end
  end

  %% Read data from file or device
  mIn = fileReader();

  % ---- 禁用 FIR 分段卷积部分 ----
  % 不再进行 x_ring / y_ring / H1..H4 / UPConv 等计算
  % 直接将输入 mIn 送入 IIR 滤波器

  %% ---- IIR 处理部分 ----
  % 此部分保留并成为唯一的滤波操作
  mOut_IIR = zeros(size(mIn)); % 输出大小与输入一致

  for iCRx = 1:2
    for iCTx = 1:iNoTx
      [mOut_IIR(:,iCTx), mIIRReg(:,iCRx,iCTx)] = ...
          filter(squeeze(mIIR_B(iCRx,iCTx,:)), squeeze(mIIR_A(iCRx,iCTx,:)), ...
                 mIn(:,iCTx), mIIRReg(:,iCRx,iCTx));
    end
  end

  % 将滤波输出设为主输出
  mOut = mOut_IIR;

  %% Headphone equalization（可选，默认关闭）
  if 0%bHPEQ
    for iCRx=1:2
      [mOut(:,iCRx),mHPReg(:,iCRx)] = filter(vHPIR,1,mOut(:,iCRx),mHPReg(:,iCRx));
    end
  end

  %% 预放大（可调节输出音量）
  mOut = 2 * mOut;  % 可改为 1.5 或 3 根据主观听感调节音量
  fMaxAmpl = max(max(abs(mOut(:))), fMaxAmpl);

  if bShowDisplay && mod(iCount,iNoIterShowDisplay)==1
    if fMaxAmpl>0.5
      disp(['Critical maximal amplitude: ',num2str(fMaxAmpl)]);
    end
  end

  %% 播放输出
  nUnderrun = play(deviceWriterActive, mOut);
  if nUnderrun > 0
    fprintf('Audio writer underrun by %d samples.\n', nUnderrun);
  end
end
