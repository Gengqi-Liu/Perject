Error using double
Conversion to double from cell is not possible.

Error in sgengqi (line 313)
    b = double(squeeze(mIIR_B(iCRx,iCTx,:)));
    
    
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
    if iCalCount > round(5*fSamplFreq/frameLength)
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
        mReg   = zeros(length(vIR)-1,iNoTx);
        mOut   = zeros(frameLength,iNoTx);
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
  mIn       = fileReader();

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
    iCMod                 = 1+mod(iCount-1,length(vAngleHorSave));
    vAngleHorSave(iCMod)  = fAngleHor; % save angle for debugging
   
    [~,iAngleIndOld]      = min(abs(vAngle-vAngleHorPred(1+mod(iCMod-2,length(vAngleHorSave)))));
    [~,iAngleIndCur]      = min(abs(vAngle-vAngleHorPred(iCMod)));

    %% Update ring buffer
    iCircPt_x             = mod(iCount-1,N_RB_x)+1;
    iCircPt_y             = mod(iCount-1,N_RB_y)+1;
    x_ring(:,iCircPt_x,:) = mIn;
    iCircPt_y_Del                 = mod(iCircPt_y-1-1,N_RB_y)+1;
    y_ring_old(:,iCircPt_y_Del,:) = 0;
    y_ring_cur(:,iCircPt_y_Del,:) = 0;

    %% SEGMENT 1
    if mod(iCount,Bi(1)/B==0) % when it is available, here every block
      iC1 = iC1 + 1;
      vInd1_x = iCircPt_x;
      mIn     = x_ring(:,vInd1_x,:);
      vIndUpdate = [iAngleIndOld,iAngleIndCur]; % will be called in any iteration
      % H1(:,:,:,:,vIndUpdate) = interpElevation(H1tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer); 
      [y_part1_old,x_in_buf1_old,mFDL_buf1_old] = UPConv(mIn,x_in_buf1_old,mFDL_buf1_old,vTxInd,H1,iC1,iAngleIndOld,Bi(1),Pi(1));
      [y_part1_cur,x_in_buf1_cur,mFDL_buf1_cur] = UPConv(mIn,x_in_buf1_cur,mFDL_buf1_cur,vTxInd,H1,iC1,iAngleIndCur,Bi(1),Pi(1));
      vInd1_y = iCircPt_y;
      y_ring_old(:,vInd1_y,:) = y_ring_old(:,vInd1_y,:) + reshape(y_part1_old,B,Bi(1)/B,[]);
      y_ring_cur(:,vInd1_y,:) = y_ring_cur(:,vInd1_y,:) + reshape(y_part1_cur,B,Bi(1)/B,[]);
    end
    %% SEGMENT 2
    if numel(Bi)>1 && mod(iCount,Bi(2)/B)==0 % when it is available
      iC2 = iC2 + 1;
      vInd2_x = mod(iCircPt_x+(-Bi(2)/B+1:0)-1,N_RB_x)+1;
      mIn     = reshape(x_ring(:,vInd2_x,:),Bi(2),[]);
      vIndUpdate = [iAngleIndOld,iAngleIndCur];
      % H2(:,:,:,:,vIndUpdate) = interpElevation(H2tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer); 
      [y_part2_old,x_in_buf2_old,mFDL_buf2_old] = UPConv(mIn,x_in_buf2_old,mFDL_buf2_old,vTxInd,H2,iC2,iAngleIndOld,Bi(2),Pi(2));
      [y_part2_cur,x_in_buf2_cur,mFDL_buf2_cur] = UPConv(mIn,x_in_buf2_cur,mFDL_buf2_cur,vTxInd,H2,iC2,iAngleIndCur,Bi(2),Pi(2));
      vInd2_y = mod(iCircPt_y+vBlockDelay(2)+(0:Bi(2)/B-1)-1,N_RB_y)+1;
      y_ring_old(:,vInd2_y,:) = y_ring_old(:,vInd2_y,:) + reshape(y_part2_old,B,Bi(2)/B,[]);
      y_ring_cur(:,vInd2_y,:) = y_ring_cur(:,vInd2_y,:) + reshape(y_part2_cur,B,Bi(2)/B,[]);
    end
    %% SEGMENT 3
    if numel(Bi)>2 && mod(iCount,Bi(3)/B)==0+vSchedOffset(3) % when it is available
      iC3 = iC3 + 1;
      vInd3_x = mod(iCircPt_x-vSchedOffset(3)+(-Bi(3)/B+1:0)-1,N_RB_x)+1;
      mIn     = reshape(x_ring(:,vInd3_x,:),Bi(3),[]);
      vIndUpdate = [iAngleIndOld,iAngleIndCur];
      % H3(:,:,:,:,vIndUpdate) = interpElevation(H3tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer); 
      [y_part3_old,x_in_buf3_old,mFDL_buf3_old]  = UPConv(mIn,x_in_buf3_old,mFDL_buf3_old,vTxInd,H3,iC3,iAngleIndOld,Bi(3),Pi(3));
      [y_part3_cur,x_in_buf3_cur,mFDL_buf3_cur]  = UPConv(mIn,x_in_buf3_cur,mFDL_buf3_cur,vTxInd,H3,iC3,iAngleIndCur,Bi(3),Pi(3));
      vInd3_y = mod(iCircPt_y-vSchedOffset(3)+vBlockDelay(3)+(0:Bi(3)/B-1)-1,N_RB_y)+1;
      y_ring_old(:,vInd3_y,:) = y_ring_old(:,vInd3_y,:) + reshape(y_part3_old,B,Bi(3)/B,[]);
      y_ring_cur(:,vInd3_y,:) = y_ring_cur(:,vInd3_y,:) + reshape(y_part3_cur,B,Bi(3)/B,[]);
    end
    %% SEGMENT 4
    if numel(Bi)>=4 && mod(iCount,Bi(4)/B)==0+vSchedOffset(4) % when it is available
      iC4 = iC4 + 1;
      vInd4_x = mod(iCircPt_x-vSchedOffset(4)+(-Bi(4)/B+1:0)-1,N_RB_x)+1;
      mIn     = reshape(x_ring(:,vInd4_x,:),Bi(4),[]);
      % H4(:,:,:,:,vIndUpdate) = interpElevation(H4tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer);       
      [y_part4_old,x_in_buf4_old,mFDL_buf4_old]  = UPConv(mIn,x_in_buf4_old,mFDL_buf4_old,vTxInd,H4,iC4,iAngleIndOld,Bi(4),Pi(4));
      [y_part4_cur,x_in_buf4_cur,mFDL_buf4_cur]  = UPConv(mIn,x_in_buf4_cur,mFDL_buf4_cur,vTxInd,H4,iC4,iAngleIndCur,Bi(4),Pi(4));
      vInd4_y = mod(iCircPt_y-vSchedOffset(4)+vBlockDelay(4)+(0:Bi(4)/B-1)-1,N_RB_y)+1;
      y_ring_old(:,vInd4_y,:) = y_ring_old(:,vInd4_y,:) + reshape(y_part4_old,B,Bi(4)/B,[]);
      y_ring_cur(:,vInd4_y,:) = y_ring_cur(:,vInd4_y,:) + reshape(y_part4_cur,B,Bi(4)/B,[]);
    end

    %% Take block from ring buffer
    mOut_Old  = squeeze(y_ring_old(:,iCircPt_y,:));
    mOut_Cur  = squeeze(y_ring_cur(:,iCircPt_y,:));
    % Fading: Combine old and current output
    mOut      = mWeightsDown.*mOut_Old + mWeightsUp.*mOut_Cur;


%% Apply IIR post-filtering (FIR → IIR hybrid stage)
% -------------------------------------------------------------------------
% Each receiver (2) × each transmitter (iNoTx=6)
mOut_IIR = zeros(size(mOut)); % same size as mOut (frameLength × iNoTx)
for iCRx = 1:2
  for iCTx = 1:iNoTx
    b = double(squeeze(mIIR_B(iCRx,iCTx,:)));
    a = double(squeeze(mIIR_A(iCRx,iCTx,:)));
    x = double(mOut(:,iCTx));
    [mOut_IIR(:,iCTx), mIIRReg(:,iCRx,iCTx)] = filter(b, a, x, mIIRReg(:,iCRx,iCTx));
  end
end

% Replace original mOut with filtered version
mOut = mOut_IIR;

    %% headphone equalization
    if 0%bHPEQ
      for iCRx=1:2
        [mOut(:,iCRx),mHPReg(:,iCRx)] = filter(vHPIR,1,mOut(:,iCRx),mHPReg(:,iCRx));
      end
    end

    % %% Nicos equalization
    % for iCRx=1:2
    %   for iCSec=1:size(mEQ_B,1)
    %     [mOut(:,iCRx),mEQReg(:,iCSec,iCRx)] = filter(mEQ_B(iCSec,:),mEQ_A(iCSec,:),mOut(:,iCRx),mEQReg(:,iCSec,iCRx));
    %   end
    % end
    % for iCRx=1:2
    %   [mOut(:,iCRx),mEQReg2(:,iCRx)] = filter(vEQ_B,vEQ_A,mOut(:,iCRx),mEQReg2(:,iCRx));
    % end


    % Signalverarbeitung (z. B. Block mit N Samples)
    for iCRx = 1:2
%         mOut(:,iCRx) = Fsos(msOut(:,iCRx));  % Filterung mit Zustand
    end

    %% Preamplifier
    mOut      = 2*mOut;
    fMaxAmpl  = max(max(abs(mOut(:))),fMaxAmpl);
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
  nUnderrun = play(deviceWriterActive,mOut);
  vUnderrun(mod(iCount-1,length(vRunTime))+1) = false;
  if nUnderrun > 0
    fprintf('Audio writer queue was underrun by %d samples.\n',...
      nUnderrun);
    vUnderrun(mod(iCount-1,length(vRunTime))+1) = true;
  end
end
