addpath('utils');
SetParameters;
deviceWriterActive = deviceWriterAll;

%% Initialize NUPOLS
% sFilterName = 'AKG712_2.0_3D_Home'; % 1: -30° nach unten, 2: 0°, 3: 30° nach oben
% sFilterName = 'HD800_5.1_3D_Home'; % 1: -30° nach unten, 2: 0°, 3: 30° nach oben
% sFilterName = 'Netflix_HD800'; %_with_center
% sFilterName = 'Room_Home_221025_left_right_center_HP_HD800_221025';

% sFilterName = 'Room_Home_221025_5_1_HP_HD800_221025';
%
% load(['filters/',sFilterName,'.mat']);
%
% % temporarily
% mIRInt(:,1,3,:) = TimeShiftNotCyclic(mIRInt(:,1,3,:),-5);
% mIRInt(:,2,3,:) = TimeShiftNotCyclic(mIRInt(:,2,3,:),-5);

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
% filter_order = 12;
% crossover_frequency = 5000;

% [B_highpass, A_highpass] = butter( filter_order, crossover_frequency/fSamplFreq*2, 'high' );
% [B_lowpass,  A_lowpass ] = butter( filter_order, crossover_frequency/fSamplFreq*2, 'low'  );
%
% [output_low,  state_lowpass_1st ] = filter( B_lowpass,  A_lowpass,  input_pcm_samples, state_lowpass_1st );
% [output_low,  state_lowpass_2nd ] = filter( B_lowpass,  A_lowpass,  output_low,        state_lowpass_2nd );
% [output_high, state_highpass_1st] = filter( B_highpass, A_highpass, input_pcm_samples, state_highpass_1st);
% [output_high, state_highpass_2nd] = filter( B_highpass, A_highpass, output_high,       state_highpass_2nd);


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
vAngleHorSave = zeros(1,1e5,'single');
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
mOut          = zeros(frameLength,iNoTx);

%% Preparing IIR filter addition
iIIRDelay     = (size(H1,1)-1)*size(H1,2)+(size(H2,1)-1)*size(H2,2)+...
  (size(H3,1)-1)*size(H3,2)+(size(H4,1)-1)*size(H4,2); % 42240;
mInDelayed    = zeros(iIIRDelay+frameLength,iNoTx);
mIRIntDelayed = mIRInt(iIIRDelay+1:end,:,:,(size(mIRInt,4)+1)/2);
% mIRIntIIR_b   = mIRIntDelayed;
% mIRIntIIR_a   = ones(1,2,iNoTx);
iFiltOrd = round(size(mIRIntDelayed,1)/2)+1
mIRIntIIR_b = zeros(5001,2,iNoTx);
mIRIntIIR_a = zeros(101,2,iNoTx);
for iCTx=1:iNoTx
  for iCRx=1:2
    [mIRIntIIR_b(:,iCRx,iCTx),mIRIntIIR_a(:,iCRx,iCTx)] = prony(mIRIntDelayed(:,iCRx,iCTx),5000,100);
  end
end
mRegIIR       = zeros(size(mIRIntIIR_b,1)-1,2,iNoTx);
mOutIIRTp     = zeros(frameLength,2,iNoTx);
%% 


% For fading implementation: Output is calculated twice (old/current) and
% mixed
vWeightsUp    = (1:frameLength).'/frameLength;
mWeightsUp    = repmat(vWeightsUp,1,2);
mWeightsDown  = 1-mWeightsUp;


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

  %% Toggle between headphone and loudspeakers
  if fAngleVer>-55
    bHeadphone = true;
  else
    bHeadphone = false;
  end

  %% Read data
  mIn       = fileReader();

  %% static IIR filter addition
  vIndWrite = mod((iCount-1)*frameLength:iCount*frameLength-1,size(mInDelayed,1))+1;
  mInDelayed(vIndWrite,:) = mIn;
  vIndRead = mod(vIndWrite-iIIRDelay,size(mInDelayed,1))+1;
  mInIIR = mInDelayed(vIndRead,:);
  for iCRx=1:2
    for iCTx=1:iNoTx
      [mOutIIRTp(:,iCRx,iCTx),mRegIIR(:,iCRx,iCTx)] = ...
        filter(mIRIntIIR_b(:,iCRx,iCTx),mIRIntIIR_a(:,iCRx,iCTx),mInIIR(:,iCTx),mRegIIR(:,iCRx,iCTx));
    end
  end
  mOutIIR = sum(mOutIIRTp,3);

  %%

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
    %     [mIn,mIn2] = crossFilt(mIn);
    w = 1;
    % phi = atan(sqrt(2)/2); % Gerzon: 35.26°
    % phi = asin(2/3); % Pekonen: 41.8°
    % phi = atan(sqrt(2)); % Gerzon: 54.74°
    % phi = 90*pi/180; % inactive
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
  else
    %% Headphone is ON
    if bShowDisplay && mod(iCount,iNoIterShowDisplay)==1
      PrintStatus(sRoomName,sHeadphoneName,'ON',N,iNoTx,fAngleHor,fAngleVer,...
        vAngle,vRunTime,vUnderrun,fUpdateTime,iCount);
    end
  end
  %% speaker output is continuously calculated
  for iCTx=1:iNoTx
    [mOut(:,iCTx),mReg(:,iCTx)] = filter(vIR,1,mIn(:,iCTx),mReg(:,iCTx));
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
    vIndUpdate = [iAngleIndOld,iAngleIndCur];
    %       H1(:,:,:,:,vIndUpdate) = interpElevation(H1tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer);
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
    %       H2(:,:,:,:,vIndUpdate) = interpElevation(H2tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer);
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
    %       H3(:,:,:,:,vIndUpdate) = interpElevation(H3tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer);
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
    %       H4(:,:,:,:,vIndUpdate) = interpElevation(H4tp(:,:,:,:,vIndUpdate,:),vAngleVer,fAngleVer);
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
  mOut(:,iNoTx+1:iNoTx+2) = mWeightsDown.*mOut_Old + mWeightsUp.*mOut_Cur;

  %% IIR filter addition
  mOut(:,iNoTx+1:iNoTx+2) = mOutIIR + mOut(:,iNoTx+1:iNoTx+2);

  %% headphone equalization
  if bHPEQ
    for iCRx=1:2
      [mOut(:,iNoTx+iCRx),mHPReg(:,iCRx)] = filter(vHPIR,1,mOut(:,iNoTx+iCRx),mHPReg(:,iCRx));
    end
  end

  %% Preamplifier
  mOut(:,iNoTx+1:iNoTx+2) = 2*mOut(:,iNoTx+1:iNoTx+2);
  fMaxAmpl  = max(max(abs(mOut(:,iNoTx+1:iNoTx+2))),fMaxAmpl);
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

  %% Write data to output buffer
  if bHeadphone
    mOut(:,1:iNoTx) = 0;
  else
    mOut(:,iNoTx+1:iNoTx+2) = 0;
  end

  nUnderrun = play(deviceWriterActive,mOut);
  vUnderrun(mod(iCount-1,length(vRunTime))+1) = false;
  if nUnderrun > 0
    fprintf('Audio writer queue was underrun by %d samples.\n',...
      nUnderrun);
    vUnderrun(mod(iCount-1,length(vRunTime))+1) = true;
  end
end