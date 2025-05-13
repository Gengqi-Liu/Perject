addpath('utils');
%% First BRIR
load('utils/ActiveFilter.mat');
frameLength     = 2^9;
iFiltLen        = 2^16;
bInterpolate    = false;

%% 1st BRIR
sRoomName1        = convertStringsToChars(sRoomName);
sHeadphoneName1   = sHeadphoneName;
sDirName1         = [sSetup,'/Room_',sRoomName1,'_HP_',sHeadphoneName1];
[mIRInt1,vAngle]  = PrepareRealTimeProc(sDirName1,iFiltLen,bInterpolate);
fWeight1          = 0.5;
mIRInt_GPU_1      = gpuArray(single(fWeight1*mIRInt1));
clear mIRInt2;

%% 2nd BRIR
sRoomName2        = 'Home_3_Sub_200303';
sHeadphoneName2   = 'HD800_3_Sub_200303';
sDirName2         = [sSetup,'/Room_',sRoomName2,'_HP_',sHeadphoneName2];
[mIRInt2,~]       = PrepareRealTimeProc(sDirName2,iFiltLen,bInterpolate);
fWeight2          = 1;
mIRInt_GPU_2      = gpuArray(single(fWeight2*mIRInt2));
clear mIRInt2;

%% 3rd BRIR
sRoomName3        = 'Home_2_200303';
sHeadphoneName3   = 'HD800_2_200303';
sDirName3         = [sSetup,'/Room_',sRoomName3,'_HP_',sHeadphoneName3];
[mIRInt3,~]       = PrepareRealTimeProc(sDirName3,iFiltLen,bInterpolate);
fWeight3          = 1;
mIRInt_GPU_3      = gpuArray(single(fWeight3*mIRInt3));
clear mIRInt3;

%% 0th BRIR (pass through)
sRoomName0        = 'Passthrough';
sHeadphoneName0   = '-';
mIRInt_GPU_0      = gpuArray(single(zeros(size(mIRInt_GPU_1))));

%% At first, activate filter 1
fid               = fopen('FilterNumber.txt','w'); fwrite(fid,'1'); fclose(fid);
mIRInt_GPU        = mIRInt_GPU_1;
iNoTx             = size(mIRInt1,3);
vTxInd            = 1:iNoTx;
iNoAngles         = length(vAngle);
fAngleDist        = diff(vAngle(1:2));

%% Set audio interface
sAudioInterface   = 'Babyface';
SetAudioInterface;

%% Real-time processing
% for speaker implementation
iDelay            = 50;
vIR               = [1,zeros(1,iDelay+1)];
mReg              = zeros(length(vIR)-1,iNoTx);

% for fading implementation
vWeightsUp        = [1:frameLength].'/frameLength;
% vWeightsUp(1:frameLength) = sin([1:frameLength].'/(frameLength)*pi/2).^2;
mWeightsUp        = repmat(vWeightsUp,1,2);
mWeightsDown      = 1-mWeightsUp;

%% new
N = iFiltLen;
B = frameLength;
L = B;
P = N/L;
K = 2*B;

% buffer initialization
mFDL_buf_GPU    = gpuArray(single(zeros(B+1,P,length(vTxInd))));
x_in_buf_GPU    = gpuArray(single(zeros(2*B,length(vTxInd))));
Y_out_GPU       = gpuArray(single(zeros(2*B,2)));
vFDLInd         = 1:P;
iFiltLen        = 2^nextpow2(iFiltLen);
mTFInt          = single(fft(reshape(mIRInt_GPU,B,iFiltLen/B,2,iNoTx,[]),K,1));
mTFInt(B+2:end,:,:,:,:) = [];
mTFInt_GPU      = gpuArray(mTFInt);
clear mTFInt;
mTFInt_GPU_old  = gpuArray(single(zeros(B+1,iFiltLen/B,2,iNoTx)));
mTFInt_GPU_cur  = gpuArray(single(zeros(B+1,iFiltLen/B,2,iNoTx)));
mIRInt_GPU_cur  = gpuArray(single(zeros(iFiltLen,2,iNoTx)));

% debugging
vAngleHorSave   = zeros(1,5e5);
vAngleHorSerial = zeros(1,5e5);
vAngleHorNew    = zeros(1,5e5);
vAngleVerNew    = zeros(1,5e5);
vRunTime1       = zeros(1,5e5);
vRunTime2       = zeros(1,5e5);

disp('Real-time convolving starts')
iCount        = 0;
iCountNoVal   = 0;
fAngleHor     = 0;
fAngleVer     = 0;
bHeadphone    = true;

% tic;
vAngleHorAll = [];
fMax          = 0;

%% currently testing
fHeadRadius     = 0.09;%0.08;
fDist2Source    = 2;
vAngleSource    = [30,-30,0,0,105,-105];
load('utils/mITD_table.mat','mITD_table');
%%
while iCount<10000000

  iCount = iCount + 1;
  
  %% Read angles from file
  fid       = fopen('utils/Angles.txt', 'r');
  % Matlab headtracker
  vAngleNew = str2num(fscanf(fid,'%c',[8,2]).'); % %s für Matlab headtracker??
  fclose(fid);  
  % update angles if possible
  if numel(vAngleNew)>0
    fAngleHor = vAngleNew(1);
    if numel(vAngleNew)>1
      fAngleVer = vAngleNew(2);
    end
  end
  
  %% Toggle between headphone and loudspeakers
  if mod(iCount,5)==1
    if bHeadphone
      if fAngleVer<-50
        bHeadphone    = false;
        release(deviceWriterActive);
        deviceWriterActive = deviceWriterSpeaker;
        mReg   = zeros(length(vIR)-1,iNoTx);
        mOut   = zeros(frameLength,iNoTx);
      end
    else % Speakers are on
      if fAngleVer>-55
        bHeadphone    = true;
        release(deviceWriterActive);
        deviceWriterActive = deviceWriterHeadphone;     
        mFDL_buf_GPU  = gpuArray(single(zeros(B+1,P,length(vTxInd))));
        x_in_buf_GPU  = gpuArray(single(zeros(2*B,length(vTxInd))));
        Y_out_GPU     = gpuArray(single(zeros(2*B,2)));
        vFDLInd       = 1:P;
        iCount        = 1;
      end
    end
  end
  
  %% Read data
  mIn       = fileReader();
  mIn_GPU   = gpuArray(mIn);
  
  %% Check whether headphone is on
  if ~bHeadphone 
    %% Headphone is off
    if mod(iCount,20)==1
      clc;
      fprintf('Binaural Synthesizer\n');
      fprintf('====================\n\n');
      fprintf('Room:       %s\n',sRoomName);
      fprintf('Headphone:  %s\n\n',sHeadphoneName);
      fprintf('Headphone is OFF!\n');
    end
    %% Time-delay of input signal
    for iCTx=1:iNoTx
      [mOut(:,iCTx),mReg(:,iCTx)] = filter(vIR,1,mIn(:,iCTx),mReg(:,iCTx));
    end
  else    
    %% Headphone is on
    if mod(iCount,20)==1
      clc;
      fprintf('Binaural Synthesizer\n');
      fprintf('====================\n\n');
      fprintf('Room:       %s\n',sRoomName);
      fprintf('Headphone:  %s\n\n',sHeadphoneName);
      fprintf('Headphone is ON!\n\n');
      if abs(fAngleHor) > max(abs(vAngle))
        fprintf('Horizontal orientation: %6.1f° (Out of range!)\n',fAngleHor);
      else
        fprintf('Horizontal orientation: %6.1f°\n',fAngleHor);
      end
      fprintf('Vertical orientation:   %6.1f°\n\n',fAngleVer);
    end
    %% 
    if mod(iCount,20)==1
      if exist('FilterNumber.txt','file')
        fid             = fopen('FilterNumber.txt', 'r');
        sNumber         = fscanf(fid,'%c',[1]);
        mIRInt_GPU      = eval(['mIRInt_GPU_',sNumber]);
        sRoomName       = eval(['sRoomName',sNumber]);
        sHeadphoneName  = eval(['sHeadphoneName',sNumber]);
        fclose(fid);
      end
    end
    %% Update filter
    % old filter
    mTFInt_GPU_old          = mTFInt_GPU_cur;
    % calculate current filter
    % 1st step: linear interpolation
    vAngleHorSave(iCount)   = fAngleHor;
    vAngleDiff              = fAngleHor-vAngle;
    [fAngleMin,iAngleInd]   = min(vAngleDiff(vAngleDiff>=0));
    if numel(iAngleInd)
      fLinB                 = fAngleMin/fAngleDist;
      fLinA                 = 1-fLinB;
      mIRInt_GPU_cur        = fLinA*mIRInt_GPU(:,:,:,iAngleInd)+...
        fLinB*mIRInt_GPU(:,:,:,min(iNoAngles,iAngleInd+1));
    else
      mIRInt_GPU_cur        = mIRInt_GPU(:,:,:,1);
    end
    % 2nd step: time shift
    for iCTx=vTxInd
      [fITD_calc,mTD_calc]  = CalcITD(fHeadRadius,fDist2Source,vAngleSource(iCTx),fAngleHor);
      fITD_calc             = interp1(mITD_table(1,:),mITD_table(2,:),vAngleSource(iCTx)-fAngleHor,'pchip');
      vTD_calc              = mTD_calc(1)+[0,fITD_calc];
      for iCRx=1:2
        fTimeDelay          = vTD_calc(iCRx);
        iTimeDelay          = -round(fTimeDelay);
        if iTimeDelay>0
          vIndTimeDel                 = [1+iTimeDelay:iFiltLen,1:iTimeDelay];
          mIRInt_GPU_cur(:,iCRx,iCTx) = mIRInt_GPU_cur(vIndTimeDel,iCRx,iCTx);
        elseif iTimeDelay<0
          vIndTimeDel                 = [iFiltLen+iTimeDelay:iFiltLen,1:iFiltLen+iTimeDelay-1];
          mIRInt_GPU_cur(:,iCRx,iCTx) = mIRInt_GPU_cur(vIndTimeDel,iCRx,iCTx);
        end
      end
    end
    mTFInt_GPU_cur                = single(fft(reshape(mIRInt_GPU_cur,B,iFiltLen/B,2,iNoTx),K,1));
    mTFInt_GPU_cur(B+2:end,:,:,:) = [];    

    %% Frequency-domain real-time convolution    
    % fill input buffer    
    x_in_buf_GPU(1:B,:,:)     = x_in_buf_GPU(B+1:end,:,:);    
    x_in_buf_GPU(B+1:end,:,:) = mIn_GPU(:,vTxInd);
    
    % FFT
    mInFFT_GPU                = fft(x_in_buf_GPU,K,1);

    % FDL
    iCircPt                   = mod(P-iCount,P)+1;
    mFDL_buf_GPU(:,iCircPt,:) = mInFFT_GPU(1:B+1,:);
    vFDLInd                   = circshift(vFDLInd,1);

    % OLD angle: product of spectra, sum subspectra and IFFT
    Y_out_GPU(1:B+1,1)        = sum(sum(mFDL_buf_GPU(:,vFDLInd,:).*squeeze(mTFInt_GPU_old(:,:,1,vTxInd)),3),2);
    Y_out_GPU(1:B+1,2)        = sum(sum(mFDL_buf_GPU(:,vFDLInd,:).*squeeze(mTFInt_GPU_old(:,:,2,vTxInd)),3),2);
    y_out_tp_GPU              = ifft(Y_out_GPU,K,1,'symmetric');
    mOutOld_GPU               = y_out_tp_GPU(B+1:end,:);    
    
    % CUR angle: product of spectra, sum subspectra and IFFT
    Y_out_GPU(1:B+1,1)        = sum(sum(mFDL_buf_GPU(:,vFDLInd,:).*squeeze(mTFInt_GPU_cur(:,:,1,vTxInd)),3),2);
    Y_out_GPU(1:B+1,2)        = sum(sum(mFDL_buf_GPU(:,vFDLInd,:).*squeeze(mTFInt_GPU_cur(:,:,2,vTxInd)),3),2); 
    y_out_tp_GPU              = ifft(Y_out_GPU,K,1,'symmetric');
    mOutCur_GPU               = y_out_tp_GPU(B+1:end,:);
    
    % combine old and current contribution
    mOut_GPU                  = mWeightsDown.*mOutOld_GPU + mWeightsUp.*mOutCur_GPU;
    if iNoTx>2
%       mOut_GPU = mOut_GPU + 0.5*mIn_GPU(:,4);
    end

    % preamplifier
    mOut                      = 2*gather(mOut_GPU);
    
    fMax                      = max(max(abs(mOut(:))),fMax);
    if mod(iCount,20)==1
      if fMax>0.5
        disp(['Critical maximal amplitude: ',num2str(fMax)]);
      end
    end
         
    if ~str2double(sNumber)
      mOut = 0.35*gather(single(mIn_GPU));
    end
%     vRunTime2(iCount) = toc;
    
    %% Save for debugging
%     mOutSave(:,iCount,:) = mOut;
  end
  %% 
  nUnderrun = play(deviceWriterActive,mOut);
  if nUnderrun > 0
    fprintf('Audio writer queue was underrun by %d samples.\n',...
      nUnderrun);
  end
end
release(fileReader);
release(deviceWriter);