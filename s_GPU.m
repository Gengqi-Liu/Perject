addpath('utils');
if ~exist('mIRInt','var')
  load('utils/ActiveFilter.mat');  
  disp(sBRIRName);
  sDirName        = [sSetup,'/',sBRIRName];
  frameLength     = 2^8;
  iFiltLen        = 2^16;
  bInterpolate    = true;
  [mIRInt,vAngle] = PrepareRealTimeProc(sDirName,iFiltLen,bInterpolate);  
end
load Netflix_HD800
iNoTx   = size(mIRInt,3);
if iNoTx==2
  vTxInd = [1:2];
elseif iNoTx==6
%   vTxInd = [1:3,5:6];
  vTxInd = [1:6];
end
fid     = fopen('utils/AngleRange.txt','w');
fprintf(fid,'%f',max(abs(vAngle)));
fclose(fid);

%% Read file parameters
iNoIter       = 1;% frameLength;

%% Set audio interface
sAudioInterface = 'Babyface';
SetAudioInterface;

%% Real-time processing
% for speaker implementation
iDelay        = 50;
vIR           = [1,zeros(1,iDelay+1)];
mReg          = zeros(length(vIR)-1,iNoTx);

% for fading implementation
vWeightsUp    = [1:frameLength/iNoIter].'/(frameLength/iNoIter);
% vWeightsUp(1:frameLength/iNoIter) = sin([1:frameLength/iNoIter].'/(frameLength/iNoIter)*pi/2).^2;
mWeightsUp    = repmat(vWeightsUp,1,2);
mWeightsDown  = 1-mWeightsUp;

%% new
N = iFiltLen;
B = frameLength;
L = B;
P = N/L;
K = 2*B;

% buffer initialization
mFDL_buf_GPU  = gpuArray(single(zeros(B+1,P,length(vTxInd))));
x_in_buf_GPU  = gpuArray(single(zeros(2*B,length(vTxInd))));
Y_out_GPU     = gpuArray(single(zeros(2*B,2)));
vFDLInd       = 1:P;
iFiltLen      = 2^nextpow2(iFiltLen);
mTFInt        = single(fft(reshape(mIRInt,B,iFiltLen/B,2,iNoTx,[]),K,1));
mTFInt(B+2:end,:,:,:,:) = [];
mTFInt_GPU    = gpuArray(mTFInt);
clear mTFInt;

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
bHeadphone    = true;

% tic;
vAngleHorAll = [];
fMax          = 0;

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
  else
    %       disp('No value');
  end
  %     if iCount>1
  %       fAngleHor = fAngleHor+0.025*(fAngleHor-vAngleHorSave(iCount-1));
  %     end
  
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
%   mIn(:,2)  = 1.1*mIn(:,2); % balance
  
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
    %% Frequency-domain real-time convolution
    vAngleHorSave(iCount)     = fAngleHor;
    [~,iAngleIndOld]          = min(abs(vAngle-vAngleHorSave(max(iCount-1,1))));
    [~,iAngleIndCur]          = min(abs(vAngle-fAngleHor));
    
    %%
    
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
    Y_out_GPU(1:B+1,1)        = sum(sum(mFDL_buf_GPU(:,vFDLInd,:).*squeeze(mTFInt_GPU(:,:,1,vTxInd,iAngleIndOld)),3),2);
    Y_out_GPU(1:B+1,2)        = sum(sum(mFDL_buf_GPU(:,vFDLInd,:).*squeeze(mTFInt_GPU(:,:,2,vTxInd,iAngleIndOld)),3),2);
    y_out_tp_GPU              = ifft(Y_out_GPU,K,1,'symmetric');
    mOutOld_GPU               = y_out_tp_GPU(B+1:end,:);    
    
    % CUR angle: product of spectra, sum subspectra and IFFT
    Y_out_GPU(1:B+1,1)        = sum(sum(mFDL_buf_GPU(:,vFDLInd,:).*squeeze(mTFInt_GPU(:,:,1,vTxInd,iAngleIndCur)),3),2);
    Y_out_GPU(1:B+1,2)        = sum(sum(mFDL_buf_GPU(:,vFDLInd,:).*squeeze(mTFInt_GPU(:,:,2,vTxInd,iAngleIndCur)),3),2); 
    y_out_tp_GPU              = ifft(Y_out_GPU,K,1,'symmetric');
    mOutCur_GPU               = y_out_tp_GPU(B+1:end,:);
    
    % combine old and current contribution
    mOut_GPU    = mWeightsDown.*mOutOld_GPU + mWeightsUp.*mOutCur_GPU;
    if iNoTx>2
%       mOut_GPU = mOut_GPU + 0.5*mIn_GPU(:,4);
    end

    % preamplifier
    mOut      = 2*gather(mOut_GPU);
    
    fMax      = max(max(abs(mOut(:))),fMax);
    if mod(iCount,20)==1
      if fMax>0.5
        disp(['Critical maximal amplitude: ',num2str(fMax)]);
      end
    end
         
%     mOut(:,2) = 1.1*mOut(:,2);
%     mOut = 0.1*mIn;
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