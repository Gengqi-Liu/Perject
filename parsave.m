function parsave(temp_file, iCRx, iCTx, iAng, localB, localA)
% 并行安全地保存单个角度结果到主文件中
lockfile = [temp_file '.lock'];

% 简易锁机制，避免并行冲突
while isfile(lockfile)
    pause(0.05);
end
fclose(fopen(lockfile, 'w'));

try
    save_flag = false;
    if isfile(temp_file)
        load(temp_file, 'mIIR_B', 'mIIR_A', 'doneMask');
        save_flag = true;
    else
        save_flag = false;
    end

    % 如果第一次存储则初始化
    if ~save_flag
        mIIR_B = zeros(2,6,101,length(localB));
        mIIR_A = zeros(2,6,101,length(localA));
        doneMask = false(2,6,101);
    end

    % 写入结果
    mIIR_B(iCRx,iCTx,iAng,:) = localB;
    mIIR_A(iCRx,iCTx,iAng,:) = localA;
    doneMask(iCRx,iCTx,iAng) = true;

    % 临时保存
    save(temp_file,'mIIR_B','mIIR_A','doneMask','-v7.3');

    % 备份到主文件（原子替换）
    movefile(temp_file,'IIR_filters_allAngles.mat','f');

catch ME
    warning(['parsave() failed: ', ME.message]);
end

% 删除锁文件
if isfile(lockfile)
    delete(lockfile);
end
end
