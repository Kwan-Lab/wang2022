clear all

overAllPath = 'F:\pupildata';
dataPath = 'F:\For Hongli\M3_phase2_algo2_WithPupil_1911011021';
timePath = 'F:\pupildata\TestTime';
filename = 'M3_phase2_algo2_WithPupil_1911011021';
% option = '2';
% if option is 1, concatenate video as different trials
% else, is option is 2, concatenate video as a whole session
cd(overAllPath);
if ~exist(timePath,'dir')
    mkdir(timePath);
end

cd(dataPath);

videoFile= dir('*.avi');
matFile = dir('*.mat');
logFile = dir('*.log');

delayTime = load(matFile.name);

% for ii = 3:length(dataDir)
    
    timeFilePath = [timePath, '\', filename];
    disp('working on...');
    disp(timeFilePath);
    if ~exist(timeFilePath,'dir')
        mkdir(timeFilePath);
    end
    fullpath = [dataPath];
    cd(fullpath);
%     matFiles = dir('*.mat');
%     logFile = dir('*.log');
    
    if isempty(logFile)
        display(fullpath)
    else
    
%         [logData ] = MP_parseLogfile( fullpath, logFile.name );
%     
%         [ sessionData, trialData ] = MP_getSessionData( logData );
%     
%         [ trials ] = MP_getTrialMasks( trialData );
        phase = 8;
        [ logData ] = MP_parseLogfileMixStructure (dataPath, logFile.name);
        
        [ sessionData, trialData ] = MP_getSessionData_bandit( logData,phase );
        % load the video
        
        v = VideoReader(videoFile.name);
        numFrames = v.Duration*v.FrameRate;
        
        for ii = 1:numFrames
            frames{ii} = read(v, ii);
            if mod(ii,1000)==0
                ii
            end
        end
        % frames = read(v, [1, inf]);
        % [x,y,~,lenFrames] = size(frames);
        %grayFrames = zeros(x,y,lenFrames);
        lenFrames = numFrames;
%         for ii =1:numFrames
%             grayFrames{ii} = rgb2gray(frames{ii});
%         end
        grayFrames = frames;
        clear frames
        
        % calculate the frames
        videoFrames = length(grayFrames);
        
        logFrameTrial = zeros(1, length(trialData.trigger));
        
        for ii = 1:length(trialData.trigger)
            if ii < length(trialData.trigger)
                trialTime = trialData.triggerTimes(ii+1) - trialData.triggerTimes(ii);
                logFrameTrial(ii) = trialTime * 20;
            else
                trialTime = 0; % skip it for now
            end
        end
        
        % get the trial frames from matlab
        matFrameTrial = delayTime.framePerTrial(delayTime.framePerTrial~=0);
        % ignore 2 for now
        matFrameTrial_no2 = matFrameTrial(matFrameTrial ~= 2 & matFrameTrial ~= 1 & matFrameTrial ~= 3);
        
        figure; plot(matFrameTrial_no2);
        hold on; plot(logFrameTrial);  % looks fine
        
        
%         if length(matFiles) > length(trialData.triggerTimes)+1
%             display("Warning: video recordings are more than trials")
%             display(length(matFiles)-length(trialData.triggerTimes))
%             iter = length(trialData.triggerTimes)+1;
%         else
%             iter = length(matFiles);
%         end
%         
%         triggerDelayTime = [];
        % saveDelayTime = [];
        IRMask = zeros(1, lenFrames);
        IRSum = zeros(1, lenFrames);
        IRPositiveFrame = cell(1);
        IRPositiveInd = [];
        % get the TriggerDelayTime and the number of frames first
        cc = 1;
        
        for kk = 10000:500:20000
            figure;imagesc(grayFrames{kk});
            % #8 and 17 are IR positive
        end
        for ii =1:lenFrames
            IRSum(ii) = sum(sum(grayFrames{ii}));
             
                
            if IRSum(ii) > 5000000
                IRMask(ii) = 1;
                
                IRPositiveFrame{cc} = grayFrames{ii};
                IRPositiveInd = [IRPositiveInd, ii];
                        if IRSum(ii) < 6000000
                            figure;imagesc(grayFrames{ii})
                        end
                        cc = cc + 1;
                    
            end
        end
  
    end
    
    % don't mark the IR positive frames yet
    
    % calculate time based on the delayTime
    
    % get the first frame
%     startFrame = ceil(delayTime.delayTime/0.05) + 1;
%     frameTime = zeros(1, lenFrames);
%     frameTime(startFrame) = (startFrame-1) * 0.05 -  delayTime.delayTime +trialData.triggerTimes; % time 0 is the same as the trigger time in logfile
%     for kk = 1:lenFrames
%         frameTime(kk) = frameTime(startFrame) - 0.05*(startFrame - kk);
%     end
    frameTime = zeros(1, lenFrames);
    for kk = 1:lenFrames
        frameTime(kk) = trialData.triggerTimes(1)-delayTime.delayTime + 0.05*(kk-1);
    end
    % now align to the logfile time
    figure; plot(frameTime)
    hold on; plot(trialData.cueTimes);
    % align the times to cue time = 0
    alignedTime = zeros(length(trialData.cueTimes), 5);
    for jj = 1:length(trialData.cueTimes)
        alignedTime(jj, 1) = trialData.IRTimes(jj*3 - 2) - trialData.cueTimes(jj); % IR1 time
        alignedTime(jj,2) = 0; % cuetimes
        alignedTime(jj, 3) = trialData.IRTimes(jj*3 -1) - trialData.cueTimes(jj); % IR2 time
        alignedTime(jj, 4) = trialData.outcomeTimes(jj) - trialData.cueTimes(jj); % outcome time
        alignedTime(jj, 5) = trialData.IRTimes(jj * 3) - trialData.cueTimes(jj); % Ir3 time
    end
    
    % aligned the time to the IRFlash in the video
    % based on the IR time, find the video frames
    IRPosFrames = [];
    for ii = 1:lenFrames
        for kk = 1:length(trialData.IRTimes)
            if mod(kk,3) == 1
                if frameTime(ii) >= trialData.IRTimes(kk) & frameTime(ii) <= trialData.IRTimes(kk) + 0.8
                    IRPosFrames = [IRPosFrames, ii];
                end
            else
                if frameTime(ii) >= trialData.IRTimes(kk) & frameTime(ii) <= trialData.IRTimes(kk) + 0.1
                    IRPosFrames = [IRPosFrames, ii];
                end
            end
        end
    end
    
    % check the frames every trial
    timeTrial = diff(trialData.triggerTimes);
    frameTrial = floor(timeTrial * 20);
    frameRecord = [delayTime.framePerTrial(1:164),delayTime.framePerTrial(166:500)];
    % check whether are they actually positive
    figure;plot(IRSum(IRPosFrames));
    figure; scatter(1:length(IRPosFrames), IRPosFrames);
    hold on; scatter(1:length(IRPositiveInd), IRPositiveInd);
    
    
    figure;plot(frameTime(1:75000),IRSum(1:75000));
    hold on;
    
    startTime = frameTime(2000); endTime = frameTime(3000);
    IRPosTimes = frameTime(IRPosFrames);
    IRMask = zeros(1, lenFrames);
    IRMask(IRPosFrames) = 10000000;
    plot(frameTime(2000:3000), IRMask(2000:3000));
    
    figure;plot(frameTime(1000:2000),IRSum(1000:2000));
    hold on;
    plot(frameTime(1000:2000), IRMask(1000:2000));
    
    figure;plot(frameTime(5000:10000),IRSum(5000:10000));
    hold on;
    plot(frameTime(5000:10000), IRMask(5000:10000));
    
    figure;plot(frameTime(1:end),IRSum(1:end));
    hold on;
    plot(frameTime(1:end), IRMask(1:end));
    
    
    figure;plot(frameTime(55000:56000),IRSum(55000:56000));
    hold on;
    plot(frameTime(55000:56000), IRMask(55000:56000));
    
    figure;plot(frameTime(75000:76000),IRSum(75000:76000));
    hold on;
    plot(frameTime(75000:76000), IRMask(75000:76000));
    
     figure;plot(frameTime(end-1000:end),IRSum(end-1000:end));
    hold on;
    plot(frameTime(end-1000:end), IRMask(end-1000:end));
%         fileTrigger = [timeFilePath,'\',matFiles(1).name(1:end-8),'_triggerDelayTime.mat'];
%         fileSave = [timeFilePath,'\', matFiles(1).name(1:end-8), '_saveDelayTime.mat'];
%         fileFrames = [timeFilePath,'\',matFiles(1).name(1:end-8),'_Frames.mat'];
%         save(fileTrigger, 'triggerDelayTime');
%         save(fileSave, 'saveDelayTime');
%         save(fileFrames, 'numFrames');
%         
%         calculate the time based on the triggerDelayTime and number of
%         Frames
%             get logfile trials information
%         sessionTime = [];
%         trialMasks = [];
%         for tt = 2:iter
%             timeMat = zeros(1, numFrames(tt-1));
%             trialMask = ones(1, numFrames(tt-1)) * (tt-1);
%             for jj = 1:numFrames(tt-1)
%                 for new version
%                 timeMat(jj) = trialData.triggerTimes(tt-1)+triggerDelayTime(tt-1) + 0.1 * (jj - 1); 
%                 
%                 for old version
%                 timeMat(jj) = trialData.triggerTimes(tt-1)+triggerDelayTime(tt-1) + 0.1 * (jj - 1); 
%                 I=mat2gray(data{jj});
%                 J = filter2(fspecial('sobel'),I);
%                 K = mat2gray(J);
%                 imshow(K);
%             end
%             sessionTime = [sessionTime, timeMat];
%             trialMasks = [trialMasks, trialMask];
%         end
        
        
        %% run some plot
        
        % get the time difference from the last frame to the next first
        % frame
        gapTime = zeros(1, length(trialData.cue)-1);
        trialNum = 1;
        for vv = 2:length(trialMasks)
            if trialMasks(vv) ~= trialMasks(vv-1)
                gapTime(trialNum) = sessionTime(vv) - sessionTime(vv-1);
                trialNum = trialNum + 1;
            end
        end
        figure;plot(gapTime)
        figure;scatter(gapTime, numFrames(1:49));
        % check the IRtimes and cue times
        TimeDifference = zeros(7, length(trialData.triggerTimes));
        for uu = 1:length(trialData.triggerTimes)
            TimeDifference(1, uu) = trialData.IRTimes(uu * 3 -2) - trialData.triggerTimes(uu);
            % TimeDifference(2, uu) = trialData.IRTimes(uu * 2) - trialData.triggerTimes(uu);
            % TimeDifference(3, uu) = trialData.IRTimes(uu * 2) - trialData.IRTimes(uu * 2 - 1);
            TimeDifference(2, uu) = trialData.cueTimes(uu) - trialData.IRTimes(uu * 3 - 2);
            TimeDifference(3, uu) = trialData.IRTimes(uu * 3 - 1) - trialData.cueTimes(uu);
            TimeDifference(4, uu) = trialData.outcomeTimes(uu) - trialData.IRTimes(uu * 3 - 1);
            TimeDifference(5, uu) = trialData.saveTimes(uu) - trialData.outcomeTimes(uu);
            TimeDifference(6, uu) = trialData.IRTimes(uu * 3) - trialData.saveTimes(uu);
        end
        figure;
        subplot(3, 2, 1);
        plot(TimeDifference(1, :));
        title('IR1-trigger')
        subplot(3, 2, 2);
        plot(TimeDifference(2, :));
        title('cue-IR1');
        subplot(3, 2, 3);
        plot(TimeDifference(3, :));
        title('IR2-cue');
        subplot(3, 2, 4);
        plot(TimeDifference(4, :));
        title('outcome-IR2');
        subplot(3, 2, 5);
        plot(TimeDifference(5, :));
        title('save-outcome');
        subplot(3, 2, 6);
        plot(TimeDifference(6, :));
        title('IR3-save');
        
%         figure;
%         plot(trialData.IRTimes(1:3:14));
%         hold on; plot(trialData.cueTimes(1:5));
        % find the IR Time and frames
        IRFrameTimes = sessionTime(IRMask==1);
        %mark whether a frame is missed or not
        IRStart = zeros(1, length(trialData.cue));
        IRNL = zeros(1, length(trialData.cue));
        repeated = 0;
        diffStart = [];
        diffNL = [];
        IRFrameTime = [];
        for hh = 1:length(trialData.trigger)
            lowerBound1 = trialData.IRTimes(3 * hh - 2)-0.8;
            upperBound1 = trialData.IRTimes(3 * hh - 2) - 0.5;  % a 50 ms +/- for fluctuation
            lowerBound2 = trialData.IRTimes(3 * hh - 1)- 0.8;
            upperBound2 = trialData.IRTimes(3 * hh - 1) - 0.5;
            lowerBound3 = trialData.IRTimes(3 * hh)- 0.8;
            upperBound3 = trialData.IRTimes(3 * hh) - 0.5;
            [x1,y1] = find(IRFrameTimes < upperBound1 & IRFrameTimes > lowerBound1);
            [x2,y2] = find(IRFrameTimes < upperBound2 & IRFrameTimes > lowerBound2);
            if x1
                IRStart(hh) = 1;
                diffStart = [diffStart, trialData.IRTimes(2*hh -1) - IRFrameTimes(y1(1))];
                IRFrameTime = [IRFrameTime, IRFrameTimes(y1)];
            end
            if x2
                IRNL(hh) = 1;
                diffNL = [diffNL, trialData.IRTimes(2*hh) - IRFrameTimes(y2(1))];
            end
            
            if length(x2) > 1 
                repeated = repeated + 1;
            end
            if length(x1) > 1
                repeated = repeated + 1;
            end
        end
        
        figure;
        subplot(2,1,1)
        plot(diffStart);
        subplot(2,1,2)
        plot(diffNL)
 
        %%
        filename = [timeFilePath,'\',matFiles(1).name(1:end-8),'_Time.mat'];

        save(filename, 'sessionTime');
 
