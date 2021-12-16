% analyze pupil raw data using code adopted from
% Preprocessing Pupil Size Data. Guideline and Code.   Mariska Kret & Elio Sjak-Shie. 2018.
setup_figprop;

rawDataDir = 'E:\data\matching_pennies\874pupil';
rawDataPath = 'E:\data\matching_pennies\874pupil\M874_phase2_algo2_WithPupil_1910221112DeepCut_resnet50_PupillometryNov4shuffle1_150000.csv';
savebehfigpath = 'E:\data\matching_pennies\874pupil\beh'
rawDataPath2 = 'E:\data\matching_pennies\874pupil\M874_phase2_algo2_WithPupil_1910221112DeepCut_resnet50_PupillometryNov4shuffle1_150000_iter1.csv';

%rawTimePath1 = 'F:\pupildata\MP_pupilTime_1.1\M862_phase2_algo2_WithPupil_1905140746\M862_phase2_algo2_WithPupil_1905140746_Time.mat';
%rawTimePath2 = 'F:\pupildata\MP_pupilTime\M862_phase2_algo2_WithPupil_1905140746\M862_phase2_algo2_WithPupil_1905140746_Time.mat';

%load csv files
M = csvread(rawDataPath,3);
upX = M(:, 2); upY = M(:,3);
downX = M(:,5); downY = M(:,6);
leftX = M(:,8); leftY = M(:,9);
rightX = M(:,11); rightY = M(:,12);
centerX = M(:,14); centerY = M(:,15);

diaVer = sqrt((upX-downX).^2+(upY-downY).^2);
diaHor = sqrt((leftX-rightX).^2+(leftY-rightY).^2);

M2 = csvread(rawDataPath2,3);
upX2 = M2(:, 2); upY2 = M2(:,3);
downX2 = M2(:,5); downY2 = M2(:,6);
leftX2 = M2(:,8); leftY2 = M2(:,9);
rightX2 = M2(:,11); rightY2 = M2(:,12);
centerX2 = M(:,14); centerY2 = M2(:,15);

diaVer2 = sqrt((upX2-downX2).^2+(upY2-downY2).^2);
diaHor2 = sqrt((leftX2-rightX2).^2+(leftY2-rightY2).^2);

%sessionTime1= load(rawTimePath1);
%sessionTime2 = load(rawTimePath2);
figure;
plot(diaVer);
hold on; plot(diaHor);

figure;
plot(diaHor);hold on; plot(diaHor2)


%% align the time
cd(rawDataDir)
matFile = dir('*.mat');
logFile = dir('*.log');

delayTime = load(matFile(1).name);


if isempty(logFile)
    display('Logfile not found')
else
    
    [logData ] = MP_parseLogfile( rawDataDir, logFile.name );
    
    [ sessionData, trialData ] = MP_getSessionData( logData );
    
    [ trials ] = MP_getTrialMasks( trialData );
    
        %% setting up data for plotting
    
    stats.playerlabel{1} = 'Mouse';
    stats.playerlabel{2} = 'Algo2';
    
    %choice: left=-1; right=1; miss=NaN
    stats.c=nan(sessionData.nTrials,2);
    stats.c(trials.left,1)=-1;
    stats.c(trials.right,1)=1;
    stats.c(trials.compleft,2)=-1;
    stats.c(trials.compright,2)=1;
 
    %outcome: reward=1; no reward:0; miss=NaN
    stats.r=nan(sessionData.nTrials,1);
    stats.r(trials.reward)=1;
    stats.r(trials.noreward)=0;
    
    %% analysis of behavioral performance
    
    % plot behavior in raw format
    tlabel=strcat('Subject=',char(logData.subject),', Time=',char(logData.dateTime(1)),'-',char(logData.dateTime(2)));
    time_range=[-2 6];
%    plot_session_beh_vert(trialData,trials,[],tlabel,time_range);

    % plot choice behavior - whole sessions
    cd(savebehfigpath);
    tlabel=strcat('Subject=',logData.subject,', Time=',logData.dateTime(1),'-',logData.dateTime(2));
    
    plot_session_game(stats,sessionData.nTrials,tlabel);
    
    % Number of video frames is equal to the csv pupil coordinate length
    lenFrames = length(diaHor);

        
        % calculate the frames
        
    logFrameTrial = zeros(1, length(trialData.cue));
        
    for ii = 1:length(trialData.cue)
        if ii < length(trialData.cue)
            trialTime = trialData.triggerTimes(ii+1) - trialData.triggerTimes(ii);
            logFrameTrial(ii) = trialTime * 20;
        else
            trialTime = 0; % skip it for now
        end
     end
        
        % get the trial frames from matlab
   matFrameTrial = delayTime.framePerTrial(delayTime.framePerTrial~=0);
        % ignore 2 for now
   matFrameTrial_no2 = matFrameTrial(matFrameTrial ~= 2);
        
   figure; plot(matFrameTrial_no2);
   hold on; plot(logFrameTrial);  % looks fine
end

% get the frame time
frameTime = zeros(1, lenFrames);
for kk = 1:lenFrames
    frameTime(kk) = trialData.triggerTimes(1)-delayTime.delayTime + 0.05*(kk-1);
end

%% plot some frames
figure;
xAxis = 1:500;
figure; 
subplot(2,2,1);
scatter(xAxis, diaHor(xAxis),'.');
subplot(2,2,2);
scatter(xAxis+10000, diaHor(xAxis+10000), '.');
subplot(2,2,3)
scatter(xAxis+20000, diaHor(xAxis+20000),'.');
subplot(2,2,4);
scatter(xAxis+30000, diaHor(xAxis+30000),'.');

%% analyze

% set the blink to NaN, looks like the normal pupil diameter won't exceed
% 70 (horizontal)

dia_noBlink = diaHor;
dia_noBlink(dia_noBlink > 70 | dia_noBlink < 25) = NaN;

% get z-score
dia_z = (dia_noBlink - nanmean(dia_noBlink)) / nanstd(dia_noBlink);
figure; plot(dia_z)

% use running average to get z-score / 10 mins: 
windowSize = 10*60*20;
dia_runningz = zeros(1, length(dia_z));
for jj = 1:windowSize:length(dia_noBlink)
    if jj+windowSize-1<length(dia_z)
        dia_runningz(jj:jj+windowSize-1) = (dia_noBlink(jj:jj+windowSize-1) - nanmean(dia_noBlink(jj:jj+windowSize-1))) / nanstd(dia_noBlink(jj:jj+windowSize-1));
    else
        dia_runningz(jj:length(dia_z)) = (dia_noBlink(jj:length(dia_z)) - nanmean(dia_noBlink(jj:length(dia_z)))) / nanstd(dia_noBlink(jj:length(dia_z)));
    end
end

figure; plot(dia_runningz)
% overall average across session, aligned to cue
xSession = -1:0.05:5.95;
pupilDia_session = zeros(length(trialData.cueTimes),140);
for kk = 1:length(trialData.cueTimes)
    timeLower = trialData.cueTimes(kk)-1;
    timeHigher = trialData.cueTimes(kk) + 6;
    frameIndex = frameTime >= timeLower & frameTime <= timeHigher;
    pupilDia_session(kk,:) = dia_z(frameIndex);
end
mean_session = nanmean(pupilDia_session);
figure; plot(xSession, mean_session);

xAxis = -3:0.05:2.95;
% align to the cue time
pupilDia_cue = zeros(length(trialData.cueTimes), 120);

for kk = 1:length(trialData.cueTimes)
    timeLower = trialData.cueTimes(kk)-3;
    timeHigher = trialData.cueTimes(kk) + 3;
    frameIndex = frameTime >= timeLower & frameTime <= timeHigher;
    pupilDia_cue(kk,:) = dia_z(frameIndex);
end

mean_cue = nanmean(pupilDia_cue);
figure;plot(xAxis,mean_cue);
%s = shadedErrorBar(xAxis,pupilDia_cue,{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
%set(s.edge,'LineWidth',1,'LineStyle',':')
title('Pupil diameter aligned to cue time')

% align to the response
pupilDia_rt = zeros(length(trialData.cueTimes), 120);

for kk = 1:length(trialData.cueTimes)
    if isnan(trialData.rt(kk))
        pupilDia_rt(kk,:) = NaN;
    else
        timeLower = trialData.rt(kk)-3;
        timeHigher = trialData.rt(kk) + 3;
        frameIndex = frameTime >= timeLower & frameTime <= timeHigher;
        pupilDia_rt(kk,:) = dia_z(frameIndex);
    end
end

mean_rt = nanmean(pupilDia_rt);
figure;plot(xAxis, mean_rt);
title('Pupil diameter aligned to response time')
% align to the reward

pupilDia_reward = zeros(length(trialData.cueTimes),120);

for kk = 1:length(trialData.cueTimes)
    if isnan(trialData.rt(kk))
        pupilDia_reward(kk,:) = NaN;
    else
        timeLower = trialData.outcomeTimes(kk)-3;
        timeHigher = trialData.outcomeTimes(kk) + 3;
        frameIndex = frameTime >= timeLower & frameTime <= timeHigher;
        pupilDia_reward(kk,:) = dia_z(frameIndex);
    end
end
mean_reward = nanmean(pupilDia_reward(trials.reward,:));
mean_noreward = nanmean(pupilDia_reward(trials.noreward,:));
mean_miss = nanmean(pupilDia_reward(trials.miss,:));

figure;plot(xAxis,mean_reward);
hold on; plot(xAxis,mean_noreward);
legend('Reward','No reward');
title('Pupil diameter aligned to outcome time')

% separate into left/right choice/reward

% left/right choice
mean_left = nanmean(pupilDia_rt(trials.left,:));
mean_right = nanmean(pupilDia_rt(trials.right,:));
figure;plot(xAxis,mean_left);
hold on; plot(xAxis,mean_right);
legend('Left choice', 'Right choice');
title('Pupil diameter aligned to response time')

% left/right reward
mean_rewardleft = nanmean(pupilDia_reward(trials.reward & trials.left,:));
mean_rewardright = nanmean(pupilDia_reward(trials.reward & trials.right,:));
figure;plot(xAxis,mean_rewardleft);
hold on; plot(xAxis,mean_rewardright);
legend('Left reward', 'Right reward');
title('Pupil diameter aligned to outcome time')

% missed trial, aligned to the cue time
mean_nomiss = nanmean(pupilDia_cue(~trials.miss,:));
mean_miss = nanmean(pupilDia_cue(trials.miss,:));
figure;plot(xAxis,mean_nomiss);
hold on; plot(xAxis,mean_miss);
legend('No miss', 'Miss');
title('Pupil diameter aligned to cue time')

%% running z-score
xSession = -1:0.05:5.95;
pupilDia_session = zeros(length(trialData.cueTimes),140);
for kk = 1:length(trialData.cueTimes)
    timeLower = trialData.cueTimes(kk)-1;
    timeHigher = trialData.cueTimes(kk) + 6;
    frameIndex = frameTime >= timeLower & frameTime <= timeHigher;
    pupilDia_session(kk,:) = dia_runningz(frameIndex);
end
mean_session = nanmean(pupilDia_session);
figure; plot(xSession, mean_session);

xAxis = -3:0.05:2.95;
% align to the cue time
pupilDia_cue = zeros(length(trialData.cueTimes), 120);

for kk = 1:length(trialData.cueTimes)
    timeLower = trialData.cueTimes(kk)-3;
    timeHigher = trialData.cueTimes(kk) + 3;
    frameIndex = frameTime >= timeLower & frameTime <= timeHigher;
    pupilDia_cue(kk,:) = dia_runningz(frameIndex);
end

mean_cue = nanmean(pupilDia_cue);
figure;plot(xAxis,mean_cue);
%s = shadedErrorBar(xAxis,pupilDia_cue,{@nanmean,@nanstd},'lineprops','-black','patchSaturation',0.2);
%set(s.edge,'LineWidth',1,'LineStyle',':')
title('Pupil diameter aligned to cue time')

% align to the response
pupilDia_rt = zeros(length(trialData.cueTimes), 120);

for kk = 1:length(trialData.cueTimes)
    if isnan(trialData.rt(kk))
        pupilDia_rt(kk,:) = NaN;
    else
        timeLower = trialData.rt(kk)-3;
        timeHigher = trialData.rt(kk) + 3;
        frameIndex = frameTime >= timeLower & frameTime <= timeHigher;
        pupilDia_rt(kk,:) = dia_runningz(frameIndex);
    end
end

mean_rt = nanmean(pupilDia_rt);
figure;plot(xAxis, mean_rt);
title('Pupil diameter aligned to response time')
% align to the reward

pupilDia_reward = zeros(length(trialData.cueTimes),120);

for kk = 1:length(trialData.cueTimes)
    if isnan(trialData.rt(kk))
        pupilDia_reward(kk,:) = NaN;
    else
        timeLower = trialData.outcomeTimes(kk)-3;
        timeHigher = trialData.outcomeTimes(kk) + 3;
        frameIndex = frameTime >= timeLower & frameTime <= timeHigher;
        pupilDia_reward(kk,:) = dia_runningz(frameIndex);
    end
end
mean_reward = nanmean(pupilDia_reward(trials.reward,:));
mean_noreward = nanmean(pupilDia_reward(trials.noreward,:));
mean_miss = nanmean(pupilDia_reward(trials.miss,:));

figure;plot(xAxis,mean_reward);
hold on; plot(xAxis,mean_noreward);
legend('Reward','No reward');
title('Pupil diameter aligned to outcome time')

% separate into left/right choice/reward

% left/right choice
mean_left = nanmean(pupilDia_rt(trials.left,:));
mean_right = nanmean(pupilDia_rt(trials.right,:));
figure;plot(xAxis,mean_left);
hold on; plot(xAxis,mean_right);
legend('Left choice', 'Right choice');
title('Pupil diameter aligned to response time')

% left/right reward
mean_rewardleft = nanmean(pupilDia_reward(trials.reward & trials.left,:));
mean_rewardright = nanmean(pupilDia_reward(trials.reward & trials.right,:));
figure;plot(xAxis,mean_rewardleft);
hold on; plot(xAxis,mean_rewardright);
legend('Left reward', 'Right reward');
title('Pupil diameter aligned to outcome time')

% missed trial, aligned to the cue time
mean_nomiss = nanmean(pupilDia_cue(~trials.miss,:));
mean_miss = nanmean(pupilDia_cue(trials.miss,:));
figure;plot(xAxis,mean_nomiss);
hold on; plot(xAxis,mean_miss);
legend('No miss', 'Miss');
title('Pupil diameter aligned to cue time')
% regression

%% spectrum

% get the trial mask
trialMask = zeros(1, length(dia_runningz));
for kk = 1:length(trialData.triggerTimes)
    if kk == 1
        trialMask(frameTime < trialData.triggerTimes(kk)) = 0;
    elseif kk == length(trialData.triggerTimes)
        trialMask(frameTime > trialData.triggerTimes(kk)) = kk;
    else
        trialMask(frameTime < trialData.triggerTimes(kk) & frameTime > trialData.triggerTimes(kk - 1)) = kk-1;
    end
end
    
data = [dia_runningz; trialMask];
params.tapers = [5,9];
params.Fs = 20;
params.err = [1, 0.05];
[S,f,Serr] = mtspectrumc(diaHor2, params );

figure; plot(f(1:10),S(1:10))

movingwin = [5 0.05];
params.Fs = 20;
params.trialave = 1;
[S,t,f] =mtspecgramc(data', movingwin, params );

% align this to the time
t = frameTime(1:length(S));
xSession = -1:0.05:5.95;
pupilDia_F = zeros(length(trialData.cueTimes)-1,140, length(f));
for kk = 1:length(trialData.cueTimes)-1
    timeLower = trialData.cueTimes(kk)-1;
    timeHigher = trialData.cueTimes(kk) + 6;
    frameIndex = t >= timeLower & t <= timeHigher;
    pupilDia_F(kk,:,:) = S(frameIndex,:);
end
mean_F = nanmean(pupilDia_F);
figure; imagesc(squeeze(mean_F));

% Plot the Fourier transform
figure;
Fs = 20;
diaTask = fillmissing(dia_runningz(trialMask ~= 0),'linear');
freqs = Fs*(0:(length(diaTask)/2))/length(diaTask);
tsfft = abs(fft(diaTask));
tsfft = tsfft(1:length(diaTask)/2+1)/length(diaTask);
plot(freqs, tsfft)
xlabel("Frequenzy (hz)")
ylabel("Power")

%% align with the trial Data
%% analyze the data ( with another model, not useful for now)

diameterUnit = 'px';
diameter.R = diaHor;
diameter.L = [];
diameter.t_ms = frameTime' * 1000;

dataModel =  RawFileModel(diameterUnit,diameter);
filename = ['E:\data\matching_pennies\874pupil\',matFile(1).name(1:38),'_testanalysis.mat'];
dataModel.saveMatFile(filename)

filepath = 'E:\data\matching_pennies\874pupil\';
filename = [matFile(1).name(1:38),'_testanalysis.mat'];
settingsOut = PupilDataModel.getDefaultSettings();
settingsOut.raw.PupilDiameter_Min = 0;
settingsOut.raw.PupilDiameter_Max = 100;
settingsOut.valid.interp_upsamplingFreq = 1;
%settingsOut.raw.
hPupilData = PupilDataModel(filepath,filename, settingsOut);

% try filter first
hPupilData.filterRawData();