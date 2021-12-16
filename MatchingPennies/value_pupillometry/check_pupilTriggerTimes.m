function check_pupilTriggerTimes (stackInfo, trialData, savefigpath)
% % check_pupilTriggerTimes %
%PURPOSE:   Plot the pupil file start times relative to behavioral trigger
%           times
%AUTHORS:   Hongli Wang
%
%INPUT ARGUMENTS
%   data:     extracted pupil diameter in csv file
%   stackInfo:   Time stamps for when Presentation sent TTL pulse to
%                   matlab
%   logfileTimes: Time stamps from presentation logfile

    
        % calculate the frames
        
logFrameTrial = zeros(1, length(trialData.cueTimes));
        
for ii = 1:length(trialData.cue)
    if ii < length(trialData.cue)
        if isfield(trialData, 'triggerTimes')
            trialTime = trialData.triggerTimes(ii+1) - trialData.triggerTimes(ii);
        else
            trialTime = trialData.cueTimes(ii+1) - trialData.cueTimes(ii);
        end
        logFrameTrial(ii) = trialTime * 20;
    else
        trialTime = 0; % skip it for now
    end
end
        
% get the trial frames from matlab
matFrameTrial = stackInfo.framePerTrial(stackInfo.framePerTrial~=0);
% ignore 2 for now (and other small numbers
matFrameTrial_no2 = matFrameTrial(matFrameTrial ~= 2 & matFrameTrial ~= 1 & matFrameTrial ~= 3);

maxLength = min(length(logFrameTrial), length(matFrameTrial_no2));
figure;
scatter(logFrameTrial(1:maxLength-1), matFrameTrial_no2(1:maxLength-1), 'black', 'filled');
hold on; plot([0 500], [0 500], 'lineWidth',1);
axis square;
xlabel('Frames per trial for triggers according to .log file');
ylabel('Frames per trial for matlab logs');
title('Value should not deviate from diagonal');

if ~exist(savefigpath)
    mkdir(savefigpath)
end
print(gcf,'-dpng',[savefigpath,'\check-triggertiming']);    %png format
saveas(gcf, [savefigpath,'/check-triggertiming'], 'fig');

close;
end