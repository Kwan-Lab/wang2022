function [sigbyTrial, tbyTrial] = align_signal_PR(t,signal,eventTime,window, cueTimes)
% % align_signal %
%PURPOSE:   Align a time-series signal (e.g., dF/F) to events (e.g. reward)
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   t:              time points associated with the signal vector
%   signal:         signal vector
%   eventTime:      time of the behavioral events (to align to)
%   window:         the time window around which to align
%
%OUTPUT ARGUMENTS
%   sigbyTrial:     signal, a matrix with size [time x event]
%   tbyTrial:       time points associated with first dimension of sigbyTrial

%%
dt=nanmean(diff(t(1,:)));
temptbyTrial = t(1,:);
% the time window around which to align

if isempty(eventTime)   %if there is no event to align to
    sigbyTrial=nan(numel(temptbyTrial),numel(eventTime));
else
    for j=1:numel(eventTime)
        if ~isnan(eventTime(j))  %if there is an event
            % find the trial
            trialInd = find(cueTimes==eventTime(j));
            if trialInd <= size(t,1)
                sigbyTrial(:,j)=signal(trialInd,:);
                tbyTrial(:,j) = t(trialInd,:)-cueTimes(trialInd);
%             else
%                 sigbyTrial(:,j)=nan(numel(temptbyTrial),1);
%                 tbyTrial(:,j) = nan(numel(temptbyTrial),1);
            end
        else   %miss trials, no response
            if trialInd <= size(t,1)
                sigbyTrial(:,j)=nan(numel(temptbyTrial),1);
                tbyTrial(:,j) = t(trialInd,:)-cueTimes(trialInd);
            end
        end
    end
end