function [sigbyTrial, tbyTrial] = align_signal(t,signal,eventTime,window)
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
dt=nanmean(diff(t));

% the time window around which to align
tbyTrial=dt*[round(window(1)/dt):round(window(end)/dt)]';

if isempty(eventTime)   %if there is no event to align to
    sigbyTrial=nan(numel(tbyTrial),numel(eventTime));
else
    for j=1:numel(eventTime)
        if ~isnan(eventTime(j))  %if there is an event
            eventIdx=sum(eventTime(j)>=t);  %when did the event occur in the time vector
            if (eventIdx+round(window(1)/dt))>0 && (eventIdx+round(window(end)/dt))<=numel(signal) %if not out of bounds
                sigbyTrial(:,j)=signal(eventIdx+round(window(1)/dt):eventIdx+round(window(end)/dt));
            else
                sigbyTrial(:,j)=nan(numel(tbyTrial),1);
            end
        else   %miss trials, no response
            sigbyTrial(:,j)=nan(numel(tbyTrial),1);
        end
    end
end