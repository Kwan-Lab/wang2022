function [ output ] = MP_get_psth ( signal, t, eventTime, psth_label, params)
% % MP_get_psth %
%PURPOSE:   Given a time-series signal and event times, find peri-event signal
%           by binning, and estimate CI using bootstrap
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   signal:         time-series signal (e.g., dF/F)
%   t:              time points corresponding to the signal
%   eventTime:      the event times
%   psth_label:     text saying what the event is about
%   params.window:   the time window around which to align signal
%   params.numBootstrapRepeat:  number of bootstrap repeats
%   params.CI:                  confidence interval
%
%OUTPUT ARGUMENTS
%   output:         structure containing peri-event time histogram
%
% To plot the output, use MP_plot_psth().

%%
if isfield(params,'numBootstrapRepeat')
    numRep=params.numBootstrapRepeat;
    CI=params.CI;
else
    numRep=0;
    CI=NaN;
end

output.psth_label=psth_label;

window=params.window;
window_dt=nanmean(diff(window));
output.t=window(1:end-1)+window_dt/2; %use center of bin as the time associated with the average signal
output.t=output.t(:);

nCells = 1;
output.signal=nan(numel(output.t),nCells);

%% find average (by binning/histogram)
%find the segment of signal and its relative time, around each event
tempSig=[]; tempTime=[]; tempEvent=[];
for j=1:numel(eventTime)
    relTime=t-eventTime(j);     %time relative to the event
    
    tempSig=[tempSig; signal(relTime>=window(1) & relTime<=window(end))];
    tempTime=[tempTime; relTime(relTime>=window(1) & relTime<=window(end))];
    tempEvent=[tempEvent; j*ones(sum(relTime>=window(1) & relTime<=window(end)),1)];    %the event# associated with this signal segment

    %mean value for entire window for each trial
    output.valbyTrial(j)=nanmean(signal(relTime>=window(1) & relTime<=window(end)));
end

for j=1:numel(window)-1 %for each bin of time, what is the average signal?
    output.signal(j)=nanmean(tempSig(tempTime>=window(j) & tempTime<(window(j)+window_dt)));    
end

%% find CI using the bootstrap method

if ~isnan(CI)
    bootSig=[];
    for i=1:numRep
        
        %draw a subset for bootstrap
        drawIndex=randsample(numel(eventTime),numel(eventTime),'true'); %each time draw a set of events with replacement
        drawSig=[]; drawTime=[];
        for k=1:numel(drawIndex)
            drawSig=[drawSig; tempSig(tempEvent==drawIndex(k))];
            drawTime=[drawTime; tempTime(tempEvent==drawIndex(k))];
        end
        
        %go through each time bin, and average
        for j=1:numel(window)-1
            bootSig(j,i)=nanmean(drawSig(drawTime>=window(j) & drawTime<(window(j)+window_dt)));
        end
    end
    
    %bootstrap mean and confidence interval
    output.bootSig = bootSig;
    output.bootavg=squeeze(nanmean(bootSig,2));
    output.bootlow=prctile(bootSig,0.5*(1-CI)*100,2);
    output.boothigh=prctile(bootSig,(1-0.5*(1-CI))*100,2);
end

end