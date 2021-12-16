function [ output ] = linear_regr_perAnimal ( signal, t, event, eventTime, params)
% % linear_regr %
%PURPOSE:   Multiple lienar regression
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   signal:         time-series signal (e.g., dF/F)
%   t:              time points corresponding to the signal
%   event:          event, dummy-coded (e.g., for choice, left=-1, right=1)
%                   %currently can handle up to 2 types of events, e.g., choice and outcome
%   eventTime:      the event times
%   trialSubset:    the subset of trials to investigate, all else set to NaN
%   params.window:  the time window around which to align signal
%   params.nback:   consider events from up to this # trials back
%
%OUTPUT ARGUMENTS
%   output:         structure containing regression analysis results
%
% To plot the output, use MP_plot_regr().

%% interpolate signal

window=params.window;

sigbyTrial_all = [];
Ind_all = [];
startIndex = zeros(1,length(eventTime));  % get the start Index of each session

Ind = 1:length(event);
% generate event mask
summation = 1;
eventMask = zeros(1, size(event,1));
for kk = 1:length(eventTime)
    startIndex(kk) = summation;
    eventMask(summation:summation+length(eventTime{kk})-1) = kk;
    summation = summation+length(eventTime{kk});
end
for ii = 1:length(eventTime)
    % interpolate the signal to a finer time scale
    if ~isempty(signal{ii})
        interdt=0.01;
        intert=[t{ii}(1):interdt:t{ii}(end)]';
        intersig=interp1(t{ii},signal{ii},intert);

% align signal to the event
% use window slightly wider than the regression, so regression analysis
% won't run into the boundaries of this variable

% since there are multiple sessions, then align the signal session by
% session
    
        [sigbyTrial, tbyTrial]=align_signal(intert,intersig,eventTime{ii},[window(1)-1 window(end)+1]);
        sigbyTrial_all = [sigbyTrial_all,sigbyTrial];
        Ind_all = [Ind_all,Ind(eventMask==ii)];
    end
    
end

% get trialIndex that has coresponding pupil data

%% the factors in the linear regression

output.numPredictor=size(event,2);  %number of predictor
output.nback = params.nback;

factors=[];
for j=1:output.numPredictor
    event(~Ind_all,j)=nan;        %if event is not part of trialSubset to analyze, set to NaN
    factors = [factors event(Ind_all,j)];
    
    if output.numPredictor == 1 || output.numPredictor == 2
        for k=1:output.nback
            % adding 2 NaNs before each session
            event_kback = [];
            for ll = 1:length(eventTime)
                if ~isempty(signal{ll})
                    event_kback=[event_kback;[NaN(k,1); event((startIndex(ll):startIndex(ll) + length(eventTime{ll})-k-1),j)]];     %event k-back
                end
            end
            factors = [factors event_kback];
        end
    end
end

if output.numPredictor==1
    terms=[zeros(1,1+output.nback) ; eye(1+output.nback)]; %bias, c(n), c(n-1), c(n-2), c(n-3), ...
elseif output.numPredictor==2
    if params.interaction == true
        %bias, c(n), c(n-1), c(n-2), ..., r(n), r(n-1), r(n-2), ..., 
        terms=[zeros(1,(1+output.nback)*output.numPredictor) ; eye((1+output.nback)*output.numPredictor)];  
        %add interaction terms: c(n)xr(n), c(n-1)xr(n-1), ...
        terms=[terms; [eye(1+output.nback) eye(1+output.nback)]];
    else
        %bias, c(n), c(n-1), c(n-2), ..., r(n), r(n-1), r(n-2), ...
        terms=[zeros(1,(1+output.nback)*output.numPredictor) ; eye((1+output.nback)*output.numPredictor)];
    end
elseif output.numPredictor > 2 % this is the linear regression for RL latent variables
    disp('Linear regression model for RL latent variables');
    terms = [zeros(1, output.numPredictor);eye(output.numPredictor,output.numPredictor)];
end

%% parameters of the regression
step_dur = nanmean(diff(window));
step_size = step_dur;       %if step size is step duration, then doing this in non-overlapping windows
output.regr_time=[window(1)+step_dur/2:step_size:window(end)-step_dur/2]';
        
%% multiple linear regression with moving window 
% to see which behavioral events explain signal variability

warning('off','MATLAB:singularMatrix');
warning('off','stats:pvaluedw:ExactUnavailable');
for jj=1:numel(output.regr_time)
    idx1=sum(tbyTrial<=output.regr_time(jj)-step_dur/2);     %signal corresponding to current time step
    idx2=sum(tbyTrial<(output.regr_time(jj)+step_dur/2));
    tempsig=squeeze(nanmean(sigbyTrial_all(idx1:idx2,:),1));
    
    try
        stats=regstats(tempsig',factors,terms);
        for kk=1:size(terms,1)
            output.coeff(jj,kk)=stats.beta(kk);        %coefficient
            output.pval(jj,kk)=stats.tstat.pval(kk);   %pvalue for coefficient
        end
    catch
        for kk=1:size(terms,1)
            output.coeff(jj,kk)=NaN;
            output.pval(jj,kk)=NaN;
        end
    end
end
warning('on','MATLAB:singularMatrix');
warning('on','stats:pvaluedw:ExactUnavailable');

end
