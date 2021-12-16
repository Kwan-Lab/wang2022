function [ output ] = linear_regr_PR ( signal, t, event, eventTime, trialSubset, params, cueTimes)
% % linear_regr %
%% change the regression equation to Sul et al.2011
%PURPOSE:   Multiple lienar regression, for pupil response specifically
%           need to deal with the overlapped areaPR
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
%   cueTimes:       start time of every trial for signal alignment


%OUTPUT ARGUMENTS
%   output:         structure containing regression analysis results
%
% To plot the output, use MP_plot_regr().

%% interpolate signal

window=params.window;

% interpolate the signal to a finer time scale
interdt=0.02;
for tt=1:size(t,1)
    intert(tt,:) = [t(tt,1):interdt:t(tt,size(t,2))]';
    intersig(tt,:) = interp1(t(tt,:),signal(tt,:),intert(tt,:));
end
% t might contain similar time point
%[t,index] = unique(t);
% intersig=interp1(t,signal,intert);
%intersig = interp1(t,signal,intert);
% align signal to the event
% use window slightly wider than the regression, so regression analysis
% won't run into the boundaries of this variable
[sigbyTrial, tbyTrial]=align_signal_PR(intert,intersig,eventTime,[window(1)-1 window(end)+1], cueTimes);

%% the factors in the linear regression

output.numPredictor=size(event,2);  %number of predictor
output.nback = params.nback;
output.interaction = params.interaction;

factors=[];
for j=1:output.numPredictor
    event(~trialSubset,j)=nan;        %if event is not part of trialSubset to analyze, set to NaN
    factors = [factors event(:,j)];
    
    if output.numPredictor == 1 || output.numPredictor == 2
        for k=1:output.nback
            event_kback=[NaN(k,1); event(1:end-k,j)];     %event k-back
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
% elseif output.numPredictor ==3
%     if params.interaction == true
%         %bias, c(n), c(n-1), c(n-2), ..., r(n), r(n-1), r(n-2), ..., etc.
%         terms=[zeros(1,(1+output.nback)*output.numPredictor) ; eye((1+output.nback)*output.numPredictor)];   
%         %add interaction terms: c(n)xr(n), c(n-1)xr(n-1), ...
%         terms=[terms; [eye(1+output.nback) eye(1+output.nback) zeros(1+output.nback)]; [eye(1+output.nback) zeros(1+output.nback) eye(1+output.nback)]];
%     else
%         %bias, c(n), c(n-1), c(n-2), ..., r(n), r(n-1), r(n-2), ..., etc.
%         terms=[zeros(1,(1+output.nback)*output.numPredictor) ; eye((1+output.nback)*output.numPredictor)];
%     end
% else
%     error('This function does not handle more than 3 predictor events plus their interactions.');
% end
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
    tempsig = [];
    for ss = 1:size(tbyTrial,2)
        idx1=sum(tbyTrial(:,ss)<=(output.regr_time(jj)-step_dur/2));     %signal corresponding to current time step
        idx2=sum(tbyTrial(:,ss)<(output.regr_time(jj)+step_dur/2));
        tempsig = [tempsig; nanmean(sigbyTrial(idx1:idx2,ss))];
    end
    %tempsig=squeeze(nanmean(sigbyTrial(idx1:idx2,:),1));
    
    try
        
        stats=regstats(tempsig',factors,terms);
        mdl = fitlm(factors,tempsig',terms);
        tbl=anova(mdl);
        for kk=1:size(terms,1)
            output.coeff(jj,kk)=stats.beta(kk);        %coefficient
            output.pval(jj,kk)=stats.tstat.pval(kk);   %pvalue for coefficient
            if kk == 1
                output.SumSq(jj,kk) = tbl.SumSq(end);  % last one is error term
            else
                output.SumSq(jj,kk) = tbl.SumSq(kk-1);
            end
        end
        output.mse = stats.mse;
        
    catch
        for kk=1:size(terms,1)
            output.coeff(jj,kk)=NaN;
            output.pval(jj,kk)=NaN;
            output.SumSq(jj,kk) = NaN;
        end
        output.mse = NaN;
    end
end
warning('on','MATLAB:singularMatrix');
warning('on','stats:pvaluedw:ExactUnavailable');

end
