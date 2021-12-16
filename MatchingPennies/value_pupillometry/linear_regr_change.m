function [ output ] = linear_regr_change ( signal,event, trialSubset, params)
% % linear_regr %
%% change the regression equation to Sul et al.2011
%PURPOSE:   Multiple lienar regression for pupil change
%
%INPUT ARGUMENTS
%   signal:        pupil change for every trial

%   event:          event, dummy-coded (e.g., for choice, left=-1, right=1)
%                   %currently can handle up to 2 types of events, e.g., choice and outcome
%   trialSubset:    the subset of trials to investigate, all else set to NaN
%   params.window:  the time window around which to align signal
%   params.nback:   consider events from up to this # trials back


%OUTPUT ARGUMENTS
%   output:         structure containing regression analysis results
%
% To plot the output, use MP_plot_regr().


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


%% multiple linear regression with moving window 
% to see which behavioral events explain signal variability

warning('off','MATLAB:singularMatrix');
warning('off','stats:pvaluedw:ExactUnavailable');
    
try
    stats = regstats(signal, factors, terms);
    
    for kk=1:size(terms,1)
        output.coeff(kk)=stats.beta(kk);        %coefficient
        output.pval(kk)=stats.tstat.pval(kk);   %pvalue for coefficient
    end
catch
    for kk=1:size(terms,1)
        output.coeff(kk)=NaN;
        output.pval(kk)=NaN;
    end
end
warning('on','MATLAB:singularMatrix');
warning('on','stats:pvaluedw:ExactUnavailable');

end
