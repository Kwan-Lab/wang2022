function stats=algo_logreg_RCUC(stats,params,x)
% % algo_logreg_RCUC %
%PURPOSE:   Simulate a player based on logistic regression
%           regressors are past rewarded choices and unrewarded choices
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       bias - bias term
%       RC - reward/choice terms
%       UC - unreward/choice terms
%   x:      which player? (1 or 2)
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

% make sure the params are column vectors
b_reward=params.RC(:);
b_unreward=params.UC(:);

% how many trials to look back for logistic regression
step_back = numel(params.RC);

%probability of choosing left based on logistic regression
if stats.currTrial == 1   %if this is the first trial
    stats.pl(1,x) = 0.5;  %choose randomly
else
    if stats.currTrial>step_back
        backIdx = step_back;  %if at least number of trials permit, use step_back
    else
        backIdx = stats.currTrial-1;        %else use what is available
    end
    
    %choice and reward from the past trials under consideration
    c = stats.c(stats.currTrial-1:-1:stats.currTrial-backIdx,x);
    r = stats.r(stats.currTrial-1:-1:stats.currTrial-backIdx,x);
    
    % create regressor vectors
    YR=zeros(length(c),1);   % rewarded choice: rewarded/right = 1; rewarded/left = -1
    YR=r.*((c==1) & (r>0)) + (-1*r).*((c==-1) & (r>0));  %allows for reward size r!=1
%    YR=1*((c==1) & (r==1)) + (-1).*((c==-1) & (r==1));
    NR=zeros(length(c),1);   % unrewarded choice: unrewarded/right = 1; unrewarded/left = -1
    NR=1*((c==1) & (r==0)) + (-1)*((c==-1) & (r==0));
    
    b_temp = params.bias;    %bias term
    b_temp = b_temp + nansum(b_reward(1:backIdx).*YR);    %reward/choice terms
    b_temp = b_temp + nansum(b_unreward(1:backIdx).*NR);  %unreward/choice terms
    stats.pl(stats.currTrial) = 1 - exp(b_temp)/(1+exp(b_temp));  %p(left) = 1 - p(right)
end

end
