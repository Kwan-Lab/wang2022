function stats=algo_logreg(stats,params,x)
% % algo_logreg %
%PURPOSE:   Simulate a player based on logistic regression
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       bias - bias term
%       reward - reward terms
%       unreward - unreward terms
%   x:      which player? (1 or 2)
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

% make sure the params are column vectors
b_reward=params.reward(:);
b_unreward=params.unreward(:);

% how many trials to look back for logistic regression
step_back = numel(params.reward);

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
    
    % create regressor vectors - rewarded/right = 1; rewarded/left = -1
    YR=zeros(length(c),1);   % rewarded
    YR=1*((c==1) & (r==1)) + (-1)*((c==-1) & (r==1));
    NR=zeros(length(c),1);   % unrewarded
    NR=1*((c==1) & (r==0)) + (-1)*((c==-1) & (r==0));
    
    b_temp = params.bias;    %bias term
    b_temp = b_temp + nansum(b_reward(1:backIdx).*YR);    %reward/right terms
    b_temp = b_temp + nansum(b_unreward(1:backIdx).*NR);  %unreward/right terms
    stats.pl(stats.currTrial) = 1 - exp(b_temp)/(1+exp(b_temp));  %p(left) = 1 - p(right)
end

end
