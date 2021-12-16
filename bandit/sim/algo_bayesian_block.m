function stats=algo_bayesian_block(stats,params)
% % algo_bayesian_block %
%PURPOSE:   Simulate player based on a bayeisan model
%           --- Agent has full knowledge of block structure of task
%           Agent knows when criterion for probabilities switch is met 
%           (which is unrealistic because it is not fully observable, but
%           it is safe to assume agent could infer from other observable,
%           (e.g., # reward consecutively obtained from same side, or bias
%           in recent choices) and here in this model we won't explicitly 
%           asked how the agent inferred that information
%
%AUTHORS:   H Wang and AC Kwan 200210
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%           1) p = probability of switch after criterion is met
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

p_crit = params(1);

probreward.high = 0.7;
probreward.low = 0.1;
criterion = 10;

if stats.currTrial == 1  %if this is the first trial
    stats.ql(stats.currTrial) = 0.5;  % ql here refers to the probability that the current trial is in state 1
    stats.qr(stats.currTrial) = 0.5;
    stats.pl(stats.currTrial) = 0.5;
    
else
    %% counting the number of correct choices (choosing the high-probability side) since the last switch
    
    % rule from ground truth
    rule = stats.rule;
    
    % how many trials since the last change in inferred current state
    x = rule';
    l = diff([0 find(x(1:end-1) ~= x(2:end)) length(x)]); %run-length encoding
    swTrial = cumsum(l);                    %trial where there were switches
    swIdx = sum(swTrial < stats.currTrial); %index corresponding to the last switch (relative to current trial)
    
    if swIdx > 0
        startTrial = swTrial(swIdx)+1;  %if there was a last switch, start counting form last switch
    else
        startTrial = 1;               %if last switch could not be found, start counting from trial 1
    end
    
    % tally of number of correct choices since last switch
    correctChoice = sum(stats.hr_side(startTrial:stats.currTrial) == stats.c(startTrial:stats.currTrial));
    
    %% update the posterior probability given choice and criterion
    % b_i(t+1) ~ p(o(t)|i) * sum_j (p(i|j)b_j(t))
    %
    % b_i(t+1) is posterior probability that animal is in sub-state i at time t+1
    % p(o(t)|i) is likelihood of observation o(t) under hypothetical sub-state i
    % p(i|j) is the probability of transitioning from sub-state j to sub-state i
    % b_j(t), summing over probabilities of being in sub-state j at time t
    %
    % see Starkweather, ... Uchida, Gershman, Nat Neurosci, 2017
    
    if correctChoice < criterion  % only recently switched probabilities, no further switch expected
        p = 0;  % criterion not met, no transition expected
    else
        p = p_crit;    % criterion met, transition expected according to exponential distribution
    end
    
    if stats.c(stats.currTrial-1) == -1 % choose left
        if stats.r(stats.currTrial-1) == 1 % get reward
            ql = probreward.high * ((1-p) * stats.ql(stats.currTrial-1) + p * stats.qr(stats.currTrial-1));
            qr = probreward.low  * (p * stats.ql(stats.currTrial-1) + (1-p) * stats.qr(stats.currTrial-1));
        else % no reward
            ql = (1 - probreward.high) * ((1-p) * stats.ql(stats.currTrial-1) + p * stats.qr(stats.currTrial-1));
            qr = (1 - probreward.low)  * (p * stats.ql(stats.currTrial-1) + (1-p) * stats.qr(stats.currTrial-1));
        end
        
    elseif stats.c(stats.currTrial-1) == 1  % choose right
        if stats.r(stats.currTrial-1) == 1  % get reward
            ql = probreward.low  * ((1-p) * stats.ql(stats.currTrial-1) + p * stats.qr(stats.currTrial-1));
            qr = probreward.high * (p * stats.ql(stats.currTrial-1) + (1-p) * stats.qr(stats.currTrial-1));
        else % no reward
            ql = (1 - probreward.low)  * ((1-p) * stats.ql(stats.currTrial-1) + p * stats.qr(stats.currTrial-1));
            qr = (1 - probreward.high) * (p * stats.ql(stats.currTrial-1) + (1-p) * stats.qr(stats.currTrial-1));
        end
        
    else % miss
        ql = ((1-p) * stats.ql_bayesian(stats.currTrial-1) + p * stats.qr_bayesian(stats.currTrial-1));
        qr = (p * stats.ql_bayesian(stats.currTrial-1) + (1-p) * stats.qr_bayesian(stats.currTrial-1));
        
    end
    
    stats.ql(stats.currTrial) = ql / (ql + qr); %probabilities of sub-states must add to 1
    stats.qr(stats.currTrial) = qr / (ql + qr);
    
    % assume if in sub-state left, then choose left always
    % if in sub-state right, then choose right always
    % for mixture, then choose probabilistically accordingly
    stats.pl(stats.currTrial) = stats.ql(stats.currTrial);
    
end

end
