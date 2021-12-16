function stats=algo_Q_RPE(stats,params)
% % algo_Q_RPE %
%PURPOSE:   Simulate player based on q-learning with reward predictione
%           errors (see Sul et at., Neuron 2010)
%AUTHORS:   AC Kwan 180423
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       params(1) = a - learning rate
%       params(2) = b - inverse temperature
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

a = params(1);
b = params(2);

if stats.currTrial == 1  %if this is the first trial
    stats.ql(stats.currTrial) = 0.5;
    stats.qr(stats.currTrial) = 0.5;
    stats.rpe(stats.currTrial) = NaN;
    stats.pl(stats.currTrial) = 0.5;
    
else
    %% update action values

    if stats.c(stats.currTrial-1)==-1   % if chose left on last trial
        stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.ql(stats.currTrial-1);
        stats.ql(stats.currTrial)=stats.ql(stats.currTrial-1)+a*stats.rpe(stats.currTrial-1);
        stats.qr(stats.currTrial)=stats.qr(stats.currTrial-1);
    elseif stats.c(stats.currTrial-1)==1   % else, chose right
        stats.rpe(stats.currTrial-1)=stats.r(stats.currTrial-1)-stats.qr(stats.currTrial-1);
        stats.ql(stats.currTrial)=stats.ql(stats.currTrial-1);
        stats.qr(stats.currTrial)=stats.qr(stats.currTrial-1)+a*stats.rpe(stats.currTrial-1);
    else  %miss trials
        stats.rpe(stats.currTrial-1) = 0;
        stats.ql(stats.currTrial)=stats.ql(stats.currTrial-1);
        stats.qr(stats.currTrial)=stats.qr(stats.currTrial-1);
    end
    
    %% softmax rule for action selection
    stats.pl(stats.currTrial)=1/(1+exp(-b*(stats.ql(stats.currTrial)-stats.qr(stats.currTrial))));
    
end

end
