function stats=algo_DQ_RPE(stats,params,x)
% % algo_DQ_RPE %
%PURPOSE:   Simulate player based on q-learning with reward predictione
%           errors, with differential learning rates
%AUTHORS:   AC Kwan 180423
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       params(1) = a - learning rate for choice yielding reward
%       params(2) = b - inverse temperature
%       params(3) = a2 - learning rate for choice yielding no reward
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

a = params.a;
b = params.b;
a2 = params.a2;

if stats.currTrial == 1  %if this is the first trial
    stats.ql(stats.currTrial,x) = 0.5;
    stats.qr(stats.currTrial,x) = 0.5;
    stats.rpe(stats.currTrial,x) = NaN;
    stats.pl(stats.currTrial,x) = 0.5;
    
else
    %% update action values
    
    if stats.c(stats.currTrial-1,x)==-1   % if chose left on last trial
        stats.rpe(stats.currTrial-1,x)=stats.r(stats.currTrial-1,x)-stats.ql(stats.currTrial-1,x);
        if stats.r(stats.currTrial-1,x) > 0
            stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x)+a*stats.rpe(stats.currTrial-1,x);
        else
            stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x)+a2*stats.rpe(stats.currTrial-1,x);
        end
        stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x);
    elseif stats.c(stats.currTrial-1,x)==1   % else, chose right
        stats.rpe(stats.currTrial-1,x)=stats.r(stats.currTrial-1,x)-stats.qr(stats.currTrial-1,x);
        stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x);
        if stats.r(stats.currTrial-1,x) > 0
            stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x)+a*stats.rpe(stats.currTrial-1,x);
        else
            stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x)+a2*stats.rpe(stats.currTrial-1,x);
        end
    else  %miss trials
        stats.rpe(stats.currTrial-1,x) = 0;
        stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x);
        stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x);
    end
    
    %% softmax rule for action selection
    stats.pl(stats.currTrial,x)=1/(1+exp(-b*(stats.ql(stats.currTrial,x)-stats.qr(stats.currTrial,x))));
    
end

end
