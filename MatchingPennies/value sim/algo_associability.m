function stats=algo_associability(stats,params,x)
% % algo_associability %
%PURPOSE:   Simulate player based on Pearce-Hall mechanism for associability-gated learning
%           (see Li et at., Nature Neuroscience 2011)
%AUTHORS:   Hongli Wang 191015
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       kappa - learning rate coefficient
%       eta - learning rate to update associabitity
%       alpha0 - initial associability
%   x:      which player? (1 or 2)
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

if stats.currTrial == 1  %if this is the first trial
    stats.ql(stats.currTrial,x) = 0;
    stats.qr(stats.currTrial,x) = 0;
    stats.rpe(stats.currTrial,x) = NaN;
    stats.pl(stats.currTrial,x) = 0.5;
    stats.alpha(stats.currTrial,x) = params.alpha0;
else
    %% update action values

    if stats.c(stats.currTrial-1,x)==-1   % if chose left on last trial
        stats.rpe(stats.currTrial-1,x)=stats.r(stats.currTrial-1,x)-stats.ql(stats.currTrial-1,x);
        stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x)+params.kappa*stats.alpha(stats.currTrial-1,x)*stats.rpe(stats.currTrial-1,x);
        stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x);
        stats.alpha(stats.currTrial,x) = (1-params.eta)*stats.alpha(stats.currTrial-1,x) + params.eta * abs(stats.rpe(stats.currTrial-1,x));
    elseif stats.c(stats.currTrial-1,x)==1   % else, chose right
        stats.rpe(stats.currTrial-1,x)=stats.r(stats.currTrial-1,x)-stats.qr(stats.currTrial-1,x);
        stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x);
        stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x)+params.kappa*stats.alpha(stats.currTrial-1,x)*stats.rpe(stats.currTrial-1,x);
        stats.alpha(stats.currTrial,x) = (1-params.eta)*stats.alpha(stats.currTrial-1,x) + params.eta * abs(stats.rpe(stats.currTrial-1,x));
    else  %miss trials
        stats.rpe = 0;
        stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x);
        stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x);
        stats.alpha = stats.alpha;
    end
    
    %% softmax rule for action selection
    stats.pl(stats.currTrial,x)=1/(1+exp(-(stats.ql(stats.currTrial,x)-stats.qr(stats.currTrial,x))));
    
end

end
