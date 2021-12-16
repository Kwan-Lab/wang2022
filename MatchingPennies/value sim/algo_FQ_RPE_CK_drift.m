function stats=algo_FQ_RPE_CK_drift(stats,params,x)
% % algo_FQ_RPE_CK %
%PURPOSE:   Simulate player based on q-learning with reward prediction
%           errors, with differential learning rates with choice tendency
%AUTHORS:   H Atilgan & AC Kwan 04302020
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       params(1) = a - learning rate for choice yielding reward
%       params(2) = beta - inverse temperature
%       params(3) = a CK- learning rate for choice tendency
%       params(4) = beta CK - inverse temperature for choice tendency
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step
%%
a = params.a;
beta = params.b;
alpha_c = params.ac;
beta_c = params.bc;
% 
% if stats.currTrial == 1  %if this is the first trial
%     stats.ql(stats.currTrial,x) = stats.ql;
%     stats.qr(stats.currTrial,x) = stats.qr;
%     stats.rpe(stats.currTrial,x) = stats.rpe;
%     stats.pl(stats.currTrial,x) = stats.pl;
%     stats.pr(stats.currTrial,x) = stats.pr;
%     % choice kernels updates
%     stats.ckl(stats.currTrial,x) = stats.ckl;
%     stats.ckr(stats.currTrial,x) = stats.ckr;
    
if stats.currTrial > 1
        %% update choice kernel
    
    if stats.c(stats.currTrial-1,x)==1   % else, chose right
        stats.ckr(stats.currTrial,x) = stats.ckr(stats.currTrial-1,x)+ alpha_c * (1-stats.ckr(stats.currTrial-1,x));
        stats.ckl(stats.currTrial,x) = (1-alpha_c)*stats.ckl(stats.currTrial-1,x);
    elseif stats.c(stats.currTrial-1,x)==-1
        stats.ckr(stats.currTrial,x) = (1-alpha_c)*stats.ckr(stats.currTrial-1,x);
        stats.ckl(stats.currTrial,x) = stats.ckl(stats.currTrial-1,x) + alpha_c * (1-stats.ckl(stats.currTrial-1,x)); % left
    else
        stats.ckl(stats.currTrial,x) = stats.ckl(stats.currTrial-1,x);
        stats.ckr(stats.currTrial,x) = stats.ckr(stats.currTrial-1,x);
    end

    %% update action values   
    if stats.c(stats.currTrial-1,x)==-1   % if chose left on last trial
        stats.rpe(stats.currTrial-1,x)=stats.r(stats.currTrial-1,x)-stats.ql(stats.currTrial-1,x);
        stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x)+a*stats.rpe(stats.currTrial-1,x);
        stats.qr(stats.currTrial,x)=(1-a)*stats.qr(stats.currTrial-1,x);
    elseif stats.c(stats.currTrial-1,x)==1   % else, chose right
        stats.rpe(stats.currTrial-1,x)=stats.r(stats.currTrial-1,x)-stats.qr(stats.currTrial-1,x);
        stats.ql(stats.currTrial,x)=(1-a)*stats.ql(stats.currTrial-1,x);
        stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x)+a*stats.rpe(stats.currTrial-1,x);

    else  %miss trials
        stats.rpe(stats.currTrial-1,x) = 0;
        stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x);
        stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x);
    end
    
    %% softmax with choice kernel
    Q = [stats.qr(stats.currTrial,x) stats.ql(stats.currTrial,x)];
    CK = [stats.ckr(stats.currTrial,x) stats.ckl(stats.currTrial,x)];
    V = beta*Q + beta_c*CK;
    p = exp(V)/sum(exp(V));
    stats.pr(stats.currTrial,x) = p(1);
    stats.pl(stats.currTrial,x) = p(2);
    
end

end