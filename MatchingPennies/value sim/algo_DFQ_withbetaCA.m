function stats=algo_DFQ_withbetaCA(stats,params,x)
% % algo_DFQ %
%PURPOSE:   Simulate player based on q-learning (differential, forgetting,
%           Q, see Ito and Doya, J Neurosci 2009)
%AUTHORS:   AC Kwan 180423
% add a beta parameter % the purpose is to validate the Katahira paper
% CA : choice autocorrelation
% reference Katahira 2015

%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       a1 - learning rate (also = 1 minus the forgetting rate)
%       a2 - related to forgetting rate for action not chosen
%       k1 - strength of reinforcement by reward
%       k2 - strength of aversion from no-reward outcome
%       beta - inverse temperature
%       tau - learning rate for choice autocorrelation
%       phi - coefficient for choice autocorrelation
%   x:      which player? (1 or 2)
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

if stats.currTrial == 1  %if this is the first trial
    stats.ql(stats.currTrial,x) = 0;
    stats.qr(stats.currTrial,x) = 0;
    stats.cl(stats.currTrial,x) = 0;
    stats.cr(stats.currTrial,x) = 0;
    stats.rpe(stats.currTrial,x) = NaN;
    stats.pl(stats.currTrial,x) = 0.5;
    
else
    %% update action values

    if stats.c(stats.currTrial-1,x)==-1     % if chose left on last trial
        stats.cl(stats.currTrial,x) = (1- params.tau) * stats.cl(stats.currTrial-1,x) + params.tau;
        stats.cr(stats.currTrial,x) = (1 - params.tau) * stats.cr(stats.currTrial-1,x);
        if stats.r(stats.currTrial-1,x)>0  %rewarded
            stats.rpe(stats.currTrial,x)=params.k1;
            stats.ql(stats.currTrial,x)=(1-params.a1)*stats.ql(stats.currTrial-1,x)+params.a1*params.k1*stats.r(stats.currTrial-1,x);
            stats.qr(stats.currTrial,x)=(1-params.a2)*stats.qr(stats.currTrial-1,x);
        elseif stats.r(stats.currTrial-1,x)==0  %no reward
            stats.rpe(stats.currTrial,x)=-1*params.k2;
            stats.ql(stats.currTrial,x)=(1-params.a1)*stats.ql(stats.currTrial-1,x)-params.a1*params.k2;
            stats.qr(stats.currTrial,x)=(1-params.a2)*stats.qr(stats.currTrial-1,x);
        end
    elseif stats.c(stats.currTrial-1,x)==1  % else, chose right
        stats.cr(stats.currTrial,x) = (1- params.tau) * stats.cr(stats.currTrial-1,x) + params.tau;
        stats.cl(stats.currTrial,x) = (1 - params.tau) * stats.cl(stats.currTrial-1,x);
        if stats.r(stats.currTrial-1,x)>0  %rewarded
            stats.rpe(stats.currTrial,x)=params.k1;
            stats.ql(stats.currTrial,x)=(1-params.a2)*stats.ql(stats.currTrial-1,x);
            stats.qr(stats.currTrial,x)=(1-params.a1)*stats.qr(stats.currTrial-1,x)+params.a1*params.k1*stats.r(stats.currTrial-1,x);
        elseif stats.r(stats.currTrial-1,x)==0  %no reward
            stats.rpe(stats.currTrial,x)=-1*params.k2;
            stats.ql(stats.currTrial,x)=(1-params.a2)*stats.ql(stats.currTrial-1,x);
            stats.qr(stats.currTrial,x)=(1-params.a1)*stats.qr(stats.currTrial-1,x)-params.a1*params.k2;
        end
    else            %no choice, then just hold on all latent variables
        stats.rpe(stats.currTrial,x)=0;
        stats.ql(stats.currTrial,x)=stats.ql(stats.currTrial-1,x);
        stats.qr(stats.currTrial,x)=stats.qr(stats.currTrial-1,x);
    end
    
    %% softmax rule for action selection
    aL = exp(params.beta*(stats.ql(stats.currTrial,x) + params.phi*stats.cl(stats.currTrial,x)));
    aR = exp(params.beta*(stats.qr(stats.currTrial,x) + params.phi*stats.cr(stats.currTrial,x)));
    stats.pl(stats.currTrial,x)=aL/(aL+aR);
    
end

end
