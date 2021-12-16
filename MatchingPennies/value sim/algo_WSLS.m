function stats=algo_WSLS(stats,params,x)
% % algo_WSLS %
%PURPOSE:   Simulate a player that does win-stay lose-switch
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       param(1) = p  - probability of doing win-stay, lose-switch
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

p = params.p;

if stats.currTrial == 1   %if this is the first trial
    stats.pl(1,x) = 0.5;  %choose randomly
else
    if p >= rand(1)  %follow WSLS
        if stats.r(stats.currTrial-1,x)>0      %if last trial was rewarded
            stats.pl(stats.currTrial,x)=(stats.c(stats.currTrial-1,x)==-1);    %stay: prob(L)=1 if last trial was left, 0 otherwise
        elseif stats.r(stats.currTrial-1,x)==0  %if last trial was unrewarded
            stats.pl(stats.currTrial,x)=1-(stats.c(stats.currTrial-1,x)==-1);  %switch
        else                                    %if last trial was a miss
            stats.pl(stats.currTrial,x)=0.5;                                   
        end
    else  %not following WSLS, then choose stochastically
        stats.pl(stats.currTrial,x)=0.5;
    end
end

end
