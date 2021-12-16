function stats=algo_stochastic(stats,params,x)
% % algo_stochastic %
%PURPOSE:   Simulate a player that makes stochastic choices
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       p  - probability of choosing left
%   x:      which player? (1 or 2)
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

stats.pl(stats.currTrial,x) = params.p;

end
