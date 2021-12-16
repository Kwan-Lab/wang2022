function stats_sim=predictAgent(player,stats)
% % predictAgent %
%PURPOSE:   Simulate agent based on actual choices and outcomes
%AUTHORS:   AC Kwan 180516
%
%INPUT ARGUMENTS
%   player:    player structure
%       label - the name of the algorithm, the name should correspond to a .m file
%       params - parameters associated with that algorithm
%   stats:      stats about the actual choice behavior
%
%OUTPUT ARGUMENTS
%   stats_sim:      simulated choice behavior and latent variables

%% initialize
rng('shuffle');

c = stats.c(:,1);       %choices during actual behavior
r = stats.r(:,1);       %outcomes during actual behavior
n = numel(c);           %number of trials

stats_sim = []; %stat for the game
stats_sim.currTrial = 1;    % starts at trial number 1
stats_sim.pl=nan(n,1);      % probability to choose left port
stats_sim.c=nan(n,1);       % choice vector
stats_sim.r=nan(n,1);       % reward vector
stats_sim.ql=nan(n,1);      % action value for left choice
stats_sim.qr=nan(n,1);      % action value for right choice
stats_sim.rpe=nan(n,1);     % reward prediction error vector

stats_sim.playerlabel = player.label;
stats_sim.playerparams = player.params;

%take the text label for the player, there should be a corresponding Matlab
%function describing how the player will play
simplay1 = str2func(player.label);

%% simualte the task

for j=1:n
    
    %what is the player's probability for choosing left?
    stats_sim.currTrial = j;
    stats_sim=simplay1(stats_sim,stats_sim.playerparams);
    
    %player chooses - fill in with actual choice
    stats_sim.c(j)=c(j);
    
    %what is the outcome - fill in with actual outcome
    stats_sim.r(j)=r(j);
    
end

end