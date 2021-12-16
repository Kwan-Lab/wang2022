function stats=simPenniesBiased(player1,player2,n,payoff)
% % simPenniesBiased %
%PURPOSE:   Simulate a biased matching pennies game between two players
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   player1:      player 1
%       label - the name of the algorithm, the name should correspond to a .m file
%       params - parameters associated with that algorithm
%       frac_opto - if a fraction of the trials have optogenetic perturbation
%       opto - parameters associated with that algorithm in perturbation trials
%   player2:      player 2
%   n:            number of trials to simulate
%   payoff:       the payoff matrix [2x2]
%
%OUTPUT ARGUMENTS
%   stats:      simulated choice behavior and latent variables

%% initialize 
rng(1);

stats = []; %stat for the game
stats.currTrial = 1;    % starts at trial number 1
stats.pl=nan(n,2);      % probability to choose left port
stats.c=nan(n,2);       % choice vector
stats.r=nan(n,2);       % reward vector
stats.ql=nan(n,2);      % action value for left choice
stats.qr=nan(n,2);      % action value for right choice
stats.rpe=nan(n,2);     % reward prediction error vector
stats.opto=nan(n,2);    % optogenetic perturbation?

stats.playerlabel{1} = player1.label;
stats.playerlabel{2} = player2.label;
stats.playerparams{1} = player1.params;
stats.playerparams{2} = player2.params;

%if optogenetic perturbation is specified, then what are the modified params
if isfield(player1,'frac_opto')
    if player1.frac_opto > 0
        stats.playerparams_opto{1}=stats.playerparams{1};  %copy from normal condition
        fnames=fieldnames(player1.opto);
        for j=1:numel(fnames)               %if there is a value associated with perturbation
            stats.playerparams_opto{1}.(fnames{j})=player1.opto.(fnames{j}); %replace with the perturbed value
        end
    end
end
if isfield(player2,'frac_opto')
    if player2.frac_opto > 0
        stats.playerparams_opto{2}=stats.playerparams{2};  %copy from normal condition
        fnames=fieldnames(player2.opto);
        for j=1:numel(fnames)               %if there is a value associated with perturbation
            stats.playerparams_opto{2}.(fnames{j})=player2.opto.(fnames{j}); %replace with the perturbed value
        end
    end
end

%take the text label for the player, there should be a corresponding Matlab
%function describing how the player will play
simplay1 = str2func(player1.label);
simplay2 = str2func(player2.label);

%% simualte the game
for j=1:n
    
    %what is the player's probability for choosing left?
    %simulate each player's algorithm
    stats.currTrial = j;
    if isfield(player1,'frac_opto')
        if player1.frac_opto > rand(1)
            stats.opto(j,1)=1; %opto stim trial
            stats=simplay1(stats,stats.playerparams_opto{1},1);
        else
            stats.opto(j,1)=0; %no opto stim trial
            stats=simplay1(stats,stats.playerparams{1},1);
        end
    else                       %no perturbation in general
        stats=simplay1(stats,stats.playerparams{1},1);
    end
    if isfield(player2,'frac_opto')
        if player1.frac_opto > rand(1)
            stats.opto(j,2)=1; %opto stim trial
            stats=simplay1(stats,stats.playerparams_opto{2},2);
        else
            stats.opto(j,2)=0; %no opto stim trial
            stats=simplay1(stats,stats.playerparams{2},2);
        end
    else                       %no perturbation in general
        stats=simplay2(stats,stats.playerparams{2},2);
    end

    %each player chooses
    for k=1:2
        if stats.pl(j,k)>rand(1)
            stats.c(j,k)=-1;
        else
            stats.c(j,k)=1;
        end
    end
    
    %what is the outcome for each player
    if stats.c(j,1) == -1 && stats.c(j,2) == -1  %left, left
        stats.r(j,1) = payoff(1,1);  
        stats.r(j,2) = -2;
    elseif stats.c(j,1) == -1 && stats.c(j,2) == 1  %left, right
        stats.r(j,1) = payoff(1,2); 
        stats.r(j,2) = -1;
    elseif stats.c(j,1) == 1 && stats.c(j,2) == -1  %right, left
        stats.r(j,1) = payoff(2,1);  
        stats.r(j,2) = 0;
    elseif stats.c(j,1) == 1 && stats.c(j,2) == 1  %right, right
        stats.r(j,1) = payoff(2,2);  
        stats.r(j,2) = -2;
    else   %this should not happen
        error('Error in simPenniesBiased');
    end
end

end