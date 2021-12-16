function stats=simBandit(player1,params,n)
% % simBandit %
%PURPOSE:   Simulate a two-armed bandit task
%AUTHORS:   AC Kwan 180426
%
%INPUT ARGUMENTS
%   player1:    player 1
%       label - the name of the algorithm, the name should correspond to a .m file
%       params - parameters associated with that algorithm
%       frac_opto - if a fraction of the trials have optogenetic perturbation
%       opto - parameters associated with that algorithm in perturbation trials
%   params:     task parameters
%      - p_pairs:   reward probability pairs for left/right side
%      - crit_hit:   switching criterion: number of times picked the high-reward side
%      - crit_geo:     switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
%   n:          number of trials to simulate
%
%OUTPUT ARGUMENTS
%   stats:      simulated choice behavior and latent variables

%% initialize
%rng()
rng('default');

stats = []; %stat for the game
stats.currTrial = 1;    % starts at trial number 1
stats.pl=nan(n,1);      % probability to choose left port
stats.c=nan(n,1);       % choice vector
stats.r=nan(n,1);       % reward vector
stats.ql=nan(n,1);      % action value for left choice
stats.qr=nan(n,1);      % action value for right choice
stats.rpe=nan(n,1);     % reward prediction error vector
stats.opto=nan(n,2);    % optogenetic perturbation?

stats.rewardprob=zeros(n,2);   % reward probabilities on left/right side
stats.rule=zeros(n,1);         % rule vector

stats.playerlabel{1} = player1.label;
stats.playerparams{1} = player1.params;

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

%take the text label for the player, there should be a corresponding Matlab
%function describing how the player will play
simplay1 = str2func(player1.label);

%initialize the task
num_p_pairs=size(params.p_pairs,1);
curr_pset=randperm(num_p_pairs);    %pseudo-random ordering of rules (i.e., sets of reward probabilities)
curr_psetIdx=1;  %start from index of 1, each block has unique set of reward probs
stats.rewardprob(1,:)=params.p_pairs(curr_pset(curr_psetIdx),:);
stats.rule(1)=curr_pset(curr_psetIdx);

%% simualte the task

hit = 0;    %for switching, count number of times high-reward side is chosen

for j=1:n
    
    %what is the player's probability for choosing left?
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
        % stats=simplay1(stats,stats.playerparams{1},1);
        stats=simplay1(stats,stats.playerparams{1},1);
    end
    
    %player chooses
    if stats.pl(j)>rand(1)
        stats.c(j)=-1;
    else
        stats.c(j)=1;
    end
    
    %what is the outcome
    if stats.c(j)==-1   % choose left
        if stats.rewardprob(j,1)>=stats.rewardprob(j,2)   % left is high reward side
            hit=hit+1;
        end
        stats.r(j)=(stats.rewardprob(j,1)>=rand(1));  %probabilistic reward
    else                % choose right
        if stats.rewardprob(j,2)>=stats.rewardprob(j,1)   % left is high reward side
            hit=hit+1;
        end
        stats.r(j)=(stats.rewardprob(j,2)>=rand(1));  %probabilistic reward
    end
    
    % did player meet criterion and therefore switch to new set of reward probabilities?
    if j<n  %if we are not at the last trial
        if (hit>=params.crit_hit && params.crit_geo>rand(1))   %yes, met criterion
            hit= 0;   % reset hit
            curr_psetIdx=curr_psetIdx+1;
            if curr_psetIdx > num_p_pairs       %if the current pseudo-random set is exahusted
                last_pset=curr_pset;
                curr_pset=randperm(num_p_pairs);
                while curr_pset(1) == last_pset(num_p_pairs)  %stipulate first transition must involve an actual change in reward probabilities
                    curr_pset=randperm(num_p_pairs);
                end
                curr_psetIdx = 1;
            end
            stats.rewardprob(j+1,:)=params.p_pairs(curr_pset(curr_psetIdx),:);
            stats.rule(j+1)=curr_pset(curr_psetIdx);
        else
            stats.rewardprob(j+1,:)=stats.rewardprob(j,:);
            stats.rule(j+1)=stats.rule(j);
        end
    end
end

end