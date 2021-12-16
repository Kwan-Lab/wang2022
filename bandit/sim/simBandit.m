function stats=simBandit(player1,params,n)
% % simBandit %
%PURPOSE:   Simulate a two-armed bandit task
%AUTHORS:   AC Kwan 180426
%
%INPUT ARGUMENTS
%   player1:    player 1
%       label - the name of the algorithm, the name should correspond to a .m file
%       params - parameters associated with that algorithm
%   params:     task parameters
%      - p_pairs:   reward probability pairs for left/right side
%      - crit_hit:   switching criterion: number of times picked the high-reward side
%      - crit_geo:     switching criterion: after hit, random number of trials defined by geometric dist. p=1/mean
%   n:          number of trials to simulate
%
%OUTPUT ARGUMENTS
%   stats:      simulated choice behavior and latent variables

%% initialize
rng('default');
rng('shuffle');

stats = []; %stat for the game
stats.currTrial = 1;    % starts at trial number 1
stats.pl=nan(n,1);      % probability to choose left port
stats.c=nan(n,1);       % choice vector
stats.r=nan(n,1);       % reward vector
stats.ql=nan(n,1);      % action value for left choice
stats.qr=nan(n,1);      % action value for right choice
stats.rpe=nan(n,1);     % reward prediction error vector

stats.rewardprob=zeros(n,2);   % reward probabilities on left/right side
stats.rule=zeros(n,1);         % rule vector
stats.hr_side=nan(n,1);      % high-probability side vector

stats.playerlabel{1} = ['algo_' player1.model_type];
stats.playerparams{1} = player1.params;

stats.rule_labels = params.rule_labels;

%take the text label for the player, there should be a corresponding Matlab
%function describing how the player will play
simplay1 = str2func(['algo_' player1.model_type]);

%initialize the task
nRules=size(params.p_pairs,1);
curr_pset=randperm(nRules);    %pseudo-random ordering of rules (i.e., sets of reward probabilities)
curr_psetIdx=1;                %start from index of 1, each block has unique set of reward probs
stats.rewardprob(1,:)=params.p_pairs(curr_pset(curr_psetIdx),:);
stats.rule(1)=curr_pset(curr_psetIdx);

%% simualte the task

hit = 0;    %for switching, count number of times high-reward side is chosen

for j=1:n
    
    %what is the player's probability for choosing left?
    stats.currTrial = j;

    %which is the high-probability side?
    if stats.rewardprob(j,1) > stats.rewardprob(j,2)
        stats.hr_side(j) = -1;
    elseif stats.rewardprob(j,1) < stats.rewardprob(j,2)
        stats.hr_side(j) = 1;
    else
        stats.hr_side(j) = 0;
    end
    
    stats=simplay1(stats,stats.playerparams{1});
    
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
            if curr_psetIdx > nRules       %if the current pseudo-random set is exahusted
                last_pset=curr_pset;
                curr_pset=randperm(nRules);
                while curr_pset(1) == last_pset(nRules)  %stipulate first transition must involve an actual change in reward probabilities
                    curr_pset=randperm(nRules);
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

stats = value_getTrialStatsMore(stats);

end