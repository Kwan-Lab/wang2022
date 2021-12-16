function [ stats ] = value_getTrialStats( trials, nRules )
% % value_getTrialStats %
%PURPOSE:   Create data structures: 'stats' containing dummy codes for task 
%           variables.
%AUTHORS:   H Atilgan and AC Kwan 191202
%
%INPUT ARGUMENTS
%   trials:     Structure containing fields, each a logical mask
%               indicating whether animal made a left choice 
%               (e.g., trials.left)
%   nRules:     Number of sets of reward probabilities
%
%OUTPUT VARIABLES
%   stats:      Structure containing fields, each dummy-coded 
%               tabulating the animal's choice 
%               (e.g., trial.c = -1 for left, 1 for right)

nTrials = numel(trials.go);

%% Tabulate trial-by-trial stats (dummy codes)

%choice: left=-1; right=1; miss=NaN
stats.c=nan(nTrials,1);
stats.c(trials.left==1)=-1;
stats.c(trials.right==1)=1;

%outcome: reward=1; no reward:0; miss=NaN
stats.r=nan(nTrials,1);
stats.r(trials.reward==1)=1;
stats.r(trials.noreward==1)=0;

%rule: the type of reward probabilities
stats.rule=nan(nTrials,1);
if nRules == 2          %two sets of reward prob
    stats.rule(trials.L70R10==1)=1; 
    stats.rule(trials.L10R70==1)=2;
    stats.rule_labels = {'0.7:0.1','0.1:0.7'};
    probList=[0.7 0.1; 0.1 0.7];
elseif nRules == 6      %six sets of reward prob
    stats.rule(trials.L70R30==1)=1; 
    stats.rule(trials.L70R10==1)=2;
    stats.rule(trials.L30R10==1)=3;
    stats.rule(trials.L30R70==1)=4;
    stats.rule(trials.L10R70==1)=5;
    stats.rule(trials.L10R30==1)=6;
    stats.rule_labels = {'0.7:0.3','0.7:0.1','0.3:0.1','0.3:0.7','01:0.7','01:0.3'};
    probList=[ 0.7 0.3; 0.7 0.1; 0.3 0.1; 0.3 0.7; 0.1 0.7; 0.1 0.3];
end

%reward probabilities for left and right sides
stats.rewardprob=nan(nTrials,2);
for j = 1:size(probList,1)      %number of reward probabilities sets
    for k = 1:size(probList,2)  %left and right
        stats.rewardprob(stats.rule==j,k) = probList(j,k);
    end
end

end

