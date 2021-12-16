function [ stats ] = value_getTrialStatsMore( stats_input )
% % value_getTrialStatsMore %
%PURPOSE:   Process 'stats to add more fields to data structure: 'stats' 
%           containing dummy codes for task variables
%AUTHORS:   H Atilgan and AC Kwan 191202
%
%INPUT ARGUMENTS
%   stats_iput: Structure containing fields, each dummy-coded 
%               tabulating the animal's choice 
%               (e.g., trial.c = -1 for left, 1 for right)
%
%OUTPUT VARIABLES
%   stats:      Structure containing fields, each dummy-coded 
%               tabulating the animal's choice, but more!
%               (e.g., trial.c = -1 for left, 1 for right)

stats = stats_input;

nTrials = size(stats.c,1);
nRules = numel(stats.rule_labels);

%% these are all stats of the session that can be derived from other stats

%which side has the high reward probability?
stats.hr_side=nan(size(stats.c));
stats.hr_side(stats.rewardprob(:,1)>stats.rewardprob(:,2))=-1; %high-reward probability side is left
stats.hr_side(stats.rewardprob(:,1)<stats.rewardprob(:,2))=1;  %high-reward probability side is right

%% Tabulate block-by-block stats

%run-length encoding to find when switches occur
x0=stats.rule';
stats.blockLength = diff([ 0 find(x0(1:end-1) ~= x0(2:end)) length(x0) ])';  %total length of each block
stats.blockRule = x0([ find(x0(1:end-1) ~= x0(2:end)) length(x0) ])';        %reward probabilities / rule associated with each block

%find all the possible transitions between sets of reward probabilities
stats.ruletransList=[];
for i=1:nRules
    for j=1:nRules
        if i~=j
            stats.ruletransList=[stats.ruletransList; i j];    %all the possible rule transitions
        end
    end
end

%go through each block, identify transition, trial to performance, and random number added
stats.blockTrans = nan(numel(stats.blockLength)-1,1);        %type of transition 
stats.blockTrialtoCrit = nan(numel(stats.blockLength)-1,1);  %trial to criterion  
stats.blockTrialRandomAdded = nan(numel(stats.blockLength)-1,1);  %random number added
stats.blockPreSwitchBetterChoiceRate = nan(numel(stats.blockLength)-1,1); %rate of choosing initial better option in the 5 trials prior to switch
stats.blockPreSwitchWorseChoiceRate = nan(numel(stats.blockLength)-1,1); %rate of choosing initial worse option in the 5 trials prior to switch
stats.blockPreSwitchBetterChoiceAtSwitch = nan(numel(stats.blockLength)-1,1); %rate of choosing initial better option for the trial prior to switch
stats.blockPreSwitchWorseChoiceAtSwitch = nan(numel(stats.blockLength)-1,1); %rate of choosing initial worse option for the trial prior to switch

for i=1:numel(stats.blockLength)-1
    %which kind of transition is it?
    %this will exclude the rule==NaN if combining trials using merge_sessions function
    idx=find(stats.ruletransList(:,1)==stats.blockRule(i) & stats.ruletransList(:,2)==stats.blockRule(i+1));
    
    if ~isempty(idx)
        stats.blockTrans(i) = idx;
        
        %find the trial in which animal's choice fulfilled the criterion - choosing high-probability side for 10 times
        thTrial = find(cumsum(stats.c((sum(stats.blockLength(1:i-1))+1):(sum(stats.blockLength(1:i))))==stats.hr_side(sum(stats.blockLength(1:i-1))+1))==10); %trial meeting switching criterion
        if isempty(thTrial) % this should not happen
            stats.blockTrialtoCrit(i) = inf;        %never met criterion
        else
            stats.blockTrialtoCrit(i) = thTrial(1); %first trial when criterion is met
        end
        
        %random number added is block length minus the trial to meet criterion
        stats.blockTrialRandomAdded(i) = stats.blockLength(i) - stats.blockTrialtoCrit(i);
        
        %rate of choosing initial better or worse option in the 5 trials prior to switch
        %--- averaging 5 trials, for cases when there are fewer switches, e.g. lesion data
        stats.blockPreSwitchBetterChoiceRate(i) = nanmean((stats.c(sum(stats.blockLength(1:i))-4:sum(stats.blockLength(1:i)))==stats.hr_side(sum(stats.blockLength(1:i-1))+1)));
        stats.blockPreSwitchWorseChoiceRate(i) = nanmean((stats.c(sum(stats.blockLength(1:i))-4:sum(stats.blockLength(1:i)))==-1*stats.hr_side(sum(stats.blockLength(1:i-1))+1)));

        %rate of choosing initial better or worse option for the trial prior to switch
        %--- looking specifically at the trial prior to switch, avoiding
        %contamination of an estimate that include trials used to reach
        %performance criterion which by definitely are more likely to be
        %correct. However this requires more data.
        stats.blockPreSwitchBetterChoiceAtSwitch(i) = (stats.c(sum(stats.blockLength(1:i)))==stats.hr_side(sum(stats.blockLength(1:i-1))+1));
        stats.blockPreSwitchWorseChoiceAtSwitch(i) = (stats.c(sum(stats.blockLength(1:i)))==-1*stats.hr_side(sum(stats.blockLength(1:i-1))+1));
        
    end
end

%% check consistency
if any(stats.blockTrialtoCrit<10)        %should need at least 10 trials to achieve choosing criterion of choosing high-prob side for 10 times
    disp('Error in value_getTrialStats: fewer than 10 trials to reach criterion for block switch ?!');
end    
if any(stats.blockTrialRandomAdded<0)  %random number added should never be fewer than 0
    disp('Error in value_getTrialStats: random number added for block length < 0 ?!');
end    

end

