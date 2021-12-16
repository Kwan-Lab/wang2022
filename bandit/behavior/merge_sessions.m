function [trialDataCombined, trialsCombined, nRules] = merge_sessions(dataIndex)
% % merge_sessions %
%PURPOSE:   Merge different sessions into one long session, with the gaps
%           between sessions filled with NaNs
%AUTHORS:   H Atilgan and AC Kwan 191204
%
%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%
%OUTPUT ARGUMENTS
%   trialDataCombined:  the concatenated 'trialData' structure
%   trialsCombined:     the concatenated 'trials' structure
%   nRules:             number of sets of reward probabilities

n_nan = 20;   %insert this many NaN in each gap

sessionLength = [];
aveReward = [];

for i=1:size(dataIndex,1)
    
    load(fullfile(dataIndex.BehPath{i}, ['bandit_',dataIndex.LogFileName{i}(end-29:end-4),'_beh.mat']));
    trials = value_getTrialMasks(trialData);
    % calculate number of nolicks
    sessionLength = [sessionLength, length(trials.left)];
    aveReward = [aveReward, mean(trials.reward)];
 
    if i==1
        trialDataCombined = trialData;
        trialsCombined = trials;
        nRules = sessionData.nRules;
    else
        fields=fieldnames(trialData);
        for j = 1:numel(fields)
            if ~strcmp(fields{j}, 'trigger') & ~strcmp(fields{j}, 'triggerTimes')  % no need for trigger
                if iscell(trialDataCombined.(fields{j}))    %licktimes are stored in cells
                    trialDataCombined.(fields{j}) = [trialDataCombined.(fields{j}); cell(n_nan,1); trialData.(fields{j})];
                else
                    trialDataCombined.(fields{j}) = [trialDataCombined.(fields{j}); nan(n_nan,1); trialData.(fields{j})];
                end
            end
        end
        
        fields=fieldnames(trials);
        for j = 1:numel(fields)
            trialsCombined.(fields{j}) = [trialsCombined.(fields{j}); nan(n_nan,1); trials.(fields{j})];
        end
        
        if nRules ~= sessionData.nRules
            error('Error in merge_sessions: the nRules for the sessions do not match');
        end
    end
     
end
trialDataCombined.sessionLength = sessionLength;
trialDataCombined.aveReward = aveReward;

end


