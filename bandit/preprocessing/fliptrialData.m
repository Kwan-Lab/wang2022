function newtrialData = fliptrialData(trialData)
% % fliptrialData %
%PURPOSE:   Flip the choice direction in the trialData. Why? Animals with
%           left lesion has right side as contralateral side. Animals with
%           right lesion has left side as contralateral side. Will need to
%           flip choice direction for one group to collate the results to
%           look at effect of lesion on sides with respect to lesion.
%AUTHORS:   H Atilgan and AC Kwan 191208
%
%INPUT ARGUMENTS
%   trialData:     the trialData structure
%
%OUTPUT ARGUMENTS
%   newtrialData:  trialData, with the choice direction flipped
%

%% copy directly except for a few entries that require flipping
newtrialData = trialData;
newtrialData = rmfield(newtrialData,'outcome');
newtrialData = rmfield(newtrialData,'rule');
newtrialData = rmfield(newtrialData,'response');
newtrialData = rmfield(newtrialData,'leftlickTimes');
newtrialData = rmfield(newtrialData,'rightlickTimes');

%% for lick times, simply swap

newtrialData.leftlickTimes = trialData.rightlickTimes;
newtrialData.rightlickTimes = trialData.leftlickTimes;

%% for events, need to be more careful depending on the event type

[ STIM, RESP, OUTCOME, RULE, EVENT ] = value_getPresentationCodes(trialData.presCodeSet);

% flip left and right responses
temp_resp = trialData.response;
temp_resp(trialData.response == RESP.LEFT) = RESP.RIGHT;
temp_resp(trialData.response == RESP.RIGHT) = RESP.LEFT;
newtrialData.response = temp_resp;

% flip left and right outcomes
temp_outcome = trialData.outcome;
temp_outcome(trialData.outcome == OUTCOME.REWARDLEFT) = OUTCOME.REWARDRIGHT;
temp_outcome(trialData.outcome == OUTCOME.REWARDRIGHT) = OUTCOME.REWARDLEFT;
temp_outcome(trialData.outcome == OUTCOME.NOREWARDLEFT) = OUTCOME.NOREWARDRIGHT;
temp_outcome(trialData.outcome == OUTCOME.NOREWARDRIGHT) = OUTCOME.NOREWARDLEFT;
newtrialData.outcome = temp_outcome;

% flip left and right reward probabilities (rules)
if trialData.presCodeSet == 3 || trialData.presCodeSet == 8 
    temp_rule = trialData.rule;
    temp_rule(trialData.rule == RULE.L70R10) = RULE.L10R70;
    temp_rule(trialData.rule == RULE.L10R70) = RULE.L70R10;
    newtrialData.rule = temp_rule;
else
    error('Error in fliptrialData: Currently the code does not support task beyond reversal');
end

end
