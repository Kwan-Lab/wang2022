function [includeIndex, excludeIndex] = determineBehCriteria(dataIndex)
% % determineBehCriteria %
%PURPOSE:   Check each session to see if it fulfills certain performance
%           criteria to be included for subsequent analyses
%AUTHORS:   H Atilgan and AC Kwan 191203
%
%INPUT ARGUMENTS
%   dataIndex:  a table of the data files
%
%OUTPUT ARGUMENTS
%   includeIndex:  a table of the data files that fulfill criteria
%   excludeIndex:  a table of the data files that will be excluded
%

%% The criteria

% Criterion 1: number of responsive trials (chose left or right) should exceed 100
numTrial = 100;

% Criterion 2: number of switches should exceed 2
numSwitch = 2;

%% Create criteria-related database index table

nFile = size(dataIndex,1);

critIndex = table(...
    NaN(nFile,1),...
    NaN(nFile,1),...
    NaN(nFile,1)...
    );

critIndex.Properties.VariableNames = {...
    'respNum',...   number of responsive trials in the session
    'switchNum',...  number of switches in the session
    'motorBias'...  absolute difference in response time for left versus right
    };

%% Calculate session performance relevant to the criteria

disp(['-----------------------------------------------------------']);
disp(['--- Determine behavioral criteria for ' int2str(size(dataIndex,1)) ' behavioral logfiles.']);
disp(['-----------------------------------------------------------']);

for i = 1:size(dataIndex,1)
    if ~isnan(dataIndex.BehCreated(i))
        
        disp(['Processing file #' int2str(i) '.']);
        
        load(fullfile(dataIndex.BehPath{i},['bandit_',dataIndex.LogFileName{i}(end-29:end-4),'_beh.mat']));
        
        [ trials ] = value_getTrialMasks( trialData );
        
        % Calculate number of responsive trials in a session
        critIndex.trialNum(i) = sum(trials.left) + sum(trials.right);
        
        % Calculate number of switches in a session
        x0=trialData.rule';
        switchTrials = find(x0(1:end-1) ~= x0(2:end)) + 1;  %trial numbers, for very first trial after a switch
        critIndex.switchNum(i) = numel(switchTrials);  %number of switches
        
        % Calculate motor bias in terms of left vs right response times
        trialType={'go','left','right'};
        edges=[-1:0.05:2];
        valLabel='Response time (s)';
        respTime_trType=get_val_byTrialType_bandit(trialData.rt,trials,trialType,edges,valLabel);
        critIndex.motorBias(i) = abs(respTime_trType.valMedian{2} - respTime_trType.valMedian{3});
        
    else
        error(['Error in determineBehCriteria: Behavioral .mat file was not created for file #' int2str(i)]);
    end
end

combinedIndex = [dataIndex critIndex];

%% Determine which session meets performance criteria

crit1 = combinedIndex.trialNum > numTrial;
crit2 = combinedIndex.switchNum > numSwitch;
critMet = all([crit1 crit2],2);

disp(['Out of ' int2str(size(dataIndex,1)) ' sessions, ' int2str(sum(critMet)) ' met performance criteria.']);
includeIndex = combinedIndex(critMet,:);
excludeIndex = combinedIndex(~critMet,:);

end

