function MP_determineBehCriteria(dataIndex)

%% there is (possible) fatigue effects in some of the sessions in MP training
% this code is designed to find the right point where the animal stops to
% play the game by choosing only one side


% criterion: 
% 1) cut-off point should be larger than 250
% 2) the average of entropy is lower than 1
% 3) the average of entropy after the cut-off point never get back to
% original level (difference larger than 0.1?)

nFiles = size(dataIndex,1);
for ii = 1:nFiles
    fullfilePath = [dataIndex.BehPath{ii}, '\', dataIndex.LogFileName{ii}(1:end-4),'_beh.mat'];
    load(fullfilePath);
    
    trials = MP_getTrialMasks(trialData);
    % calculate the running entropy
    running_window = 30;
    running_entropy = cal_runningEntropy(trialData, running_window); 
    
    % get the cut point
    [cutPoint, ave_entropy] = cutoff(running_entropy);
    
    % save the information in .mat file
    trialData.cutPoint = cutPoint;
    trialData.aveEntropy = ave_entropy;
    trialData.runningEntropy = running_entropy;
    % save the .mat file
    save(fullfilePath,...
    'trialData','sessionData','logfileData');
% );,'trials',...
%             'lregRCUC_output','lregCRInt_output',...
%             'lick_trType','iti_trType','respTime_trType',...
%             'fitpar','bic','nlike',...
%             'entro',...
%             'stats');
    % save the figure for later check
    saveFigPath = [dataIndex.BehPath{ii}, '\cut\'];
    if ~exist(saveFigPath)
        mkdir(saveFigPath)
    end
    plot_cutoff(trialData, saveFigPath, dataIndex(ii,:));
end