function [stats_all] = MP_merge_sessions(dataIndex)
% % MP_merge_sessions %
%PURPOSE:   Merge different sessions into one long session
%
%INPUT ARGUMENTS
%   dataIndex:    a database index table for the sessions to analyze
%
%OUTPUT ARGUMENTS
%   stats_all:    concatenated choices and outcomes



stats_all.c = [];
stats_all.r = [];
stats_all.sessionLength = [];
for i=1:size(dataIndex,1)
    
    load(fullfile(dataIndex.BehPath{i},[dataIndex.LogFileName{i}(1:end-4),'_beh_cut.mat']));

    lick_trType_array{i}=lick_trType;
    
    lregRCUC_array{i}=lregRCUC_output;
    lregCRInt_array{i}=lregCRInt_output;
    
    iti_array{i}=iti_trType;
    respTime_array{i}=respTime_trType;
    trueRespTime{i} = trialData.rt;
    choiceBySession{i} = stats;
    
    if exist('nlike_array','var')
        fname = fieldnames(nlike);
        for j=1:numel(fname)  %append for each field
            nlike_array.(fname{j})=[nlike_array.(fname{j}); nlike.(fname{j})];
            bic_array.(fname{j})=[bic_array.(fname{j}); bic.(fname{j})];
        end
        fname = fieldnames(fitpar);
        for j=1:numel(fname)  %append for each field
            fitpar_array.(fname{j})=[fitpar_array.(fname{j}); fitpar.(fname{j})];
        end
    else
        nlike_array = nlike;
        bic_array = bic;
        fitpar_array = fitpar;
    end
    
    nTrial_array(i)=sum(stats.c(:,1)==-1)+sum(stats.c(:,1)==1);
    entro_array(i)=entro;
    rrate_array(i)=sum(stats.r==1)/(sum(stats.r==1)+sum(stats.r==0));
    
    stats_all.c=[stats_all.c; stats.c];
    stats_all.r=[stats_all.r; stats.r];
    stats_all.sessionLength = [stats_all.sessionLength,length(stats.c)];
    close all;
    clearvars -except i dirs dataIndex ...
        lick_trType_array lregRCUC_array lregCRInt_array iti_array respTime_array trueRespTime choiceBySession...
        nlike_array bic_array fitpar_array...
        nTrial_array entro_array rrate_array subMask...
        stats_all;
    
%     trials = MP_getTrialMasks(trialData);
%         
%     if i==1
%         trialData.aveEntropy = [trialData.aveEntropy, zeros(1, length(trialData.cue)-length(trialData.aveEntropy))];
%         trialData.runningEntropy = [trialData.runningEntropy, zeros(1, length(trialData.cue)-length(trialData.runningEntropy))];
%         trialDataCombined = trialData;
%         trialsCombined = trials;
%     else
%         fields=fieldnames(trialData);
%         for j = 1:numel(fields)
%             if ~strcmp(fields{j}, 'trigger') && ~strcmp(fields{j},'triggerTimes')  % no need to concatenate trigger information
%                 if iscell(trialDataCombined.(fields{j}))    %licktimes are stored in cells
%                     trialDataCombined.(fields{j}) = [trialDataCombined.(fields{j}), cell(1, n_nan), trialData.(fields{j})];
%                 else
%                     if strcmp(fields{j},'aveEntropy') || strcmp(fields{j},'runningEntropy')
%                         trialData.(fields{j})=[trialData.(fields{j}),zeros(1, length(trialData.cue)-length(trialData.(fields{j})))];
%                         trialDataCombined.(fields{j}) = [trialDataCombined.(fields{j}), nan(1,n_nan), trialData.(fields{j})];
%                     else
%                         trialDataCombined.(fields{j}) = [trialDataCombined.(fields{j}); nan(n_nan,1); trialData.(fields{j})];
%                     end
%                 end
%             end
%         end
%         
%         fields=fieldnames(trials);
%         for j = 1:numel(fields)
%             trialsCombined.(fields{j}) = [trialsCombined.(fields{j}); nan(n_nan,1); trials.(fields{j})];
%         end
%         
%     end
     
end

end


