function output=bandit_get_lickrate_byTrialType(trialData,trials,trialType,edges)
% % get_lickrate_byTrialType %
%PURPOSE:   Analyze lick rate for different trial types
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   trialData:  Structure generated by value_getSessionData().
%   trials:     Structure generated by value_getTrialMasks().
%   trialType:  The trial types (fieldnames of the trials variable)
%   edges:      The edges used for a histc() operation to get lick rate histogram
%
% To plot the output, use plot_lickrate_byTrialType().

%% 

leftTimes=[]; 
rightTimes=[]; 
edgeWidth=nanmean(diff(edges));

for i=1:numel(trialType)
    
    %trials belonging to specific response or outcome or rule type
    if numel(trialType{i})==1       %if only one event specified
        tempMask = (trials.(trialType{i})==1);
        trialLabel{i} = trialType{i};
    elseif numel(trialType{i})==2   %a conjunction of two events
        tempMask = (trials.(trialType{i}{1})==1) & (trials.(trialType{i}{2})==1);
        trialLabel{i} = [trialType{i}{1} ' + ' trialType{i}{2}];
    else
        error('Error in get_lickrate_byTrialType: Currently do not support conjunction of more than two events');
    end
        
    %histogram of all left lick times = left lick rate
    temp=histc([trialData.leftlickTimes{tempMask}],edges)/sum(tempMask)/edgeWidth;   %in Hz
    if ~isempty(temp)    %if there are any such trials and licks
        leftTimes{i}=temp(1:end-1)';
    else                 %otherwise fill with NaN
        leftTimes{i}=nan(size(edges(1:end-1)'));
    end
    
    %histogram of all right lick times = right lick rate    
    temp=histc([trialData.rightlickTimes{tempMask}],edges)/sum(tempMask)/edgeWidth;   %in Hz
    if ~isempty(temp)    
        rightTimes{i}=temp(1:end-1)';
    else
        rightTimes{i}=nan(size(edges(1:end-1)'));
    end
    
end

output.trialType=trialType;
output.trialLabel=trialLabel;
output.edges=edges;
output.leftTimes=leftTimes;
output.rightTimes=rightTimes;

end


