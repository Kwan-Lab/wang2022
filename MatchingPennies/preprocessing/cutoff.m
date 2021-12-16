function [cutoff, segline] = cutoff(runningEntropy)

% determine where to cut in order to get rid of the fatigue effect

% criterion: 
% 1) cut-off point should be larger than 250
% 2) the average of entropy is lower than 1
% 3) the average of entropy after the cut-off point never get back to
% original level 

try
    [TF, S1, S2] = ischange(runningEntropy, 'linear', 'MaxNumChanges',2);
catch
    warning('There is no result!');
    TF = NaN;
    S1 = NaN;
    S2 = NaN;
end

if ~isnan(TF)
    segline = S1.*(1:length(TF)) + S2;
    
    
    xAxis = 1:length(segline);
    
    % find the point where the average entropy go below 1
    
    IndAbove  =segline<1;  % critierion 2
    TF_one  = xAxis(ischange(double(IndAbove))); % that's include the change point both go below 1 and go above 1
    TF_belowone = [];
    for k = 1:length(TF_one)
        if segline(TF_one(k)-1) >1
            TF_belowone = [TF_belowone, TF_one(k)];
        end
    end
    
    tfPoint = xAxis(TF);
    cutoff = 0;
    
    % if isnan(TF_belowone)  % if the entropy go below 1 gradulally
    %     % these seems to be redundant? but I'll keep them for now
    %      for ii = 1:length(tfPoint)
    %         if tfPoint(ii) > 300 & segline(tfPoint(ii)) < 1% criterion 1
    %             if max(segline(tfPoint(ii)+1:end)) < segline(tfPoint(ii)-1)-0.1 % criterion 3
    %                 cutoff = tfPoint(ii);
    %                 return;
    %             end
    %         end
    %     end
    
    %else % if there is a sudden drop
    for ii = 1:length(TF_belowone)
        % check the point 1 by 1
        if TF_belowone(ii) > 250  % critierion 1
            % finding the change point before I(ii)
            if ismember(TF_belowone(ii), tfPoint)  % the point itself is  a change point
                if segline(TF_belowone(ii)) < 1% criterion 2
                    if max(segline(TF_belowone(ii)+1:end)) < segline(TF_belowone(ii)-1) % criterion 3
                        cutoff = TF_belowone(ii);
                        return;
                    end
                end
            else
                tf = find(tfPoint<TF_belowone(ii), 1,'last');
                
                tfP = tfPoint(tf);
                if max(segline(TF_belowone(ii)+1:end)) < max(segline(tfP-1),segline(tfP))  % criterion 3
                    cutoff = TF_belowone(ii);
                    return;
                end
            end
        end
    end
else
    cutoff = 0;
    segline = zeros(1, length(runningEntropy));
end  
%end

end
