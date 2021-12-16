function [running_entropy] = cal_runningEntropy(trialData, runningWindow)

%% calculate the running entropy to determine whether the animal is tired in the end of the session
% window: 30 trials

running_entropy = [];
for ii = 1:length(trialData.response)-runningWindow-2
    entropyCount = zeros(1, 8);
    ent = 0;
    for jj=ii:ii+runningWindow-1
        rawchoice = trialData.response(jj:jj+2); % check for miss
        if ~ismember(0, rawchoice)     % if there is miss in the three choices, don't count
           
            choiceSeq=(trialData.response(jj:jj+2)-2)';
        
            entInd = bin2dec(num2str(choiceSeq))+1;
            entropyCount(entInd) = entropyCount(entInd)+1;
        end
    end
    entP = entropyCount / sum(entropyCount);
    if ismember(0,entP)
        warning('certain patterns did not occur, zero probability generated');
        IndexZero=find(~entP);
        for i=1:length(IndexZero)
            entP(IndexZero(i))=realmin;
        end
    end
    for k = 1:8
        ent = ent - entP(k)*log2(entP(k));
    end
    running_entropy = [running_entropy,ent];
end

end
