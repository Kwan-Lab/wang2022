function output=choice_switch_random(stats,trials_back,L1_ranges,L2_ranges)
% % choice_switch_random %
%PURPOSE:   Analyze choice behavior around switches, separate switches
%           depending on the random number added for the block preceding
%           the switch
%AUTHORS:   H Atilgan and AC Kwan 191205
%
%INPUT ARGUMENTS
%   stats:  stats of the task
%   trials_back: how many trials around switch to plot
%   L1_ranges: consider only subset of blocks within the range, for trials to 
%                  reach criterion for the block preceding the switch (ranges
%                  are *inclusive*)
%   L2_ranges: consider only subset of blocks within the range, for random
%                  number added for the block preceding the switch (ranges
%                  are *inclusive*)
%
%OUTPUT ARGUMENTS
%   output:     numbers used to plot figure
%
% To plot the output, use plot_switch_random().

%%

if size(L1_ranges,1) ~= size(L2_ranges,1)
    error('Error in choice_switch_random: The L1_range and L2_range should have the same number of rows');
else
    numRange = size(L1_ranges,1);
end

%how many different types of rule transitions?
numruletrans=size(stats.ruletransList,1);

if numruletrans==0 % No switch occured
    disp('No switch of reward probabilities in this session');
    output = 0;
else
    chosel=zeros(1+2*trials_back,numruletrans,numRange); %prob of choosing left
    choser=zeros(1+2*trials_back,numruletrans,numRange); %prob of choosing right
    choseneither=zeros(1+2*trials_back,numruletrans,numRange); %prob of not choosing (miss)
    
    numtrans=zeros(numruletrans,numRange);    %number of that type of transition
    
    for i=1:numel(stats.blockLength)-1
        %which kind of transition is it?
        idx=stats.blockTrans(i);
        
        if ~isnan(idx)   %NaN could arise if we are using merge_session
            
            %within which range does this block fall into?
            range_idx = (stats.blockTrialtoCrit(i) >= L1_ranges(:,1)) & (stats.blockTrialtoCrit(i) <= L1_ranges(:,2)) &...
                (stats.blockTrialRandomAdded(i) >= L2_ranges(:,1)) & (stats.blockTrialRandomAdded(i) <= L2_ranges(:,2));
            
            if sum(range_idx) == 1   %if this transition falls into one of the subset of switches to be considered
                %what were the choices around that transition
                %note: output.n=0 is the first trial with the switched probabilities
                trial1=(sum(stats.blockLength(1:i))+1)-trials_back;  %trials before switch
                trial2=(sum(stats.blockLength(1:i))+1)+trials_back;  %trials after switch
                
                if trial2<=numel(stats.c)   %if trial after switch does not exceed end of session
                    numtrans(idx,range_idx)=numtrans(idx,range_idx)+1;
                    chosel(:,idx,range_idx)=chosel(:,idx,range_idx)+(stats.c(trial1:trial2)==-1);
                    choser(:,idx,range_idx)=choser(:,idx,range_idx)+(stats.c(trial1:trial2)==1);
                    choseneither(:,idx,range_idx)=choseneither(:,idx,range_idx)+isnan(stats.c(trial1:trial2));
                end
            end
        end
    end
    
    probl=nan(1+2*trials_back,numruletrans,numRange); %prob of choosing left
    probr=nan(1+2*trials_back,numruletrans,numRange); %prob of choosing right
    probneither=nan(1+2*trials_back,numruletrans,numRange); %prob of choosing neither
    for j=1:numruletrans
        for k=1:numRange
            if numtrans(j,k)>0
                probl(:,j,k)=chosel(:,j,k)/numtrans(j,k);
                probr(:,j,k)=choser(:,j,k)/numtrans(j,k);
                probneither(:,j,k)=choseneither(:,j,k)/numtrans(j,k);
            end
        end
    end
    
    n=[-trials_back:1:trials_back]';
    
    output.n=n;
    output.probl=probl;
    output.probr=probr;
    output.probneither=probneither;
    output.numtransType=numruletrans;
    output.transType=stats.ruletransList;
    output.numRange=numRange;
    output.L1_ranges=L1_ranges;
    output.L2_ranges=L2_ranges;
    
    output.L1=stats.blockTrialtoCrit(~isnan(stats.blockTrans));
    output.L2=stats.blockTrialRandomAdded(~isnan(stats.blockTrans));
    
end

end
