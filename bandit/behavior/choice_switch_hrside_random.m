function output=choice_switch_hrside_random(stats,trials_back,L1_ranges,L2_ranges)
% % choice_switch_hrside_random %
%PURPOSE:   Analyze choice behavior around switches, when the high reward
%           side has changed, separate switches depending on the random
%           number added for the block preceding
%           the switch
%AUTHORS:   H Atilgan and AC Kwan 191206
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
% To plot the output, use plot_switch_hrside_random().

%%

if size(L1_ranges,1) ~= size(L2_ranges,1)
    error('Error in choice_switch_random: The L1_range and L2_range should have the same number of rows');
else
    numRange = size(L1_ranges,1);
end

choseh=zeros(1+2*trials_back,numRange); %prob of choosing high rew side
chosel=zeros(1+2*trials_back,numRange); %prob of choosing low rew side
choseneither=zeros(1+2*trials_back,numRange); %prob of choosing neither (miss)

numtrans=zeros(1,numRange);

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
            trial1=(sum(stats.blockLength(1:i))+1)-trials_back;
            trial2=(sum(stats.blockLength(1:i))+1)+trials_back;
            if trial2<=numel(stats.c)
                numtrans(range_idx)=numtrans(range_idx)+1;
                choseh(:,range_idx)=choseh(:,range_idx)+(stats.c(trial1:trial2) == stats.hr_side(trial1));
                chosel(:,range_idx)=chosel(:,range_idx)+(stats.c(trial1:trial2) == -1*stats.hr_side(trial1));
                choseneither(:,range_idx)=choseneither(:,range_idx)+isnan(stats.c(trial1:trial2));
            end
        end
    end
end

probh=nan(1+2*trials_back,numRange); %prob of choosing initial best option
probl=nan(1+2*trials_back,numRange); %prob of choosing initial worse option
probneither=nan(1+2*trials_back,numRange); %prob of missing
for j=1:numRange
    probh(:,j)=choseh(:,j)/numtrans(j);
    probl(:,j)=chosel(:,j)/numtrans(j);
    probneither(:,j)=choseneither(:,j)/numtrans(j);
end

n=[-trials_back:1:trials_back]';

output.n=n;
output.probh=probh;
output.probl=probl;
output.probneither=probneither;

output.numRange=numRange;
output.L1_ranges=L1_ranges;
output.L2_ranges=L2_ranges;

output.L1=stats.blockTrialtoCrit(~isnan(stats.blockTrans));
output.L2=stats.blockTrialRandomAdded(~isnan(stats.blockTrans));

end


