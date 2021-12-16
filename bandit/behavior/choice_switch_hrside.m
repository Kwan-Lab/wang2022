function output=choice_switch_hrside(stats,trials_back)
% % choice_switch_hrside %
%PURPOSE:   Analyze choice behavior around switches, when the high reward
%           side has changed
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   stats:  stats of the task
%   trials_back: how many trials around switch to plot
%
%OUTPUT ARGUMENTS
%   output:     numbers used to plot figure
%
% To plot the output, use plot_switch_hrside().

%%

choseh=zeros(1+2*trials_back,1); %prob of choosing high rew side
chosel=zeros(1+2*trials_back,1); %prob of choosing low rew side
choseneither=zeros(1+2*trials_back,1); %prob of choosing neither (miss)

numtrans=0;

for i=1:numel(stats.blockLength)-1
    
    %which kind of transition is it?
    idx=stats.blockTrans(i);
    
    if ~isnan(idx)   %NaN could arise if we are using merge_session
        
        trial1=(sum(stats.blockLength(1:i))+1)-trials_back;
        trial2=(sum(stats.blockLength(1:i))+1)+trials_back;
        if trial2<=numel(stats.c)
            numtrans=numtrans+1;
            choseh=choseh+(stats.c(trial1:trial2) == stats.hr_side(trial1));
            chosel=chosel+(stats.c(trial1:trial2) == -1*stats.hr_side(trial1));
            choseneither=choseneither+isnan(stats.c(trial1:trial2));
        end
        
    end
end

probh=choseh/numtrans;
probl=chosel/numtrans;
probneither=choseneither/numtrans;

n=[-trials_back:1:trials_back]';

output.n=n;
output.probh=probh;
output.probl=probl;
output.probneither=probneither;

end


    