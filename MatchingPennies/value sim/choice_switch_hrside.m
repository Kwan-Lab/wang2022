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

%which side has the high reward probability?
hr_side=nan(size(stats.c));   
hr_side(stats.rewardprob(:,1)>stats.rewardprob(:,2))=-1; %high-reward side is left
hr_side(stats.rewardprob(:,1)<stats.rewardprob(:,2))=1;  %high-reward side is right

%run-length encoding to find when switches occur
x0=hr_side';
l = diff([ 0 find(x0(1:end-1) ~= x0(2:end)) length(x0) ]);
v = x0([ find(x0(1:end-1) ~= x0(2:end)) length(x0) ]);

choseh=zeros(1+2*trials_back,1); %prob of choosing high rew side
chosel=zeros(1+2*trials_back,1); %prob of choosing low rew side
numtrans=0;

for i=1:numel(l)-1    
    trial1=(sum(l(1:i))+1)-trials_back;
    trial2=(sum(l(1:i))+1)+trials_back;
    if trial2<=numel(stats.c)
        numtrans=numtrans+1;
        choseh=choseh+(stats.c(trial1:trial2) == hr_side(trial1));
        chosel=chosel+(stats.c(trial1:trial2) == -1*hr_side(trial1));    
    end
end

probh=choseh/numtrans;
probl=chosel/numtrans;

n=[-trials_back:1:trials_back]';

output.n=n;
output.probh=probh;
output.probl=probl;

end


    