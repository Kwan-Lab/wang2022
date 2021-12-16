function output=choice_block(stats,trials_forw)
% % choice_block %
%PURPOSE:   Analyze choice behavior along a block
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   stats:  stats of the task
%   trials_forw: how many trials into the block to plot
%
%OUTPUT ARGUMENTS
%   output:     numbers used to plot figure
%
% To plot the output, use plot_block().

%%
%run-length encoding to find when switches occur
x0=stats.rule';
l = diff([ 0 find(x0(1:end-1) ~= x0(2:end)) length(x0) ]);
v = x0([ find(x0(1:end-1) ~= x0(2:end)) length(x0) ]);

%number of reward probability pairs
num_probpairs=size(unique(stats.rewardprob,'rows'),1);

%find the unique reward probabilities
blockType=1:num_probpairs;
numblockType=num_probpairs;

chosel=zeros(1+trials_forw,numblockType); %prob of choosing left
choser=zeros(1+trials_forw,numblockType); %prob of choosing right
rewarded=zeros(1+trials_forw,numblockType); %prob of rewarded
numtrans=zeros(1,numblockType);    %number of that type of transition

for i=1:numel(l)  
    %which kind of transition is it?
    idx=find(blockType==v(i));
    
    %what were the choices around that transition
    if i==1
        trial1=1;
    else
        trial1=(sum(l(1:i-1))+1);
    end
    trial2=trial1+trials_forw;
    if trial2<=numel(stats.c)
        numtrans(idx)=numtrans(idx)+1;
        chosel(:,idx)=chosel(:,idx)+(stats.c(trial1:trial2)==-1);
        choser(:,idx)=choser(:,idx)+(stats.c(trial1:trial2)==1);    
        rewarded(:,idx)=rewarded(:,idx)+(stats.r(trial1:trial2)==1);
    end
end
        
probl=nan(1+trials_forw,numblockType); %prob of choosing left
probr=nan(1+trials_forw,numblockType); %prob of choosing right
probreward=nan(1+trials_forw,numblockType); %prob of rewarded
for j=1:numblockType
    if numtrans(j)>0
        probl(:,j)=chosel(:,j)/numtrans(j);
        probr(:,j)=choser(:,j)/numtrans(j);
        probreward(:,j)=rewarded(:,j)/numtrans(j);
    end
end

n=[0:trials_forw]';

output.n=n;
output.probl=probl;
output.probr=probr;
output.probreward=probreward;
output.numblockType=numblockType;
output.blockType=blockType;

end


    