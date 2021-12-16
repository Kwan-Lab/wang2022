function output=choice_switch(stats,trials_back)
% % choice_switch %
%PURPOSE:   Analyze choice behavior around switches of all possible types
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   stats:  stats of the task
%   trials_back: how many trials around switch to plot
%
%OUTPUT ARGUMENTS
%   output:     numbers used to plot figure
%
% To plot the output, use plot_switch().

%%
%run-length encoding to find when switches occur
x0=stats.rule';
l = diff([ 0 find(x0(1:end-1) ~= x0(2:end)) length(x0) ]);
v = x0([ find(x0(1:end-1) ~= x0(2:end)) length(x0) ]);

%number of reward probability pairs
num_probpairs=size(unique(stats.rewardprob,'rows'),1);

%find the unique reward probabilities transitions
%reward prob left/pre-switch, right/pre-switch -> left/post-switch, right/post-switch
transType=[];
for i=1:num_probpairs
    for j=1:num_probpairs
        if i~=j
            transType=[transType; i j]; %all the possible rule transitions
        end
    end    
end
numtransType=size(transType,1);

chosel=zeros(1+2*trials_back,numtransType); %prob of choosing left
choser=zeros(1+2*trials_back,numtransType); %prob of choosing right
numtrans=zeros(1,numtransType);    %number of that type of transition

for i=1:numel(l)-1    
    %which kind of transition is it?
    idx=find(transType(:,1)==v(i) & transType(:,2)==v(i+1));
    
    %what were the choices around that transition
    trial1=(sum(l(1:i))+1)-trials_back;
    trial2=(sum(l(1:i))+1)+trials_back;
    if trial2<=numel(stats.c)
        numtrans(idx)=numtrans(idx)+1;
        chosel(:,idx)=chosel(:,idx)+(stats.c(trial1:trial2)==-1);
        choser(:,idx)=choser(:,idx)+(stats.c(trial1:trial2)==1);    
    end
end

probl=nan(1+2*trials_back,numtransType); %prob of choosing left
probr=nan(1+2*trials_back,numtransType); %prob of choosing right
for j=1:numtransType
    if numtrans(j)>0
        probl(:,j)=chosel(:,j)/numtrans(j);
        probr(:,j)=choser(:,j)/numtrans(j);
    end
end

n=[-trials_back:1:trials_back]';

output.n=n;
output.probl=probl;
output.probr=probr;
output.numtransType=numtransType;
output.transType=transType;

end


    