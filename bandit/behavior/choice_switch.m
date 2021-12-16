function output=choice_switch(stats,trials_back)
% % choice_switch %
%PURPOSE:   Analyze choice behavior around switches
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
%how many different types of rule transitions?
numruletrans=size(stats.ruletransList,1);

if numruletrans==0 % No switch occured
    disp('No switch of reward probabilities in this session');
    output = 0;
else
    chosel=zeros(1+2*trials_back,numruletrans); %prob of choosing left
    choser=zeros(1+2*trials_back,numruletrans); %prob of choosing right
    choseneither=zeros(1+2*trials_back,numruletrans); %prob of not choosing (miss)
    
    numtrans=zeros(1,numruletrans);    %number of that type of transition
    
    for i=1:numel(stats.blockLength)-1
        %which kind of transition is it?
        idx=stats.blockTrans(i);
        
        if ~isnan(idx)   %NaN could arise if we are using merge_session
            
            %what were the choices around that transition
            %note: 0 is the first trial with the switched probabilities
            trial1=(sum(stats.blockLength(1:i))+1)-trials_back;  %trials before switch
            trial2=(sum(stats.blockLength(1:i))+1)+trials_back;  %trials after switch
            
            if trial2<=numel(stats.c)   %if trial after switch does not exceed end of session
                numtrans(idx)=numtrans(idx)+1;
                chosel(:,idx)=chosel(:,idx)+(stats.c(trial1:trial2)==-1);
                choser(:,idx)=choser(:,idx)+(stats.c(trial1:trial2)==1);
                choseneither(:,idx)=choseneither(:,idx)+isnan(stats.c(trial1:trial2));
            end
        end
    end
    
    probl=nan(1+2*trials_back,numruletrans); %prob of choosing left
    probr=nan(1+2*trials_back,numruletrans); %prob of choosing right
    probneither=nan(1+2*trials_back,numruletrans); %prob of choosing neither
    for j=1:numruletrans
        if numtrans(j)>0
            probl(:,j)=chosel(:,j)/numtrans(j);
            probr(:,j)=choser(:,j)/numtrans(j);
            probneither(:,j)=choseneither(:,j)/numtrans(j);
        end
    end
    
    n=[-trials_back:1:trials_back]';
    
    output.n=n;
    output.probl=probl;
    output.probr=probr;
    output.probneither=probneither;
    output.numtransType=numruletrans;
    output.transType=stats.ruletransList;
end

end


