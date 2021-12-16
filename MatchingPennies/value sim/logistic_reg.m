function [output, negloglike]=logistic_reg(stats,step_back) 
% % logistic_reg %
%PURPOSE:   Logistic regression analysis of choice behavior
%AUTHORS:   AC Kwan 170518
%changed a little to fit my code / Hongli, 180507
%
%INPUT ARGUMENTS
%   stats:  stats of the game
%   x:      which player? (1 or 2)
%   step_back:  how many trials back to analyze
%
%OUTPUT ARGUMENTS
%   output:     numbers used to plot figure
%   negloglike: negative log likelihood
%
% To plot the output, use plot_logreg().

%%
c = stats.c;  %take player x's choice history, use all including NaN (do not concatenate or skip trials)
r = stats.r;  %take player x's reward history, use all including NaN (do not concatenate or skip trials)

% create regressor vectors - rewarded/right = 1; rewarded/left = -1
YR=zeros(length(c),1);   % rewarded
YR=1*((c==1) & (r==1)) + (-1)*((c==-1) & (r==1));

NR=zeros(length(c),1);   % unrewarded
NR=1*((c==1) & (r==0)) + (-1)*((c==-1) & (r==0));

% generate regressor matrix
rmat=zeros(length(NR)-step_back,2*step_back);
for i=1+step_back:length(NR)
    for j=1:step_back
        rmat(i-step_back,j)=YR(i-j);
        rmat(i-step_back,j+step_back)=NR(i-j);
    end
end

%which trial goes into the regression fit?
crit1 = [zeros(step_back,1); ones(numel(c)-step_back,1)]; %criterion 1: not the first several trials with no choice/outcome history
crit2 = ~isnan(c);  %criterion 2: a choice was made
goodTrial = crit1 & crit2;

c_fit=(c==1);    %convert to success = right for glmfit, binomial option
c_fit=c_fit(goodTrial);                       %choice, for trials that go into fit
rmat_fit=rmat(goodTrial(step_back+1:end),:);  %history, for trials that go into fit

%logistic regression
[b,~,stats] =glmfit(rmat_fit, c_fit, 'binomial', 'link', 'logit');
pval=stats.p;

output.n=-1:-1:-step_back;

output.b_bias=b(1);
output.b_reward=b(2:step_back+1);
output.b_unreward=b(step_back+2:2*step_back+1);

%% negative log likelihood

%probability of choosing left based on logistic regression
pr(1) = 0.5;
for k=2:numel(c)
    if k>step_back
        backIdx = step_back;  %if at least number of trials permit, use step_back
    else
        backIdx = k-1;        %else use what is available
    end
        
    b_temp = output.b_bias;    %bias term
    b_temp = b_temp + nansum(output.b_reward(1:backIdx).*YR(k-1:-1:k-backIdx));    %reward/right terms
    b_temp = b_temp + nansum(output.b_unreward(1:backIdx).*NR(k-1:-1:k-backIdx));  %unreward/right terms
    pr(k) = exp(b_temp)/(1+exp(b_temp));
end

%calculate negative log likelihood
negloglike = 0;
for k=1:numel(c)
    if c(k)==1  %choose right
        negloglike = negloglike - log(pr(k));
    elseif c(k)==-1  %choose left
        negloglike = negloglike - log(1-pr(k));
    end
end

