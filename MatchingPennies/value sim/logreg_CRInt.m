function [output, negloglike, bic, nlike]=logreg_CRInt(stats,x,step_back) 
% % logreg_CRInt %
%PURPOSE:   Logistic regression analysis of choice behavior
%           regressors are choice, and choice x reward
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   stats:  stats of the game
%   x:      which player? (1 or 2)
%   step_back:  how many trials back to analyze
%
%OUTPUT ARGUMENTS
%   output:     numbers used to plot figure
%   negloglike: negative log likelihood
%   bic:        Bayesian information criterion
%   nlike:      normalized likelihood
%
% To plot the output, use plot_logreg().

%%
c = stats.c(:,x);  %take player x's choice history, use all including NaN (do not concatenate or skip trials)
r = stats.r(:,x);  %take player x's reward history, use all including NaN (do not concatenate or skip trials)

% create regressor vectors
% double-check, these parameters should yield uncorrelated regressors
% corrcoef(Ch,RC)

Ch=zeros(length(c),1);   % right = 1; left = -1
Ch=1*(c==1) + (-1)*(c==-1);

Re=zeros(length(c),1);   % reward = 1; unrewarded = -1
Re=r.*(r>0) + (-1)*(r==0);  %allows for reward size r!=1
% Re=1*(r==1) + (-1)*(r==0); 

RC=zeros(length(c),1);   % choice x reward interaction (right/reward=1; right/unreward=-1; left/reward=-1; left/unreward=1)
RC=Ch.*Re;

clear Re;

% generate regressor matrix
rmat=zeros(length(Ch)-step_back,2*step_back);
for i=1+step_back:length(Ch)
    for j=1:step_back
        rmat(i-step_back,j)=Ch(i-j);
        rmat(i-step_back,j+step_back)=RC(i-j);
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
[b,~,regstats] =glmfit(rmat_fit, c_fit, 'binomial', 'link', 'logit');

output.n=-1:-1:-step_back;

output.b_bias=b(1);
output.pval_bias=regstats.p(1);

output.b_coeff(:,1)=b(2:step_back+1);
output.pval_coeff(:,1)=regstats.p(2:step_back+1);
output.b_label{1}='Choice';

output.b_coeff(:,2)=b(step_back+2:2*step_back+1);
output.pval_coeff(:,2)=regstats.p(step_back+2:2*step_back+1);
output.b_label{2}='Choice x Reward';

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
    b_temp = b_temp + nansum(output.b_coeff(1:backIdx,1).*Ch(k-1:-1:k-backIdx));  %choice terms
    b_temp = b_temp + nansum(output.b_coeff(1:backIdx,2).*RC(k-1:-1:k-backIdx));  %choice x reward terms
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

%% BIC, bayesian information criterion
%BIC = -logL + klogN
%L = negative log-likelihood, k = number of parameters, N = number of trials
%... larger value = worse because more parameters is worse, obviously
bic = negloglike + numel(b)*log(sum(~isnan(c)));

%% Normalized likelihood 
%(Ito and Doya, PLoS Comp Biol, 2015)
nlike = exp(-1*negloglike)^(1/sum(~isnan(c)));

