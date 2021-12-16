function [output, negloglike, bic, nlike]=bandit_logreg_RCUC_LR(stats,step_back) 
% % logreg_RCUC_LR %
%PURPOSE:   Logistic regression analysis of choice behavior
%           regressors are past left rewarded choices, right rewarded choices
%           left unrewarded choices, and right unrewarded choices,
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   stats:  stats of the game
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

% create regressor vectors: 
% YR == rewarded choices, NR == unrewarded choices
% reason for using these two vectors is that they are uncorrelated: corrcoef(YR,NR)
nTrial = numel(stats.c);

YR_right=zeros(length(stats.c),1);   % rewarded/right = 1
YR_right=stats.r.*((stats.c==1) & (stats.r>0));       %allows for reward size r!=1

YR_left=zeros(length(stats.c),1);    % rewarded/left = -1
YR_left=(-1)*stats.r.*((stats.c==-1) & (stats.r>0));  %allows for reward size r!=1

NR_right=zeros(length(stats.c),1);   % unrewarded/right = 1
NR_right=1*((stats.c==1) & (stats.r==0));

NR_left=zeros(length(stats.c),1);    % unrewarded/left = -1
NR_left=(-1)*((stats.c==-1) & (stats.r==0));

% generate regressor matrix
rmat=zeros(nTrial-step_back,4*step_back);
for i=1+step_back:nTrial
    for j=1:step_back
        rmat(i-step_back,j)=YR_right(i-j);
        rmat(i-step_back,j+step_back)=YR_left(i-j);
        rmat(i-step_back,j+2*step_back)=NR_right(i-j);
        rmat(i-step_back,j+3*step_back)=NR_left(i-j);
    end
end

%which trial goes into the regression fit?
crit1 = [zeros(step_back,1); ones(nTrial-step_back,1)]; %criterion 1: not the first several trials with no choice/outcome history
crit2 = ~isnan(stats.c);  %criterion 2: a choice was made
goodTrial = crit1 & crit2;

c_fit=(stats.c==1);    %convert to success = right for glmfit, binomial option
c_fit=c_fit(goodTrial);                       %choice, for trials that go into fit
rmat_fit=rmat(goodTrial(step_back+1:end),:);  %history, for trials that go into fit

%logistic regression
opts = statset('glmfit');
opts.MaxIter = 1000; % default value for glmfit is 100.
[b,~,regstats] =glmfit(rmat_fit, c_fit, 'binomial', 'link', 'logit', 'options', opts);

output.n=-1:-1:-step_back;

output.b_bias=b(1);
output.pval_bias=regstats.p(1);

output.b_coeff(:,1)=b(2:step_back+1);
output.pval_coeff(:,1)=regstats.p(2:step_back+1);
output.b_label{1}='Rewarded right choice';

output.b_coeff(:,2)=b(step_back+2:2*step_back+1);
output.pval_coeff(:,2)=regstats.p(step_back+2:2*step_back+1);
output.b_label{2}='Rewarded left choice';

output.b_coeff(:,3)=b(2*step_back+2:3*step_back+1);
output.pval_coeff(:,3)=regstats.p(2*step_back+2:3*step_back+1);
output.b_label{3}='Unrewarded right choice';

output.b_coeff(:,4)=b(3*step_back+2:4*step_back+1);
output.pval_coeff(:,4)=regstats.p(3*step_back+2:4*step_back+1);
output.b_label{4}='Unrewarded left choice';

%% negative log likelihood

%probability of choosing left based on logistic regression
pr(1) = 0.5;
for k=2:nTrial
    if k>step_back
        backIdx = step_back;  %if at least number of trials permit, use step_back
    else
        backIdx = k-1;        %else use what is available
    end
        
    b_temp = output.b_bias;    %bias term
    b_temp = b_temp + nansum(output.b_coeff(1:backIdx,1).*YR_right(k-1:-1:k-backIdx)); %reward/right terms
    b_temp = b_temp + nansum(output.b_coeff(1:backIdx,2).*YR_left(k-1:-1:k-backIdx));  %reward/left terms
    b_temp = b_temp + nansum(output.b_coeff(1:backIdx,3).*NR_right(k-1:-1:k-backIdx)); %unreward/right terms
    b_temp = b_temp + nansum(output.b_coeff(1:backIdx,4).*NR_left(k-1:-1:k-backIdx));  %unreward/left terms
    pr(k) = exp(b_temp)/(1+exp(b_temp));
end

%calculate negative log likelihood
negloglike = 0;
for k=1:nTrial
    if stats.c(k)==1  %choose right
        negloglike = negloglike - log(pr(k));
    elseif stats.c(k)==-1  %choose left
        negloglike = negloglike - log(1-pr(k));
    end
end

%% BIC, bayesian information criterion
%BIC = -logL + klogN
%L = negative log-likelihood, k = number of parameters, N = number of trials
%... larger value = worse because more parameters is worse, obviously
bic = negloglike + numel(b)*log(sum(~isnan(stats.c)));

%% Normalized likelihood 
%(Ito and Doya, PLoS Comp Biol, 2015)
nlike = exp(-1*negloglike)^(1/sum(~isnan(stats.c)));

