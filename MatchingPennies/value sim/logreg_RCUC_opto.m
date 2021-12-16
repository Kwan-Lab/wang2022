function [output, negloglike]=logreg_RCUC_opto(stats,x,step_back)
% % logreg_RCUC %
%PURPOSE:   Logistic regression analysis of choice behavior
%           regressors are past rewarded choices and unrewarded choices
%           (fraction of trials have perturbed learning parameters)
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
%
% To plot the output, use plot_logreg().

%%
c = stats.c(:,x);  %take player x's choice history, use all including NaN (do not concatenate or skip trials)
r = stats.r(:,x);  %take player x's reward history, use all including NaN (do not concatenate or skip trials)

% create regressor vectors: rewarded/right = 1; rewarded/left = -1
YR=zeros(length(c),1);
YR=r.*((c==1) & (r>0)) + (-1*r).*((c==-1) & (r>0));  %allows for reward size r!=1
% YR=1*((c==1) & (r==1)) + (-1).*((c==-1) & (r==1));

NR=zeros(length(c),1);   % unrewarded/right = 1; unrewarded/left = -1
NR=1*((c==1) & (r==0)) + (-1)*((c==-1) & (r==0));

% generate regressor matrix
rmat=zeros(length(NR)-step_back,2*step_back);
for i=1+step_back:length(NR)
    for j=1:step_back
        rmat(i-step_back,j)=YR(i-j);
        rmat(i-step_back,j+step_back)=NR(i-j);
    end
end

%% logistic regression for trials without perturbation (control)
opto = stats.opto(:,x);  %take player x's perturbation history, use all including NaN (do not concatenate or skip trials)

%which trial goes into the regression fit?
crit1 = [zeros(step_back,1); ones(numel(c)-step_back,1)]; %criterion 1: not the first several trials with no choice/outcome history
crit2 = ~isnan(c);  %criterion 2: a choice was made
crit3 = (opto==0);  %criterion 3: no perturbation
goodTrial = crit1 & crit2 & crit3;

c_ctrl=(c==1);    %convert to success = right for glmfit, binomial option
c_ctrl=c_ctrl(goodTrial);                      %choice, for trials that go into fit
rmat_ctrl=rmat(goodTrial(step_back+1:end),:);  %history, for trials that go into fit

[b_ctrl,~,stats] =glmfit(rmat_ctrl, c_ctrl, 'binomial', 'link', 'logit');
output.pval_ctrl=stats.p;

%% logistic regression for trials with perturbation

%which trial goes into the regression fit?
crit1 = [zeros(step_back,1); ones(numel(c)-step_back,1)]; %criterion 1: not the first several trials with no choice/outcome history
crit2 = ~isnan(c);  %criterion 2: achoice was made
crit3 = (opto==1);  %criterion 3: with perturbation!
goodTrial = crit1 & crit2 & crit3;

c_opto=(c==1);    %convert to success = right for glmfit, binomial option
c_opto=c_opto(goodTrial);                      %choice, for trials that go into fit
rmat_opto=rmat(goodTrial(step_back+1:end),:);  %history, for trials that go into fit

[b_opto,~,stats] =glmfit(rmat_opto, c_opto, 'binomial', 'link', 'logit');
output.pval_opto=stats.p;

%%
output.n=-1:-1:-step_back;

output.b_bias_ctrl=b_ctrl(1);
output.b_coeff_ctrl(:,1)=b_ctrl(2:step_back+1);
output.b_label{1}='Rewarded choice (ctrl)';
output.b_coeff_ctrl(:,2)=b_ctrl(step_back+2:2*step_back+1);
output.b_label{2}='Unrewarded choice (ctrl)';

output.b_bias_opto=b_opto(1);
output.b_coeff_opto(:,1)=b_opto(2:step_back+1);
output.b_label{3}='Rewarded choice (opto)';
output.b_coeff_opto(:,2)=b_opto(step_back+2:2*step_back+1);
output.b_label{4}='Unrewarded choice (opto)';

%% negative log likelihood

%probability of choosing left based on logistic regression
pr(1) = 0.5;
for k=2:numel(c)
    if k>step_back
        backIdx = step_back;  %if at least number of trials permit, use step_back
    else
        backIdx = k-1;        %else use what is available
    end
    
    if opto(k)==0
        b_temp = output.b_bias_ctrl;    %bias term
        b_temp = b_temp + nansum(output.b_coeff_ctrl(1:backIdx,1).*YR(k-1:-1:k-backIdx));  %reward/right terms
        b_temp = b_temp + nansum(output.b_coeff_ctrl(1:backIdx,2).*NR(k-1:-1:k-backIdx));  %unreward/right terms
    elseif opto(k)==1
        b_temp = output.b_bias_opto;    %bias term
        b_temp = b_temp + nansum(output.b_coeff_opto(1:backIdx,1).*YR(k-1:-1:k-backIdx));  %reward/right terms
        b_temp = b_temp + nansum(output.b_coeff_opto(1:backIdx,2).*NR(k-1:-1:k-backIdx));  %unreward/right terms
    end
    
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

