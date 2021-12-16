function [qpar, negloglike, bic, nlike]=fit_fun(stats,fit_fun,initpar,lb,ub)    
% % fit_fun %
%PURPOSE:   Fit the choice behavior to a model using BADs optimizer
%AUTHORS:   H Atilgan & AC Kwan 03172020
%
%INPUT ARGUMENTS
%   stats:      stats of the game
%   fit_fun:    the model to fit, e.g., Q_RPEfun
%   initpar:    initial values for the parameters
%   lb:         lower bound values for the parameters (if used)
%   ub:         upper bound values for the parameters (if used)
%
%OUTPUT ARGUMENTS
%   qpar:       extracted parameters
%   negloglike: negative log likelihood
%   bic:        Bayesian information criterion
%   nlike:      normalized likelihood

%%
maxit  =1e6;
maxeval=1e6;
% op=optimset('fminsearch');
% op.MaxIter=maxit;
% op.MaxFunEvals=maxeval;

% % Bads optimizer params - do not make any changes
op = bads('defaults');
op.MaxIter=maxit;
op.MaxFunEvals=maxeval;
% op.uncertaintyHadling =1; % to tell BADS that optimization is noisy
% op.NoiseSize = 0.5;
 op.NoiseFinalSamples =50; % default is 10

func_handle = str2func(fit_fun);
if ~exist('lb','var')
    [qpar, negloglike, exitflag]= bads(func_handle,initpar,[],[],[],[],[],op,stats);
   % Default version: bads(fun,X0,LB,UB,PLB,PUB,NONBCON,OPTIONS,data)
    %[qpar, negloglike, exitflag]=fminsearch(func_handle, initpar, op, stats);
else
     [qpar, negloglike, exitflag]= bads(func_handle,initpar,lb,ub,lb,ub,[],op,stats);
    %[qpar, negloglike, exitflag]=fmincon(func_handle, initpar, [], [], [], [], lb, ub, [], op, stats);
end

if exitflag==0
    qpar=nan(size(qpar));   %did not converge to a solution, so return NaN
    negloglike=nan;
end

%% BIC, bayesian information criterion
%BIC = -logL + klogN
%L = negative log-likelihood, k = number of parameters, N = number of trials
%larger BIC value is worse because more parameters is worse, obviously
bic = negloglike + numel(initpar)*log(sum(~isnan(stats.c)));

%% Normalized likelihood 
%(Ito and Doya, PLoS Comp Biol, 2015)
nlike = exp(-1*negloglike)^(1/sum(~isnan(stats.c)));

end
