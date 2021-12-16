function [qpar, negloglike, bic, nlike]=fit_fun(stats,fit_fun,initpar,lb,ub)    
% % fit_fun %
%PURPOSE:   Fit the choice behavior to a model using maximum likelihood estimate
%AUTHORS:   AC Kwan 170518
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
maxit=1e6;
maxeval=1e6;
op=optimset('fminsearch');
op.MaxIter=maxit;
op.MaxFunEvals=maxeval;

func_handle = str2func(fit_fun);
if strcmp(fit_fun, 'funHybrid_FQBayesian')
    if ~exist('lb','var')
        [qpar, negloglike, exitflag]=fminsearch(func_handle, initpar, op, stats);
    else
        [qpar, negloglike, exitflag]=fmincon(func_handle, initpar, [], [], [], [], lb, ub, [], op, stats);
    end
elseif strcmp(fit_fun, 'funHybrid_FQBayesian_all')
    if ~exist('lb','var')
        [qpar, negloglike, exitflag]=fminsearch(func_handle, initpar, op, stats);
    else
        [qpar, negloglike, exitflag]=fmincon(func_handle, initpar, [], [], [], [], lb, ub, [], op, stats);
    end
elseif strcmp(fit_fun, 'funHybrid_FQBayesian_FQ')
    if ~exist('lb','var')
        [qpar, negloglike, exitflag]=fminsearch(func_handle, initpar, op, stats);
    else
        [qpar, negloglike, exitflag]=fmincon(func_handle, initpar, [], [], [], [], lb, ub, [], op, stats);
    end
elseif strcmp(fit_fun, 'funFQ_RPE_CK_drift') || strcmp(fit_fun, 'funFQ_RPE_drift')
    if ~exist('lb','var')
        [qpar, negloglike, exitflag]=fminsearch(func_handle, initpar, op, stats);
    else
        [qpar, negloglike, exitflag]=fmincon(func_handle, initpar, [], [], [], [], lb, ub, [], op, stats);
    end
else
    if ~exist('lb','var')
        [qpar, negloglike, exitflag]=fminsearch(func_handle, initpar, op, [stats.c stats.r]);
    else
        [qpar, negloglike, exitflag]=fmincon(func_handle, initpar, [], [], [], [], lb, ub, [], op, [stats.c stats.r]);
    end
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
