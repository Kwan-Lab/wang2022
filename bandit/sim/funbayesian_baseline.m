function negloglike=funbayesian_baseline(xpar,stats)
% % funbayesian_baseline %
%PURPOSE:   Function for maximum likelihood estimation, called by fit_fun().
%
%INPUT ARGUMENTS
%   xpar:       1 bayesian param
%   stats:      info about animal/agent's performance
%
% this model could be used around switch, say +- 10 trials
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized

%% parameters to estimate

% bayesian parameters
p=xpar(1);

% task parameters
probreward.high = 0.7;
probreward.low = 0.1;
criterion = 10;

%% initialize

nt=size(stats.c,1);
negloglike=0;

% initialize bayesian posteriors
pleft = 0.5;   %probability of being in the left sub-state
pright = 0.5;

%% start the loop
for k=1:nt
    
    %smallest positive normalized floating point number, because otherwise log(zero) is -Inf
    if pleft==0
        pleft=realmin;
    end
    if pright==0
        pright=realmin;
    end
    
    % calculate log likelihood
    if stats.c(k)==1
        logp=log(pright);
    elseif stats.c(k)==-1
        logp=log(pleft);
    else
        logp=0;
    end
    negloglike=negloglike-logp;
        
    %% updating the bayesian posterior belief
    
    % calculate posterior probabilities
    if stats.c(k) == -1 % choose left
        if stats.r(k) == 1 % get reward
            pleft_post = probreward.high * ((1-p) * pleft + p * pright);
            pright_post = probreward.low  * (p * pleft + (1-p) * pright);
        else % no reward
            pleft_post = (1 - probreward.high) * ((1-p) * pleft + p * pright);
            pright_post = (1 - probreward.low)  * (p * pleft + (1-p) * pright);
        end
        
    elseif stats.c(k) == 1  % choose right
        if stats.r(k) == 1  % get reward
            pleft_post = probreward.low  * ((1-p) * pleft + p * pright);
            pright_post = probreward.high * (p * pleft + (1-p) * pright);
        else % no reward
            pleft_post = (1 - probreward.low)  * ((1-p) * pleft + p * pright);
            pright_post = (1 - probreward.high) * (p * pleft + (1-p) * pright);
        end
        
    else  % miss
        pleft_post = ((1-p) * pleft + p * pright);
        pright_post = (p * pleft + (1-p) * pright);
        
    end
    
    pleft = pleft_post / (pleft_post + pright_post);
    pright = 1 - pleft;

end

end

