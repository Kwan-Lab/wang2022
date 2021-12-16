function negloglike=funbayesian_twoparams(xpar,stats)
% % funbayesian_twoparams %
%PURPOSE:   Function for maximum likelihood estimation, called by fit_fun().
%
%INPUT ARGUMENTS
%   xpar:       2 bayesian params
%   stats:      info about animal/agent's performance
%
% this model could be used around switch, say +- 10 trials
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized

%% parameters to estimate

% bayesian parameters
p_trans=xpar(1);
p_crit=xpar(2);

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

% how many correct choices have occurred, after the last rule switch
correctChoice = 0;

%% identify when the switches occur, for Bayesian part

% rule from ground truth (not observable)
rule = stats.rule;

% find the switch trials
x = rule';
l = diff([0 find(x(1:end-1) ~= x(2:end)) length(x)]); %run-length encoding
swTrial = cumsum(l)+1;                  %trial where there were switches

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
    
    % calculate the criterion
    if ismember(k, swTrial)  %reset if a rule switch occurred
        correctChoice = 0;
    else
        if stats.hr_side(k) == stats.c(k)
            correctChoice = correctChoice + 1;
        end
    end
    
    % based on whether criterion is met, what is the transition probability?
    if correctChoice < criterion
        p = p_trans; % criterion not met, transition expected at baseline rate
    else
        p = p_crit;  % criterion met, transition expected according to param
    end
 
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

