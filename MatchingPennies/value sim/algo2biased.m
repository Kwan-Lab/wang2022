function stats=algo2biased(stats,params,x)
% % algo2biased %
%PURPOSE:   Simulate computer opponent, algorithm 2
%AUTHORS:   AC Kwan 180423
%
%INPUT ARGUMENTS
%   stats:  stats of the game thus far
%   params: parameters that define the player's strategy
%       trial_back      - number of trial back for calculating conditional probabilities
%       trial_history   - trials older than this number are not considered
%       pLeft           - test deviation from this probability to choose left
%   x:      which player? (1 or 2)
%
%OUTPUT ARGUMENTS
%   stats:  updated with player's probability to choose left for next step

% reference: Daeyeol Science 2014
%%

%consider the opponent's choice history
if x==1 %if player 1
    c = stats.c(~isnan(stats.c(:,2)),2);  %take player 2's choice history
    r = stats.r(~isnan(stats.r(:,2)),2);  %take player 2's reward history
elseif x==2 %if player 2
    c = stats.c(~isnan(stats.c(:,1)),1);  %take player 1's choice history
    r = stats.r(~isnan(stats.r(:,1)),1);  %take player 1's reward history
end

%consider only up to a certain number of trials for history
if numel(c) > params.trial_history
    c = c(end-params.trial_history+1:end);
    r = r(end-params.trial_history+1:end);
end

if numel(c) > params.trial_back  %if there is enough choice history
    %% trial back n = 0
    nom = sum(c==-1);
    denom = sum(c==-1)+sum(c==1);
    pval = myBinomTest(nom,denom,params.pLeft);
    pl = nom/denom;
    
    %% trial back n > 0
    for n = 1:params.trial_back
        % considering opponent's past choice(s)
        template = c(end-n+1:end);            %subject's past choice(s)
        idx = strfind(c(1:end-1)',template'); %search choice history for their occurrences
        
        nom = sum(c(idx+n)==-1);
        denom = sum(c(idx+n)==1)+sum(c(idx+n)==-1);
        pval = [pval myBinomTest(nom,denom,params.pLeft)];
        pl = [pl nom/denom];

        % considering both choice and reward
        template_r = r(end-n+1:end);            %subject's past reward(s)
        idx_r = strfind(r(1:end-1)',template_r'); %search reward history for their occurrences
        
        idx_cr = intersect(idx,idx_r);
        
        nom = sum(c(idx_cr+n)==-1);
        denom = sum(c(idx_cr+n)==1)+sum(c(idx_cr+n)==-1);
        pval = [pval myBinomTest(nom,denom,params.pLeft)];
        pl = [pl nom/denom];
    end
    
    %% which conditional probability is statistically significant?
    pl = pl(pval<0.05);
    
    %% player will choose based on conditional probability is furthest away from the specified probability
    if ~isempty(pl)
        [~,max_idx]=max(abs(pl - params.pLeft));
        dev = pl(max_idx)-params.pLeft;   %deviation from specified probability
        %stats.pl(stats.currTrial,x) = params.pLeft - dev;  %algorithm will strive for a non-match
        stats.pl(stats.currTrial,x) = 1-pl(max_idx);
        if stats.pl(stats.currTrial,x) > 1  %check to make sure we are within bounds
            stats.pl(stats.currTrial,x) = 1;
        elseif stats.pl(stats.currTrial,x) < 0
            stats.pl(stats.currTrial,x) = 0;
        end

%-------daeyeol
%         if pl(max_idx) > params.pLeft
%             stats.pl(stats.currTrial,x) = 0;
%         else
%             stats.pl(stats.currTrial,x) = 1;
%         end
    else
        stats.pl(stats.currTrial,x) = 1 - params.pLeft; %note! this only make sense for a diagonal-element only, zero-sum matrix
    end
    
else  %if there is not enough choice history
    stats.pl(stats.currTrial,x) = 1 - params.pLeft; %note! this only make sense for a diagonal-element only, zero-sum matrix
end

end
