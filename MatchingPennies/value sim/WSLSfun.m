function negloglike=WSLSfun(xpar,dat)
% % WSLSfun %
%PURPOSE:   Function for maximum likelihood estimation, called by
%           fit_fun().
%
%INPUT ARGUMENTS
%   xpar:       p (probability of following win-stay, lose-switch
%   dat:        data
%               dat(:,1) = choice vector
%               dat(:,2) = reward vector
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized

%%
p=xpar(1);
nt=size(dat,1);
negloglike=0;

for k=1:nt
  
    %probability of choosing right
    if k==1  %first trial
        pright=0.5;
    elseif dat(k-1,2)==1 && dat(k-1,1)==1  %last trial = win + right
        pright=p;
    elseif dat(k-1,2)==1 && dat(k-1,1)==-1 %last trial = win + left
        pright=1-p;
    elseif dat(k-1,2)==0 && dat(k-1,1)==1  %last trial = lose + right
        pright=1-p;
    elseif dat(k-1,2)==0 && dat(k-1,1)==-1 %last trial = lose + left
        pright=p;
    else                                   %last trial is a miss
        pright=0.5;
    end
    pleft=1-pright;
        
    if pright==0, 
        pright=realmin;   % Smallest positive normalized floating point number, because otherwise log(zero) is -Inf
    end        
    if pleft==0, 
        pleft=realmin;
    end;            

    %compare with actual choice to calculate log-likelihood
    if dat(k,1)==1
        logp=log(pright);
    elseif dat(k,1)==-1
        logp=log(pleft);
    else
        logp=0;
    end;   
    negloglike=negloglike-logp;  % calculate log likelihood
    
end
