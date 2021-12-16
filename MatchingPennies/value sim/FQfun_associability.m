function negloglike=FQfun_associability(xpar,dat)
% % FQfun_associability %
%PURPOSE:   Function for maximum likelihood estimation, called by
%           fit_fun().
% reference: Li 2011, Nature Neuroscience.
% parameter

%INPUT ARGUMENTS
%   xpar:       kappa, eta, alpha0
%   dat:        data
%               dat(:,1) = choice vector
%               dat(:,2) = reward vector
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized

%%
kappa=xpar(1);
eta = xpar(2);
alpha0 = xpar(3);

nt=size(dat,1);
negloglike=0;

v_right=0;
v_left=0;      
alpha = 0;
for k=1:nt
  
    pright=exp(v_right)/(exp(v_right)+exp(v_left));
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
    
    if dat(k,1)==1      % chose right, update the v_right and alpha
        v_right = v_right+alpha*kappa*(dat(k,2)-v_right);
        alpha = (1-eta)*alpha + eta * abs(dat(k,2) - v_right);
    elseif dat(k,1)==-1    % choose left
        v_left = v_left+alpha*kappa*(dat(k,2)-v_left);
        alpha = (1-eta)*alpha + eta * abs(dat(k,2) - v_left);
    end   
end
