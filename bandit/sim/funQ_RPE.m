function negloglike=funQ_RPE(xpar,stats)
% % funQ_RPE %
%PURPOSE:   Function for maximum likelihood estimation, called by
%           fit_fun().
%
%INPUT ARGUMENTS
%   xpar:       alpha, beta
%   stats:      info about animal/agent's performance
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized

%%
alpha=xpar(1);
beta=xpar(2);

nt=size(stats.c,1);
negloglike=0;

v_right=0.5;
v_left=0.5;      

for k=1:nt
  
    pright=exp(beta*v_right)/(exp(beta*v_right)+exp(beta*v_left));
    pleft=1-pright;
        
    if pright==0 
        pright=realmin;   % Smallest positive normalized floating point number, because otherwise log(zero) is -Inf
    end        
    if pleft==0 
        pleft=realmin;
    end            
  
    %compare with actual choice to calculate log-likelihood
    if stats.c(k)==1
        logp=log(pright);
    elseif stats.c(k)==-1
        logp=log(pleft);
    else
        logp=0;
    end  
    negloglike=negloglike-logp;  % calculate log likelihood
    
    % update value for the performed action
    if stats.c(k)==1      %chose right
        v_right=v_right+alpha*(stats.r(k)-v_right);   
    elseif stats.c(k)==-1 %chose left
        v_left=v_left+alpha*(stats.r(k)-v_left);
    end   
end

end
