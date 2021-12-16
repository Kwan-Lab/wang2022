function negloglike=FQfun_withbeta(xpar,dat)
% % FQfun_withbeta %
%PURPOSE:   Function for maximum likelihood estimation, called by
%           fit_fun().
% reference: Katahira 2015. Set kapa1 = 1, kapa2 and beta as a free
% parameter

%INPUT ARGUMENTS
%   xpar:       alpha1, alpha2, kappa1, kappa2
%   dat:        data
%               dat(:,1) = choice vector
%               dat(:,2) = reward vector
%
%OUTPUT ARGUMENTS
%   negloglike:      the negative log-likelihood to be minimized

%%
a1=xpar(1);
a2=a1;
% k1=xpar(2);
k1 = 2;  % set kappa1 = 2
k2=xpar(2);
beta = xpar(3);
nt=size(dat,1);
negloglike=0;

v_right=0;
v_left=0;      

for k=1:nt
  
    pright=exp(beta*v_right)/(exp(beta*v_right)+exp(beta*v_left));
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
    
    if dat(k,1)==1      % chose right
        if dat(k,2)==1  %reward
            v_right = (1-a1)*v_right+a1*k1;
            v_left = (1-a2)*v_left;
        elseif dat(k,2)==0            %no reward
            v_right = (1-a1)*v_right-a1*k2;
            v_left = (1-a2)*v_left;
        end
    elseif dat(k,1)==-1    % choose left
        if dat(k,2)==1     %reward
            v_left = (1-a1)*v_left+a1*k1;
            v_right = (1-a2)*v_right;
        elseif dat(k,2)==0            %no reward
            v_left = (1-a1)*v_left-a1*k2;
            v_right = (1-a2)*v_right;
        end
    end   
end
