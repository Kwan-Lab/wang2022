function negloglike=FQfun_withbetaCA(xpar,dat)
% % FQfun_withbeta %
%PURPOSE:   Function for maximum likelihood estimation, called by
%           fit_fun().
% reference: Katahira 2015. Set kapa1 = 1, kapa2 and beta as a free
% parameter: CA: choice autocorrelation

%INPUT ARGUMENTS
%   xpar:       alpha1, kappa2, beta, tau, phi
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
%deltaQ = xpar(4);
tau = xpar(4);
phi = xpar(5);
nt=size(dat,1);
negloglike=0;

v_right=0;
v_left=0;      

c_left = 0;
c_right = 0;

for k=1:nt
  
    pright=exp(beta*(v_right + phi*c_right))/(exp(beta*(v_right + phi*c_right))+exp(beta*(v_left+ phi*c_left)));
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
        c_right = (1-tau)*c_right + tau;
        c_left = (1-tau)*c_left;
        if dat(k,2)==1  %reward
            v_right = (1-a1)*v_right+a1*k1;
            v_left = (1-a2)*v_left;
            
        elseif dat(k,2)==0            %no reward
            v_right = (1-a1)*v_right-a1*k2;
            v_left = (1-a2)*v_left;
        end
    elseif dat(k,1)==-1    % choose left
        c_right = (1-tau)*c_right;
        c_left = (1-tau)*c_left + tau;
        if dat(k,2)==1     %reward
            v_left = (1-a1)*v_left+a1*k1;
            v_right = (1-a2)*v_right;
        elseif dat(k,2)==0            %no reward
            v_left = (1-a1)*v_left-a1*k2;
            v_right = (1-a2)*v_right;
        end
    end   
end
