function reg_results = getBootstrp(reg_all,ifpval, alpha)

% if ifpval = 0.05:  bootstrap the p value to get baseline significant
% sessions, significant level
%    ifpval = 0:  bootstrap the coefficient to get average and errorbar of the coefficients
%    alpha:       if ifpval != 0, alpha determines the significance, ifpval
%    determines the significance of the coefficient

numRep = 1000;
bootSig = [];
CI = 1-alpha;

if ifpval == 0
% bootstrap the coefficient
    
    all_coeff = reg_all.coeff;
% using bootstrap to get p value

    for kk = 1:numRep
        if length(size(all_coeff)) == 3  % coefficient for regression of pupil average
            drawIndex=randsample(size(all_coeff,3),size(all_coeff,3),'true'); %each time draw a set of events with replacement
        elseif length(size(all_coeff)) == 2 % coefficient for regression of pupil response
            drawIndex=randsample(size(all_coeff,1),size(all_coeff,1),'true');
        end
        drawSig=[]; drawTime=[];
        for k=1:numel(drawIndex)
            if length(size(all_coeff)) == 3  % coefficient for regression of pupil average
                drawSig=cat(3,drawSig, all_coeff(:,:,drawIndex(k)));
            elseif length(size(all_coeff)) == 2 % coefficient for regression of pupil response
                drawSig = cat(3, drawSig, all_coeff(drawIndex(k),:));
            end
        end
    
    %go through each time bin, and average
        bootSig = cat(3, bootSig,nanmean(drawSig,3));
    end

%bootstrap mean and confidence interval
    reg_all.bootSig = bootSig;
    reg_all.coeff_bootave=nanmean(bootSig,3);
    reg_all.bootlow=prctile(bootSig,0.5*(1-CI)*100,3);
    reg_all.boothigh=prctile(bootSig,(1-0.5*(1-CI))*100,3);
    reg_results = reg_all;

else
    % bootstrap the p value, for pupil average regression only
    % loop through field
    
    pval = reg_all;
    field = fieldnames(reg_all);
    for uu = 1:length(field)
        Sig = [];
        for kk = 1:numRep
            drawIndex=randsample(size(pval.(field{uu}),3),size(pval.(field{uu}),3),'true'); %each time draw a set of events with replacement
            drawSig=[]; drawTime=[];
            for k=1:numel(drawIndex)
                drawSig=cat(3,drawSig, pval.(field{uu})(:,:,drawIndex(k)));
            end
            
           % get the significant sessions
           sigSessions = sum(drawSig(:,uu+1,:) < ifpval,3) / size(drawSig,3) * 100;
           Sig = [Sig, sigSessions];
           
        end
        
        % get median of Sig, and 95% CI
        reg_results.(field{uu}).bootmedian = nanmedian(Sig, 2);
        reg_results.(field{uu}).bootlow = prctile(Sig,0.5*(1-CI)*100,2);
        reg_results.(field{uu}).boothigh = prctile(Sig,(1-0.5*(1-CI))*100,2);
        
        %get the corresponding p value
        pvalPred = pval.(field{uu});
        percentage = 100*sum(pvalPred(:,uu+1,:)<ifpval,3)/size(pvalPred,3);
        reg_results.(field{uu}).pval = zeros(1,length(percentage));
        for vv = 1:length(percentage)
            reg_results.(field{uu}).pval(vv) = sum(Sig(vv,:)> percentage(vv))/size(Sig,2);
        end
    end
end


end

