function MP_plot_regrcoef_pupil_both(input,pvalThresh,tlabel,xtitle)
% % plot_regr %
%PURPOSE:   Plot the average coefficients from multiple linear regression
%AUTHORS:   AC Kwan 170515
%
%INPUT ARGUMENTS
%   input:        Structure generated by linear_regr().
%   pvalThresh:   Threshold value to deem whether p-values are significant
%   tlabel:       Text to put as title of the plot.
%   xtitle:       Text to put as the label for x-axis.

%% setup
t=input.left.regr_time;
dt=nanmean(diff(t));

nCells=numel(input.left);
for j=1:nCells
    if isfield(input.left, 'pval')
        pval.left(:,:,j)=input.left.pval;
        pval.right(:,:,j)=input.right.pval;
    end
    coeff.left(:,:,j) = nanmean(input.left.coeff,3);
    coeff.right(:,:,j) = nanmean(input.right.coeff,3);
end

nPredictor=input.left.numPredictor;

nback=input.left.nback;

if (input.left.interaction == true)
    if nPredictor == 2
        panelv = nPredictor + 1;    %plot extra row for the interaction terms
        nInteraction = 1;
    else
        panelv = nPredictor + 2;    %plot extra row for the interaction terms
        nInteraction = 2;
    end
else
    panelv = nPredictor;
    nInteraction = 0;
end

%% plot results

h = figure;

% determine the y limits, ignore the bias term
ymax = max(max(max(coeff.left(:,2:end))),max(max(coeff.right(:,2:end))));
ymin = min(min(min(coeff.left(:,2:end))),min(min(coeff.right(:,2:end))));

% ceil to 0.5
upperLimit = ceil(ymax*2)/2;
lowerLimit = floor(ymin*2)/2;
for l=1:nPredictor
    for k=1:1+nback
        currPredictor=1+k+(l-1)*(1+nback); %first +1 because first term is bias
        
        % determine whether to put x/y labels
        if nPredictor == 8  % for dQ and chosen Q
            subplot(2,4,currPredictor-1);hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 9 % for control plot
            subplot(3,3, currPredictor-1);hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 7  % for RPE/chosenQ regression
            subplot(2,4,currPredictor-1); hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 2 & nback == 0% for pos/neg RPE
            subplot(1,2, currPredictor-1); hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 14 % for future choice and reward
            subplot(4,4, currPredictor-1); hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 5  % for reward 
            subplot(2,3, currPredictor-1); hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        elseif nPredictor == 12 % for dQ and dC
            subplot(4,3, currPredictor-1); hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        else
            subplot(panelv,1+nback,currPredictor-1); hold on;
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback);
        end
        % patch([t(1) t(end) t(end) t(1)],[0 0 100*pvalThresh 100*pvalThresh],[0.7 0.7 0.7],'EdgeColor','none');
        plot(t,coeff.left(:,currPredictor,:),'k.-','MarkerSize',15);
        hold on;
        plot(t,coeff.right(:,currPredictor,:),'k--','MarkerSize',15);
        
        xlim([t(1) t(end)]);
        title(tlabel{currPredictor-1});
        if if_xlabel == 1
            xlabel(xtitle);
        end
        if if_ylabel == 1
            ylabel('Coefficients');
        end
        yl = ylim;
        % get significant point
        if isfield(input.left, 'pval') & ~isfield(input.left,'bootlow')  % if there is a pvalue and there is no bootstrap
            for ii = 1:length(pval(:,currPredictor,:))
                if pval(ii,currPredictor,:) < pvalThresh
                    hold on;
                    scatter(t(ii),1,25,'filled','black');
                end
            end
           
        else
            % this is average across session, obtained a CI with bootstrap
            % plot the shaded errors
            hold on;
            gray=[0.7 0.7 0.7];
            errorshade(t,input.left.bootlow(:,currPredictor),input.left.boothigh(:,currPredictor),gray);
            hold on;
            plot(t,coeff.left(:,currPredictor,:),'k.-','MarkerSize',15);
            
            hold on;
            gray=[0.7 0.7 0.7];
            errorshade(t,input.right.bootlow(:,currPredictor),input.right.boothigh(:,currPredictor),gray);
            hold on;
            plot(t,coeff.right(:,currPredictor,:),'k--','MarkerSize',15);
            
        end
        % plot the vertical line aligned to the cue
        
        % find the right limit
        ylim([lowerLimit upperLimit]);
        plot([0 0],[lowerLimit upperLimit],'k','LineWidth',1);
    end
end

if nInteraction > 0
    for l = 1:nInteraction
        for k=1:1+nback
            
            currPredictor=1+nPredictor*(1+nback)+(l-1)*(1+nback)+k;
            
            subplot(panelv,1+nback,currPredictor-1); hold on;
            
            plot(t,coeff.left(:,currPredictor,:),'k.-','MarkerSize',15);
            hold on;
            plot(t,coeff.right(:,currPredictor,:),'k--','MarkerSize',15);
            [if_xlabel,if_ylabel] = if_label(currPredictor, nPredictor, nback)
            xlim([t(1) t(end)]);
            title(tlabel{currPredictor-1});
            if if_xlabel == 1
                xlabel(xtitle);
            end
            if if_ylabel == 1
                ylabel('Coefficients');
            end
            yl = ylim;
            % get significant point
            if isfield(input, 'pval')
                for ii = 1:length(pval(:,currPredictor,:))
                    if pval(ii,currPredictor,:) < pvalThresh
                        hold on;
                        scatter(t(ii),1,25,'filled','black');
                    end
                end
              
            else
                hold on;
                errorshade(t,input.left.bootlow(:,currPredictor),input.left.boothigh(:,currPredictor),gray);
                hold on;
                plot(t,coeff.left(:,currPredictor,:),'k.-','MarkerSize',15);
                
                 hold on;
                gray=[0.7 0.7 0.7];
                errorshade(t,input.right.bootlow(:,currPredictor),input.right.boothigh(:,currPredictor),gray);
                hold on;
                plot(t,coeff.right(:,currPredictor,:),'k--','MarkerSize',15);
            end
             % plot the vertical line aligned to the cue
        	ylim([lowerLimit upperLimit]);
            plot([0 0],[lowerLimit upperLimit],'k','LineWidth',1);
        end
    end
end



end