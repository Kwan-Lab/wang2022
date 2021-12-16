function plot_session_qparam(stats,x,n_plot)
% % plot_session_qparam %
%PURPOSE:   Plot the q-learning parameters across a behavioral session
%AUTHORS:   AC Kwan 170518
%
%INPUT ARGUMENTS
%   stats:  stats of the game
%   x:      which player? (1 or 2)
%   n_plot:     plot up to this number of trials

%%

figure;

for j=1:4
    switch j
        case 1
            val = stats.ql(:,x);
            col = 'r';
            ytitle = 'Q_L';
        case 2
            val = stats.qr(:,x);
            col = 'b';
            ytitle = 'Q_R';
        case 3
            val = stats.ql(:,x)-stats.qr(:,x);
            col = 'k';
            ytitle = 'Q_L - Q_R';
        case 4
            val = stats.rpe(:,x);
            col = 'k';
            ytitle = 'RPE';
    end
    
    subplot(4,1,j);
    plot(val,col,'LineWidth',2);
    ylabel(ytitle);
    ylim([floor(nanmin(val)) ceil(nanmax(val))]);
    xlim([0 n_plot]);
    set(gca,'xticklabel',[]);
    
    if j ==1
        title(['Player ' int2str(x) ': ' stats.playerlabel{x}],'interpreter','none');
    end
    if j < 4
        set(gca,'xticklabel',[]);
    else
        xlabel('Trials');
    end
end

print(gcf,'-dpng',['session-qparam-player' int2str(x)]);    %png format
saveas(gcf,['session-qparam-player' int2str(x)], 'fig');
saveas(gcf, ['session-qparam-player' int2str(x)],'svg');

% plot  the latent variable in the same plot along with choice and reward
% get actual probability to choose L? gaussian moving average filter
choice = stats.c;
choice(choice==1) = 0;
choice(choice==-1) = 1;

PL = smoothdata(choice,'gaussian',15);

% try binomial expansion to mimic a gaussian kernel
% h = [1/2 1/2];
% binomialCoeff = conv(h,h);
% for n = 1:10
%     binomialCoeff = conv(binomialCoeff,h);
% end
% 
% fDelay = (length(binomialCoeff)-1)/2;
% binomialMA = filter(binomialCoeff, 1, choice);

if isfield(stats ,'ckl')
    
%% the following commented codes were specifically for bandit example plot   
%     figure;
%     xAxis = 469:968;  % match the example session of the bandit task
%     subplot(3,1,2)
%     plot(xAxis-469,stats.ql(xAxis,x),'r');
%     hold on; plot(xAxis-469,stats.qr(xAxis,x),'b');
%     set(gca,'xticklabel',[]);
%     ylabel('Action values');
%     Lgd = legend('Q_L','Q_R');
%     set(Lgd,'EdgeColor','none');
%     ylim([0 1]);
%     set(gca,'xticklabel',[]);
%     set(gca,'box','off') 
%     
%     subplot(3,1,3)
%     plot(xAxis-469,stats.ckl(xAxis,x),'r');
%     hold on; plot(xAxis-469,stats.ckr(xAxis,x),'b');
%     Lgd = legend('C_L','C_R');
%     set(Lgd,'EdgeColor','none');
%     xlabel('Trials');
%     ylabel('Choice kernels');
%     set(gca,'xticklabel',[]);
%     set(gca,'box','off') 
%     
%     subplot(3,1,1)
%     plot(xAxis-469,stats.pl(xAxis,x),'Black');
%     ylabel('P_L');
%     hold on; plot(xAxis-469,PL(xAxis),'color',[0.7,0.7,0.7])
    %hold on; plot(xAxis, binomialMA(1:500));
    % plot the choice information
    
%     for ii = xAxis(1):xAxis(end)
%         if stats.c(ii,1) == 1 && stats.r(ii,1) == 1
%             hold on; plot([ii-468,ii-468],[1.05,1.25],'b');  %reward
%         elseif stats.c(ii,1) == 1 && stats.r(ii,1) == 0
%             hold on; plot([ii-468,ii-468],[1.05,1.15],'b');  %no reward
%         elseif stats.c(ii,1) == -1 && stats.r(ii,1) == 1
%              hold on; plot([ii-468,ii-468],[-0.05,-0.25],'r');  %reward
%         elseif stats.c(ii,1) == -1 && stats.r(ii,1) == 0
%              hold on; plot([ii-468,ii-468],[-0.05,-0.15],'r');  %no reward
%         end
%     end
%     %set(gca,'xticklabel',[]);
%     set(gca,'box','off') 
%     
    %% matching pennies plot
%     figure;
%     xAxis = 200:699;  
%     subplot(3,1,2)
%     plot(xAxis-199,stats.ql(xAxis,x),'r');
%     hold on; plot(xAxis-199,stats.qr(xAxis,x),'b');
%     set(gca,'xticklabel',[]);
%     ylabel('Action values');
%     Lgd = legend('Q_L','Q_R');
%     set(Lgd,'EdgeColor','none');
%     ylim([0 1]);
%     set(gca,'box','off') 
%     
%     subplot(3,1,3)
%     plot(xAxis-199,stats.ckl(xAxis,x),'r');
%     hold on; plot(xAxis-199,stats.ckr(xAxis,x),'b');
%     xlabel('Trials');
%     ylabel('Choice-autocorrelation');
%     Lgd = legend('C_L','C_R');
%     set(Lgd,'EdgeColor','none');
%     set(gca,'box','off') 
%     
%     subplot(3,1,1)
%     plot(xAxis-199,stats.pl(xAxis,x),'Black');
%     ylabel('P_L');
%     hold on; plot(xAxis-199,PL(xAxis),'color',[0.7,0.7,0.7])
%     %hold on; plot(xAxis, binomialMA(1:500));
%     % plot the choice information
% 
%     for ii = xAxis(1):xAxis(end)
%         if stats.c(ii,1) == 1 && stats.r(ii,1) == 1
%             hold on; plot([ii-199,ii-199],[1.05,1.25],'b');  %reward
%         elseif stats.c(ii,1) == 1 && stats.r(ii,1) == 0
%             hold on; plot([ii-199,ii-199],[1.05,1.15],'b');  %no reward
%         elseif stats.c(ii,1) == -1 && stats.r(ii,1) == 1
%              hold on; plot([ii-199,ii-199],[-0.05,-0.25],'r');  %reward
%         elseif stats.c(ii,1) == -1 && stats.r(ii,1) == 0
%              hold on; plot([ii-199,ii-199],[-0.05,-0.15],'r');  %no reward
%         end
%     end
%     set(gca,'xticklabel',[]);
%     set(gca,'box','off') 
%     
%     print(gcf,'-dpng',['session-qparam-player-choice' int2str(x)]);    %png format
%     saveas(gcf,['session-qparam-player-choice' int2str(x)], 'fig');
%     saveas(gcf, ['session-qparam-player-choice' int2str(x)],'svg');
end
%%




