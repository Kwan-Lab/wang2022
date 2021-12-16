function plot_session_sim(stats,stats2,n_plot,tlabel,tlabel2)
% % plot_session_task %
%PURPOSE:   Plot performance of a single session, comparing mice with
%           simulation
%AUTHORS:   AC Kwan 191212
%
%INPUT ARGUMENTS
%   stats:  stats of the task
%   n_plot: plot up to this number of trials
%   tlabel: Text string that will be put on top of the figure

%%
figure;

subplot(3,1,1); hold on; % Reward probabilities
plot(stats.rewardprob(:,1),'r','LineWidth',2);
hold on;
plot(stats.rewardprob(:,2),'b','LineWidth',2);
ylabel('Reward probability (%)');
legend('Left','Right');
legend box off;
xlim([0 n_plot]);
ylim([0 1]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[0 0.1 0.7 1]);
set(gca,'yticklabel',{'','10','70', ''});
title(tlabel,'interpreter','none');

subplot(3,1,2); hold on; % Probability of choosing right
plot(movmean((stats.c==1),10),'b-','Linewidth',2); %experiment, moving average of 10 trials
plot(1-stats2.pl,'m-','LineWidth',2); %model
ylabel({'P(Right)'},'interpreter','none');
xlim([0 n_plot]);
ylim([0 1]);
xlabel('Trial');
title(tlabel2,'interpreter','none');
legend('Animal','Model');
legend box off;

subplot(3,1,3); hold on; % Latent variables
if isfield(stats2,'ql')
    plot(stats2.ql,'k-','LineWidth',2);
end
if isfield(stats2,'qr')
    plot(stats2.qr,'g-','LineWidth',2);
end
ylabel({'Latent variables'},'interpreter','none');
xlim([0 n_plot]);
xlabel('Trial');
title(tlabel2,'interpreter','none');
legend({'Q_L','Q_R'});
legend box off;

end

