function plot_session_PrL(stats,stats_sim,n_plot,tlabel)
% % plot_session_task %
%PURPOSE:   Plot performance of a single session
%           comparing measured and predicted probabilities of choosing left
%AUTHORS:   AC Kwan 180516
%
%INPUT ARGUMENTS
%   stats:  stats of the task
%   stats_sim:  stats of an agent simulated based on actual choices and outcomes
%   n_plot: plot up to this number of trials
%   tlabel: Text string that will be put on top of the figure

%%
figure;

c=nan(size(stats.c(:,1)));  %miss = NaN
c(stats.c(:,1)==1)=1;   %right = 1
c(stats.c(:,1)==-1)=0;  %left = 0

subplot(4,1,1); hold on;
plot(movmean(c,10),'k','LineWidth',2);
ylabel({stats.playerlabel{1},'Prob(R)'},'interpreter','none');
ylim([0 1]);
xlim([0 n_plot]);
set(gca,'ytick',[0 1]);
set(gca,'xticklabel',[]);

subplot(4,1,2); hold on;
plot(1-stats_sim.pl(:,1),'m','LineWidth',2);
ylabel({stats_sim.playerlabel{1},'Prob(R)'},'interpreter','none');
ylim([0 1]);
xlim([0 n_plot]);
set(gca,'ytick',[0 1]);
set(gca,'xticklabel',[]);

subplot(4,1,3); hold on;
plot(movmean(c,10),'k','LineWidth',2);
plot(1-stats_sim.pl(:,1),'m','LineWidth',2);
ylabel({'Overlay','Prob(R)'},'interpreter','none');
ylim([0 1]);
xlim([0 n_plot]);
set(gca,'ytick',[0 1]);

if isfield(stats,'rewardprob')
    set(gca,'xticklabel',[]);
    
    subplot(4,1,4); hold on;
    plot(stats.rewardprob(:,1),'r','LineWidth',2);
    hold on;
    plot(stats.rewardprob(:,2),'b','LineWidth',2);
    ylabel('Reward prob.');
    legend('Left','Right');
    xlim([0 n_plot]);
    xlabel('Trial');
end

print(gcf,'-dpng','session-behPrL');    %png format
saveas(gcf, 'session-behPrL', 'fig');

end

