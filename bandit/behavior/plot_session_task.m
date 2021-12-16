function plot_session_task(stats,n_plot,tlabel)
% % plot_session_task %
%PURPOSE:   Plot performance of a single session of one-player task
%AUTHORS:   AC Kwan 180424
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

subplot(3,1,2); hold on; % Choices and outcomes
bar(-1*(stats.c(:,1)==-1 & stats.r(:,1)==1),1,'FaceColor','k','EdgeColor','none');
bar(-0.8*(stats.c(:,1)==-1 & stats.r(:,1)==1),1,'FaceColor','w','EdgeColor','none'); %use the white color to create space
bar(-0.7*(stats.c(:,1)==-1),1,'FaceColor','r','EdgeColor','none');
bar(1*(stats.c(:,1)==1 & stats.r(:,1)==1),1,'FaceColor','k','EdgeColor','none');
bar(0.8*(stats.c(:,1)==1 & stats.r(:,1)==1),1,'FaceColor','w','EdgeColor','none');
bar(0.7*(stats.c(:,1)==1),1,'FaceColor','b','EdgeColor','none');
ylabel({'Choice'},'interpreter','none');
xlim([0 n_plot]);
ylim([-1 1]);
xlabel('Trial');
set(gca,'ytick',[-1 -0.7 0.7 1]);
set(gca,'yticklabel',{'Reward','Left','Right','Reward'});
title(['Overall reward rate = ' num2str(sum(stats.r(:,1)==1)/(sum(stats.r(:,1)==1)+sum(stats.r(:,1)==0)))]);

end

