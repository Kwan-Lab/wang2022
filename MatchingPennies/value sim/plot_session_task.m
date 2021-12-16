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

subplot(4,1,1); hold on;
bar(stats.r(:,1),1,'k');
ylabel({stats.playerlabel{1},'Outcome'},'interpreter','none');
xlim([0 n_plot]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[0 1]);
set(gca,'yticklabel',{'','R'});
legend(['Reward rate = ' num2str(sum(stats.r(:,1)==1)/(sum(stats.r(:,1)==1)+sum(stats.r(:,1)==0)))]);
title(tlabel);

subplot(4,1,2); hold on;
bar(-1*(stats.c(:,1)==-1),1,'FaceColor','r','EdgeColor','none');
hold on;
bar((stats.c(:,1)==1),1,'FaceColor','b','EdgeColor','none');
ylabel({stats.playerlabel{1},'Choice'},'interpreter','none');
xlim([0 n_plot]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[-1 1]);
set(gca,'yticklabel',{'L','R'});

subplot(4,1,4); hold on;
plot(stats.rewardprob(:,1),'r','LineWidth',2);
hold on;
plot(stats.rewardprob(:,2),'b','LineWidth',2);
ylabel('Reward prob.');
legend('Left','Right');
xlim([0 n_plot]);
xlabel('Trial');

print(gcf,'-dpng','session-beh');    %png format
saveas(gcf, 'session-beh', 'fig');

end

