function plot_session_game(stats,n_plot,tlabel)
% % plot_session_game %
%PURPOSE:   Plot performance of a single session of two-player game
%AUTHORS:   AC Kwan 180424
%
%INPUT ARGUMENTS
%   stats:  stats of the game
%   n_plot: plot up to this number of trials
%   tlabel: Text string that will be put on top of the figure

setup_figprop;

%%
figure;

% subplot(4,1,1); hold on;
% bar(stats.r(:,1),1,'k');  %reward
% ylabel({stats.playerlabel{1},'Outcome'},'interpreter','none');
% xlim([0 n_plot]);
% set(gca,'xticklabel',[]);
% set(gca,'ytick',[0 1]);
% set(gca,'yticklabel',{'','R'});
% legend(['Reward rate = ' num2str(nanmean(stats.r(:,1)))],'location','NorthEast');
% title(tlabel);
%bar(0.2*stats.opto(:,1),1,'b'); %opto stim

n_plot = ceil(length(stats.c)/10)*10;
subplot(4,1,1); hold on;
bar(-1*(stats.c(:,1)==-1),1,'FaceColor','r','EdgeColor','none');
bar((stats.c(:,1)==1),1,'FaceColor','b','EdgeColor','none');
hold on;
for ii = 1:length(stats.r)
    if stats.r(ii,1) == 1
        if stats.c(ii,1) == 1
            hold on; plot([ii,ii],[1.35,1.75],'black');  %reward right
        elseif stats.c(ii,1) == -1
            hold on; plot([ii,ii],[-1.35,-1.75],'black');
        end
    end
end

ylabel({stats.playerlabel{1},'Choice'},'interpreter','none');
%Lgd = legend(['Prob(L) = ' num2str(nanmean(stats.c(:,1)==-1))],'location',[0.8,0.95,0.05,0.05]);
%legend boxoff
%set(Lgd,'EdgeColor','none');
%set(Lgd,'color','none');
%txt = ['Reward rate = ' num2str(nanmean(stats.r(:,1)))];
xlim([0 n_plot]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[-1.75 -1 1 1.75]);
set(gca,'yticklabel',{'Reward','L','R','Reward'});
pos = get(gca, 'Position');
pos(1) = 0.2;
pos(4) = 0.16;
pos(3) = 0.7;
set(gca, 'Position', pos)

subplot(4,1,2); hold on;
bar(-1*(stats.c(:,2)==-1),1,'FaceColor','r','EdgeColor','none');
bar((stats.c(:,2)==1),1,'FaceColor','b','EdgeColor','none');
ylabel({'Computer', 'Choice'},'interpreter','none');
%Lgd = legend(['Prob(L) = ' num2str(nanmean(stats.c(:,2)==-1))],'location',[0.8,0.75,0.05,0.05]);
%set(Lgd,'EdgeColor','none');
xlim([0 n_plot]);
set(gca,'ytick',[-1 1]);
set(gca,'yticklabel',{'L','R'});

xlabel('Trial');

pos = get(gca, 'Position');
pos(1) = 0.2;
pos(4) = 0.08;
pos(3) = 0.7;
set(gca, 'Position', pos)

if isfield(stats,'pl')
    set(gca,'xticklabel',[]);  %if plotting 4th panel, then set 3rd panel's xtick label off

    subplot(4,1,4); hold on;
    plot(1-stats.pl(:,2),'b','LineWidth',2);
    ylabel({stats.playerlabel{2},'Prob(R)'},'interpreter','none');
    xlabel('Trial');
    ylim([0 1]);
    xlim([0 n_plot]);
    set(gca,'ytick',[0 1]);
end

print(gcf,'-dpng','session-beh');    %png format
saveas(gcf, 'session-beh', 'fig');
saveas(gcf, 'sesson-beh','svg');
end

