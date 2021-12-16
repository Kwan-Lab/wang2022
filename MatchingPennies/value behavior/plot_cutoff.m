function plot_cutoff(trialData, saveBehFigPath, dataIndex)


setup_figprop

% get the dummy code for choice
c=double(trialData.response);
for i=1:length(c)
    if c(i)~=0
        c(i)=(c(i)-2.5)*2;
    end
end

n_plot=ceil(length(trialData.cue) / 100 * 100);
    
xAxis = 1:length(trialData.aveEntropy);
figure;
subplot(2,1,1);
plot(xAxis, trialData.runningEntropy, xAxis, trialData.aveEntropy); hold on;
if trialData.cutPoint ~= 0
    scatter(trialData.cutPoint-1, trialData.aveEntropy(trialData.cutPoint-1), 80,'black','filled');
end
ylabel('Running Entropy (30 trials)');
xlim([0 n_plot]);

subplot(2,1,2);

bar(-1*(c==-1),1,'FaceColor','r','EdgeColor','none');
hold on;
bar((c==1),1,'FaceColor','b','EdgeColor','none');
ylabel("Animal's choice");
yticks([-1 1]);
yticklabels({'Left','Right'});
xlim([0 n_plot]);
xlabel('Trials');

title('The cutoff point for fatigue');


% format long g
print(gcf, '-r0', [saveBehFigPath dataIndex.Animal{1} '_' num2str(dataIndex.DateNumber) '_fatigue' ], '-dpng'); %png format ], '-dpng'); %png format

close all

end
