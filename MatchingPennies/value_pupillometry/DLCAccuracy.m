% read in the human labelling and DLC labeling, calculate the distance
% between them

% plot: 1) distance against frame
%       2) mean and standard deviation of the distance 
setup_figprop;
savefigpath = 'J:\figures';
cd(savefigpath);
humanLabel_path = 'J:\figures\CollectedData_Hongli.csv';
dlcLabel_path = 'J:\figures\dlc_est.csv';
pupilData_path = 'J:\MatchingPennies\pupilData\data\870\M870_phase2_algo2_WithPupil_1911131005DeepCut_resnet50_PupillometryNov4shuffle1_350000.csv';
%humanCoord = loadPupil(pupilData_path);
[distance, xoff, yoff] = loadLabel(humanLabel_path,dlcLabel_path);

% plot the distance against frame
figure;
for ii = 1:5
    plot(distance(ii,:),'LineWidth',1.5); hold on;
end
% this plot would show that there are some outliers
xlabel('Frames');
ylabel('Distance between human labels and DLC labels (px)');
legend('Up', 'Down','Left','Right','Center');
print(gcf,'-dpng','distance-frame');    %png format
saveas(gcf, 'distance-frame', 'fig');
saveas(gcf, 'distance-frame','svg');

% plot the mean and standard deviation of the distance
meanDis = nanmean(distance,2);
stdDis = nanstd(distance');
name = {'Up';'Down';'Left';'Right';'Center'};
x = [1:5]; 
figure;
bar(x, meanDis,'black');
hold on; errorbar(x, meanDis, zeros(size(x)), stdDis, 'black','LineWidth',2);
set(gca,'xticklabel',name)
ylabel('x-offsets of DLC labels (px)','FontSize',28);
print(gcf,'-dpng','distance-stats');    %png format
saveas(gcf, 'distance-stats', 'fig');
saveas(gcf, 'distance-stats','svg');


for ii=1:5
    Ind = 1:size(xoff,2);
    value = xoff(ii,:);
    outx{ii} = {num2str(value(value>=30 | value <=-30))}
end
figure;
boxplot(xoff','Colors','k','Symbol','k+','Labels',{'Up','Down','Left','Right','Center'});
hold on;annotation ( 'textarrow', [0.52 0.52], [0.18 0.13],'String',outx{3})
hold on;annotation ( 'textarrow', [0.85 0.85], [0.8 0.85],'String',outx{5})
set(gca, 'box','off');
ylabel('Difference between human labels and DLC labels');
ylim([-30,30])
print(gcf,'-dpng','distance-xoffsets');    %png format
saveas(gcf, 'distance-xoffsets', 'fig');
saveas(gcf, 'distance-xoffsets','svg');

for ii=1:5
    Ind = 1:size(yoff,2);
    value = yoff(ii,:);
    outy{ii} = {num2str(value(value>=30 | value <=-30))}
end
figure;
boxplot(yoff','Colors','k','Symbol','k+','Labels',{'Up','Down','Left','Right','Center'});
hold on;annotation ( 'textarrow', [0.22 0.22], [0.8 0.85],'String',outy{1})
hold on;annotation ( 'textarrow', [0.52 0.52], [0.18 0.13],'String',outy{3})
hold on;annotation ( 'textarrow', [0.85 0.85], [0.8 0.85],'String',outy{5})
set(gca, 'box','off');
ylabel('Difference between human labels and DLC labels');
ylim([-30,30])
print(gcf,'-dpng','distance-yoffsets');    %png format
saveas(gcf, 'distance-yoffsets', 'fig');
saveas(gcf, 'distance-yoffsets','svg');

%% violin plot
x = {'Top','Bottom','Left','Right','Center'};
figure;
% plot them separately to set the colors
v=violinplot(xoff',x,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0],'ShowData',false);
hold on;
plot([0 6],[0 0], 'k--','LineWidth',2);
v(1).ViolinColor = [6 6 122]/255;
v(2).ViolinColor = [2 124 246]/255;
v(3).ViolinColor = [108 206 108]/255;
v(4).ViolinColor = [245 145 1]/255;
v(5).ViolinColor = [122 1 1]/255;

set(gca, 'box','off');
ylabel('x-offsets of DLC labels (px)','FontSize',28);
ylim([-30,30])
print(gcf,'-dpng','distance-xoffsets-violin');    %png format
saveas(gcf, 'distance-xoffsets-violin', 'fig');
saveas(gcf, 'distance-xoffsets-violin','svg');

figure;
v=violinplot(yoff',x,'ViolinAlpha',0.8,'EdgeColor',[0.5 0.5 0.5],'BoxColor',[0.5 0.5 0.5],'MedianColor',[1 0 0],'ShowData',false);
hold on;
plot([0 6],[0 0], 'k--','LineWidth',2);
v(1).ViolinColor = [6 6 122]/255;
v(2).ViolinColor = [2 124 246]/255;
v(3).ViolinColor = [108 206 108]/255;
v(4).ViolinColor = [245 145 1]/255;
v(5).ViolinColor = [122 1 1]/255;
set(gca, 'box','off');
ylabel('y-offsets of DLC labels','FontSize',28);
ylim([-30,30])
print(gcf,'-dpng','distance-yoffsets-violin');    %png format
saveas(gcf, 'distance-yoffsets-violin', 'fig');
saveas(gcf, 'distance-yoffsets-violin','svg');

